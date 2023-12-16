/*
 * Copyright (c) 2019-2022 Genome Research Ltd.
 * Author(s): James Bonfield
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *    1. Redistributions of source code must retain the above copyright notice,
 *       this list of conditions and the following disclaimer.
 *
 *    2. Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *    3. Neither the names Genome Research Ltd and Wellcome Trust Sanger
 *       Institute nor the names of its contributors may be used to endorse
 *       or promote products derived from this software without specific
 *       prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH
 * LTD OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

// As per standard rANS_static but using optional RLE or bit-packing
// techniques prior to entropy encoding.  This is a significant
// reduction in some data sets.

// top bits in order byte
#define X_PACK   0x80    // Pack 2,4,8 or infinite symbols into a byte.
#define X_RLE    0x40    // Run length encoding with runs & lits encoded separately
#define X_CAT    0x20    // Nop; for tiny segments where rANS overhead is too big
#define X_NOSZ   0x10    // Don't store the original size; used by STRIPE mode
#define X_STRIPE 0x08    // For 4-byte integer data; rotate & encode 4 streams.
#define X_EXT    0x04    // External compression codec via magic num (gz, xz, bz2)
#define X_ORDER  0x03    // Mask to obtain order

#include "config.h"

#include <bzlib.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <assert.h>
#include <string.h>
#include <sys/time.h>
#include <limits.h>

#include "htscodecs/arith_dynamic.h"
#include "htscodecs/varint.h"
#include "htscodecs/pack.h"
#include "htscodecs/utils.h"

#define MIN(a,b) ((a)<(b)?(a):(b))

/*-----------------------------------------------------------------------------
 * Memory to memory compression functions.
 *
 * These are original versions without any manual loop unrolling. They
 * are easier to understand, but can be up to 2x slower.
 */
#define MAGIC 8

unsigned int arith_compress_bound(unsigned int size, int order) {
    int N = (order>>8) & 0xff;
    if (!N) N=4;
    return (order == 0
        ? 1.05*size + 257*3 + 4
        : 1.05*size + 257*257*3 + 4 + 257*3+4) + 5 +
        ((order & X_PACK) ? 1 : 0) +
        ((order & X_RLE) ? 1 + 257*3+4: 0) +
        ((order & X_STRIPE) ? 7 + 5*N: 0);
}

#ifndef MODEL_256 // see fqzcomp_qual_fuzz.c
#define NSYM 256
#include "htscodecs/c_simple_model.h"
#endif

static
unsigned char *arith_uncompress_O0(unsigned char *in, unsigned int in_size,
                                   unsigned char *out, unsigned int out_sz) {
    RangeCoder rc;
    int i;
    unsigned int m = in[0] ? in[0] : 256;

    SIMPLE_MODEL(256,_) byte_model;
    SIMPLE_MODEL(256,_init)(&byte_model, m);

    if (!out)
        out = malloc(out_sz);
    if (!out)
        return NULL;

    RC_SetInput(&rc, (char *)in+1, (char *)in+in_size);
    RC_StartDecode(&rc);

    for (i = 0; i < out_sz; i++)
        out[i] = SIMPLE_MODEL(256, _decodeSymbol)(&byte_model, &rc);

    RC_FinishDecode(&rc);
    
    return out;
}

static
unsigned char *arith_uncompress_O1(unsigned char *in, unsigned int in_size,
                                   unsigned char *out, unsigned int out_sz) {
    RangeCoder rc;
    unsigned char *out_free = NULL;

    if (!out)
        out_free = out = malloc(out_sz);
    if (!out)
        return NULL;


    SIMPLE_MODEL(256,_) *byte_model =
        htscodecs_tls_alloc(256 * sizeof(*byte_model));
    if (!byte_model) {
        free(out_free);
        return NULL;
    }

    unsigned int m = in[0] ? in[0] : 256, i;
    for (i = 0; i < 256; i++)
        SIMPLE_MODEL(256,_init)(&byte_model[i], m);

    RC_SetInput(&rc, (char *)in+1, (char *)in+in_size);
    RC_StartDecode(&rc);

    unsigned char last = 0;
    for (i = 0; i < out_sz; i++) {
        out[i] = SIMPLE_MODEL(256, _decodeSymbol)(&byte_model[last], &rc);
        last = out[i];
    }

    RC_FinishDecode(&rc);
    
    htscodecs_tls_free(byte_model);
    return out;
}

#undef NSYM
#define NSYM 258
#include "htscodecs/c_simple_model.h"
#define MAX_RUN 4

static
unsigned char *arith_uncompress_O0_RLE(unsigned char *in, unsigned int in_size,
                                       unsigned char *out, unsigned int out_sz) {
    RangeCoder rc;
    int i;
    unsigned int m = in[0] ? in[0] : 256;
    unsigned char *out_free = NULL;

    if (!out)
        out_free = out = malloc(out_sz);
    if (!out)
        return NULL;

    SIMPLE_MODEL(256,_) byte_model;
    SIMPLE_MODEL(256,_init)(&byte_model, m);

    SIMPLE_MODEL(NSYM,_) *run_model =
        htscodecs_tls_alloc(NSYM * sizeof(*run_model));
    if (!run_model) {
        free(out_free);
        return NULL;
    }

    for (i = 0; i < NSYM; i++)
        SIMPLE_MODEL(NSYM,_init)(&run_model[i], MAX_RUN);

    RC_SetInput(&rc, (char *)in+1, (char *)in+in_size);
    RC_StartDecode(&rc);

    for (i = 0; i < out_sz; i++) {
        unsigned char last;
        last = out[i] = SIMPLE_MODEL(256, _decodeSymbol)(&byte_model, &rc);
        //fprintf(stderr, "lit %c\n", last);
        int run = 0, r = 0, rctx = out[i];
        do {
            r = SIMPLE_MODEL(NSYM, _decodeSymbol)(&run_model[rctx], &rc);
            if (rctx == last)
                rctx = 256;
            else
                rctx += (rctx < NSYM-1);
            //fprintf(stderr, "run %d (ctx %d, %d)\n", r, last, l);
            run += r;
        } while (r == MAX_RUN-1 && run < out_sz);
        while (run-- && i+1 < out_sz)
            out[++i] = last;
    }

    RC_FinishDecode(&rc);

    htscodecs_tls_free(run_model);
    return out;
}
static
unsigned char *arith_uncompress_O1_RLE(unsigned char *in, unsigned int in_size,
                                       unsigned char *out, unsigned int out_sz) {
    RangeCoder rc;
    int i;
    unsigned int m = in[0] ? in[0] : 256;
    unsigned char *out_free = NULL;

    if (!out)
        out_free = out = malloc(out_sz);
    if (!out)
        return NULL;

    SIMPLE_MODEL(256,_) *byte_model =
        htscodecs_tls_alloc(256 * sizeof(*byte_model));
    if (!byte_model) {
        free(out_free);
        return NULL;
    }
    for (i = 0; i < 256; i++)
        SIMPLE_MODEL(256,_init)(&byte_model[i], m);

    SIMPLE_MODEL(NSYM,_) *run_model =
        htscodecs_tls_alloc(NSYM * sizeof(*run_model));
    if (!run_model) {
        htscodecs_tls_free(byte_model);
        free(out_free);
        return NULL;
    }
    for (i = 0; i < NSYM; i++)
        SIMPLE_MODEL(NSYM,_init)(&run_model[i], MAX_RUN);

    RC_SetInput(&rc, (char *)in+1, (char *)in+in_size);
    RC_StartDecode(&rc);

    unsigned char last = 0;
    for (i = 0; i < out_sz; i++) {
        out[i] = SIMPLE_MODEL(256, _decodeSymbol)(&byte_model[last], &rc);
        //fprintf(stderr, "lit %c (ctx %c)\n", out[i], last);
        last = out[i];
        int run = 0, r = 0, rctx = last;

        do {
            r = SIMPLE_MODEL(NSYM, _decodeSymbol)(&run_model[rctx], &rc);
            if (rctx == last)
                rctx = 256;
            else
                rctx += (rctx < NSYM-1);
            run += r;
        } while (r == MAX_RUN-1 && run < out_sz);
        while (run-- && i+1 < out_sz)
            out[++i] = last;
    }

    RC_FinishDecode(&rc);

    htscodecs_tls_free(byte_model);
    htscodecs_tls_free(run_model);
    return out;
}

unsigned char *arith_uncompress_to(unsigned char *in,  unsigned int in_size,
                                   unsigned char *out, unsigned int *out_size) {
    unsigned char *in_end = in + in_size;
    unsigned char *out_free = NULL;
    unsigned char *tmp_free = NULL;

    if (in_size == 0)
        return NULL;

    if (*in & X_STRIPE) {
        unsigned int ulen, olen, c_meta_len = 1;
        int i;
        uint64_t clen_tot = 0;

        // Decode lengths
        c_meta_len += var_get_u32(in+c_meta_len, in_end, &ulen);
        if (c_meta_len >= in_size)
            return NULL;
        unsigned int N = in[c_meta_len++];
        if (N < 1)  // Must be at least one stripe
            return NULL;
        unsigned int clenN[256], ulenN[256], idxN[256];
        if (!out) {
            if (ulen >= INT_MAX)
                return NULL;
            if (!(out_free = out = malloc(ulen))) {
                return NULL;
            }
            *out_size = ulen;
        }
        if (ulen != *out_size) {
            free(out_free);
            return NULL;
        }

        for (i = 0; i < N; i++) {
            ulenN[i] = ulen / N + ((ulen % N) > i);
            idxN[i] = i ? idxN[i-1] + ulenN[i-1] : 0;
            c_meta_len += var_get_u32(in+c_meta_len, in_end, &clenN[i]);
            clen_tot += clenN[i];
            if (c_meta_len > in_size || clenN[i] > in_size || clenN[i] < 1) {
                free(out_free);
                return NULL;
            }
        }

        // We can call this with a larger buffer, but once we've determined
        // how much we really use we limit it so the recursion becomes easier
        // to limit.
        if (c_meta_len + clen_tot > in_size) {
            free(out_free);
            return NULL;
        }
        in_size = c_meta_len + clen_tot;

        //fprintf(stderr, "    stripe meta %d\n", c_meta_len); //c-size

        // Uncompress the N streams
        unsigned char *outN = malloc(ulen);
        if (!outN) {
            free(out_free);
            return NULL;
        }
        for (i = 0; i < N; i++) {
            olen = ulenN[i];
            if (in_size < c_meta_len) {
                free(out_free);
                free(outN);
                return NULL;
            }
            if (!arith_uncompress_to(in+c_meta_len, in_size-c_meta_len, outN + idxN[i], &olen)
                || olen != ulenN[i]) {
                free(out_free);
                free(outN);
                return NULL;
            }
            c_meta_len += clenN[i];
        }

        unstripe(out, outN, ulen, N, idxN);

        free(outN);
        *out_size = ulen;
        return out;
    }

    int order = *in++;  in_size--;
    int do_pack = order & X_PACK;
    int do_rle  = order & X_RLE;
    int do_cat  = order & X_CAT;
    int no_size = order & X_NOSZ;
    int do_ext  = order & X_EXT;
    order &= 3;

    int sz = 0;
    unsigned int osz;
    if (!no_size)
        sz = var_get_u32(in, in_end, &osz);
    else
        sz = 0, osz = *out_size;
    in += sz;
    in_size -= sz;

    if (osz >= INT_MAX)
        return NULL;

#ifdef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
        // Limit maximum size to get fast turnaround on fuzzing test cases
        if (osz > 100000)
            goto err;
#endif

    if (no_size && !out)
        return NULL; // Need one or the other

    if (!out) {
        *out_size = osz;
        if (!(out_free = out = malloc(*out_size)))
            return NULL;
    } else {
        if (*out_size < osz)
            return NULL;
        *out_size = osz;
    }

    uint32_t c_meta_size = 0;
    unsigned int tmp1_size = *out_size;
    unsigned int tmp2_size = *out_size;
    unsigned char *tmp1 = NULL, *tmp2 = NULL, *tmp = NULL;

    // Need In, Out and Tmp buffers with temporary buffer of the same size
    // as output.  Our entropy decode is either arithmetic (with/without RLE)
    // or external (bz2, gzip, lzma) but with an optional unPACK transform
    // at the end.
    //
    // To avoid pointless memcpy when unpacking we switch around which
    // buffers we're writing to accordingly.

    // Format is pack meta data if present, followed by compressed data.
    if (do_pack) {
        if (!(tmp_free = tmp = malloc(*out_size)))
            goto err;
        tmp1 = tmp;  // uncompress
        tmp2 = out;  // unpack
    } else {
        // no pack
        tmp  = NULL;
        tmp1 = out;  // uncompress
        tmp2 = out;  // NOP
    }

    
    // Decode the bit-packing map.
    uint8_t map[16] = {0};
    int npacked_sym = 0;
    uint64_t unpacked_sz = 0; // FIXME: rename to packed_per_byte
    if (do_pack) {
        c_meta_size = hts_unpack_meta(in, in_size, *out_size, map, &npacked_sym);
        if (c_meta_size == 0)
            goto err;

        unpacked_sz = osz;
        in      += c_meta_size;
        in_size -= c_meta_size;

        // New unpacked size.  We could derive this bit from *out_size
        // and npacked_sym.
        unsigned int osz;
        sz = var_get_u32(in, in_end, &osz);
        in += sz;
        in_size -= sz;
        if (osz > tmp1_size)
            goto err;
        tmp1_size = osz;
    }

    //fprintf(stderr, "    meta_size %d bytes\n", (int)(in - orig_in)); //c-size

    // uncompress RLE data.  in -> tmp1
    if (in_size) {
#ifdef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
        // Limit maximum size to get fast turnaround on fuzzing test cases
        if (tmp1_size > 100000)
            goto err;
#endif
        if (do_cat) {
            //fprintf(stderr, "    CAT %d\n", tmp1_size); //c-size
            if (tmp1_size > in_size)
                goto err;
            if (tmp1_size > *out_size)
                goto err;
            memcpy(tmp1, in, tmp1_size);
        } else if (do_ext) {
            if (BZ_OK != BZ2_bzBuffToBuffDecompress((char *)tmp1, &tmp1_size,
                                                    (char *)in, in_size, 0, 0))
                goto err;
          } else {
            // in -> tmp1
            if (do_rle) {
                tmp1 = order == 1
                    ? arith_uncompress_O1_RLE(in, in_size, tmp1, tmp1_size)
                    : arith_uncompress_O0_RLE(in, in_size, tmp1, tmp1_size);
            } else {
                //if (order == 2)
                //    tmp1 = arith_uncompress_O2(in, in_size, tmp1, tmp1_size)
                //else
                tmp1 = order == 1
                    ? arith_uncompress_O1(in, in_size, tmp1, tmp1_size)
                    : arith_uncompress_O0(in, in_size, tmp1, tmp1_size);
            }
            if (!tmp1)
                goto err;
        }
    } else {
        tmp1_size = 0;
    }

    if (do_pack) {
        // Unpack bits via pack-map.  tmp1 -> tmp2
        if (npacked_sym == 1)
            unpacked_sz = tmp1_size;
        //uint8_t *porig = unpack(tmp2, tmp2_size, unpacked_sz, npacked_sym, map);
        //memcpy(tmp3, porig, unpacked_sz);
        if (!hts_unpack(tmp1, tmp1_size, tmp2, unpacked_sz, npacked_sym, map))
            goto err;
        tmp2_size = unpacked_sz;
    } else {
        tmp2_size = tmp1_size;
    }

    if (tmp)
        free(tmp);

    *out_size = tmp2_size;
    return tmp2;

 err:
    free(tmp_free);
    free(out_free);
    return NULL;
}

unsigned char *arith_uncompress(unsigned char *in, unsigned int in_size,
                                unsigned int *out_size) {
    return arith_uncompress_to(in, in_size, NULL, out_size);
}
