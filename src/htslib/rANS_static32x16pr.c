/*
 * Copyright (c) 2017-2023 Genome Research Ltd.
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

#include "config.h"

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <limits.h>

#include "htscodecs/rANS_word.h"
#include "htscodecs/rANS_static4x16.h"
#include "htscodecs/rANS_static16_int.h"
#include "htscodecs/varint.h"
#include "htscodecs/utils.h"

#define TF_SHIFT 12
#define TOTFREQ (1<<TF_SHIFT)


// 9-11 is considerably faster in the O1 variant due to reduced table size.
// We auto-tune between 10 and 12 though.  Anywhere from 9 to 14 are viable.
#ifndef TF_SHIFT_O1
#define TF_SHIFT_O1 12
#endif
#ifndef TF_SHIFT_O1_FAST
#define TF_SHIFT_O1_FAST 10
#endif
#define TOTFREQ_O1 (1<<TF_SHIFT_O1)
#define TOTFREQ_O1_FAST (1<<TF_SHIFT_O1_FAST)


#define NX 32

unsigned char *rans_uncompress_O0_32x16(unsigned char *in,
                                        unsigned int in_size,
                                        unsigned char *out,
                                        unsigned int out_sz) {
    if (in_size < 16) // 4-states at least
        return NULL;

    if (out_sz >= INT_MAX)
        return NULL; // protect against some overflow cases

#ifdef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
    if (out_sz > 100000)
        return NULL;
#endif

    /* Load in the static tables */
    unsigned char *cp = in, *out_free = NULL;
    unsigned char *cp_end = in + in_size;
    int i;
    uint32_t s3[TOTFREQ]; // For TF_SHIFT <= 12

    if (!out)
        out_free = out = malloc(out_sz);
    if (!out)
        return NULL;

    // Precompute reverse lookup of frequency.
    uint32_t F[256] = {0}, fsum;
    int fsz = decode_freq(cp, cp_end, F, &fsum);
    if (!fsz)
        goto err;
    cp += fsz;

    normalise_freq_shift(F, fsum, TOTFREQ);

    // Build symbols; fixme, do as part of decode, see the _d variant
    if (rans_F_to_s3(F, TF_SHIFT, s3))
        goto err;

    if (cp_end - cp < NX * 4)
        goto err;

    int z;
    RansState R[NX];
    for (z = 0; z < NX; z++) {
        RansDecInit(&R[z], &cp);
        if (R[z] < RANS_BYTE_L)
            goto err;
    }

    int out_end = (out_sz&~(NX-1));
    const uint32_t mask = (1u << TF_SHIFT)-1;
    cp_end -= NX*2; // worst case for renorm bytes

    // assume NX is divisible by 4
    assert(NX%4==0);

    // Unsafe loop with no ptr overflow checking within loop itself
    for (i=0; likely(i < out_end && cp < cp_end); i+=NX) {
        for (z = 0; z < NX; z+=4) {
            uint32_t S[4];
            S[0] = s3[R[z+0] & mask];
            S[1] = s3[R[z+1] & mask];
            S[2] = s3[R[z+2] & mask];
            S[3] = s3[R[z+3] & mask];

            R[z+0] = (S[0]>>(TF_SHIFT+8)) * (R[z+0] >> TF_SHIFT)
                + ((S[0]>>8) & mask);
            R[z+1] = (S[1]>>(TF_SHIFT+8)) * (R[z+1] >> TF_SHIFT)
                + ((S[1]>>8) & mask);
            R[z+2] = (S[2]>>(TF_SHIFT+8)) * (R[z+2] >> TF_SHIFT)
                + ((S[2]>>8) & mask);
            R[z+3] = (S[3]>>(TF_SHIFT+8)) * (R[z+3] >> TF_SHIFT)
                + ((S[3]>>8) & mask);

            out[i+z+0] = S[0];
            out[i+z+1] = S[1];
            out[i+z+2] = S[2];
            out[i+z+3] = S[3];

            RansDecRenorm(&R[z+0], &cp);
            RansDecRenorm(&R[z+1], &cp);
            RansDecRenorm(&R[z+2], &cp);
            RansDecRenorm(&R[z+3], &cp);

            if (NX%8==0) {
                z += 4;
                S[0] = s3[R[z+0] & mask];
                S[1] = s3[R[z+1] & mask];
                S[2] = s3[R[z+2] & mask];
                S[3] = s3[R[z+3] & mask];

                R[z+0] = (S[0]>>(TF_SHIFT+8)) * (R[z+0] >> TF_SHIFT)
                    + ((S[0]>>8) & mask);
                R[z+1] = (S[1]>>(TF_SHIFT+8)) * (R[z+1] >> TF_SHIFT)
                    + ((S[1]>>8) & mask);
                R[z+2] = (S[2]>>(TF_SHIFT+8)) * (R[z+2] >> TF_SHIFT)
                    + ((S[2]>>8) & mask);
                R[z+3] = (S[3]>>(TF_SHIFT+8)) * (R[z+3] >> TF_SHIFT)
                    + ((S[3]>>8) & mask);

                out[i+z+0] = S[0];
                out[i+z+1] = S[1];
                out[i+z+2] = S[2];
                out[i+z+3] = S[3];

                RansDecRenorm(&R[z+0], &cp);
                RansDecRenorm(&R[z+1], &cp);
                RansDecRenorm(&R[z+2], &cp);
                RansDecRenorm(&R[z+3], &cp);
            }
        }
    }

    // Safe loop
    for (; i < out_end; i+=NX) {
        for (z = 0; z < NX; z+=4) {
            uint32_t S[4];
            S[0] = s3[R[z+0] & mask];
            S[1] = s3[R[z+1] & mask];
            S[2] = s3[R[z+2] & mask];
            S[3] = s3[R[z+3] & mask];

            R[z+0] = (S[0]>>(TF_SHIFT+8)) * (R[z+0] >> TF_SHIFT)
                + ((S[0]>>8) & mask);
            R[z+1] = (S[1]>>(TF_SHIFT+8)) * (R[z+1] >> TF_SHIFT)
                + ((S[1]>>8) & mask);
            R[z+2] = (S[2]>>(TF_SHIFT+8)) * (R[z+2] >> TF_SHIFT)
                + ((S[2]>>8) & mask);
            R[z+3] = (S[3]>>(TF_SHIFT+8)) * (R[z+3] >> TF_SHIFT)
                + ((S[3]>>8) & mask);

            out[i+z+0] = S[0];
            out[i+z+1] = S[1];
            out[i+z+2] = S[2];
            out[i+z+3] = S[3];

            RansDecRenormSafe(&R[z+0], &cp, cp_end+NX*2);
            RansDecRenormSafe(&R[z+1], &cp, cp_end+NX*2);
            RansDecRenormSafe(&R[z+2], &cp, cp_end+NX*2);
            RansDecRenormSafe(&R[z+3], &cp, cp_end+NX*2);
        }
    }

    for (z = out_sz & (NX-1); z-- > 0; )
        out[out_end + z] = s3[R[z] & mask];

    //fprintf(stderr, "    0 Decoded %d bytes\n", (int)(cp-in)); //c-size

    return out;

 err:
    free(out_free);
    return NULL;
}

//#define MAGIC2 111
#define MAGIC2 179
//#define MAGIC2 0

unsigned char *rans_uncompress_O1_32x16(unsigned char *in,
                                        unsigned int in_size,
                                        unsigned char *out,
                                        unsigned int out_sz) {
    if (in_size < NX*4) // 4-states at least
        return NULL;

    if (out_sz >= INT_MAX)
        return NULL; // protect against some overflow cases

#ifdef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
    if (out_sz > 100000)
        return NULL;
#endif

    /* Load in the static tables */
    unsigned char *cp = in, *cp_end = in+in_size, *out_free = NULL;
    unsigned char *c_freq = NULL;
    int i;

    /*
     * Somewhat complex memory layout.
     * With shift==12 (TF_SHIFT_O1) we fill out use both sfb and fb.
     * With shift==10 (...O1_FAST)  we fill out and use s3 only.
     *
     * sfb+fb is larger, therefore we allocate this much memory.
     */
    uint8_t *sfb_ = htscodecs_tls_alloc(256*
                                        ((TOTFREQ_O1+MAGIC2)*sizeof(*sfb_)
                                         +256 * sizeof(fb_t)));
    if (!sfb_)
        return NULL;

    // sfb and fb are consecutive
    uint8_t *sfb[257];
    if ((*cp >> 4) == TF_SHIFT_O1) {
        for (i = 0; i <= 256; i++)
            sfb[i]=  sfb_ + i*(TOTFREQ_O1+MAGIC2);
    } else {
        for (i = 0; i <= 256; i++)
            sfb[i]=  sfb_ + i*(TOTFREQ_O1_FAST+MAGIC2);
    }
    fb_t (*fb)[256] = (fb_t (*)[256]) sfb[256];

    // NOTE: s3 overlaps sfb/fb
    uint32_t (*s3)[TOTFREQ_O1_FAST] = (uint32_t (*)[TOTFREQ_O1_FAST])sfb_;

    if (!out)
        out_free = out = malloc(out_sz);

    if (!out)
        goto err;

    //fprintf(stderr, "out_sz=%d\n", out_sz);

    // compressed header? If so uncompress it
    unsigned char *tab_end = NULL;
    unsigned char *c_freq_end = cp_end;
    unsigned int shift = *cp >> 4;
    if (*cp++ & 1) {
        uint32_t u_freq_sz, c_freq_sz;
        cp += var_get_u32(cp, cp_end, &u_freq_sz);
        cp += var_get_u32(cp, cp_end, &c_freq_sz);
        if (c_freq_sz > cp_end - cp)
            goto err;
        tab_end = cp + c_freq_sz;
        if (!(c_freq = rans_uncompress_O0_4x16(cp, c_freq_sz, NULL,u_freq_sz)))
            goto err;
        cp = c_freq;
        c_freq_end = c_freq + u_freq_sz;
    }

    // Decode order-0 symbol list; avoids needing in order-1 tables
    cp += decode_freq1(cp, c_freq_end, shift, NULL, s3, sfb, fb);

    if (tab_end)
        cp = tab_end;
    free(c_freq);
    c_freq = NULL;

    if (cp_end - cp < NX * 4)
        goto err;

    RansState R[NX];
    uint8_t *ptr = cp, *ptr_end = in + in_size - 2*NX;
    int z;
    for (z = 0; z < NX; z++) {
        RansDecInit(&R[z], &ptr);
        if (R[z] < RANS_BYTE_L)
            goto err;
    }

    int isz4 = out_sz/NX;
    int i4[NX], l[NX] = {0};
    for (z = 0; z < NX; z++)
        i4[z] = z*isz4;

    const int low_ent = in_size < 0.2 * out_sz;

    // Around 15% faster to specialise for 10/12 than to have one
    // loop with shift as a variable.
    if (shift == TF_SHIFT_O1) {
        // TF_SHIFT_O1 = 12
        const uint32_t mask = ((1u << TF_SHIFT_O1)-1);
        for (; likely(i4[0] < isz4);) {
            for (z = 0; z < NX; z+=4) {
                uint16_t m[4], c[4];

                c[0] = sfb[l[z+0]][m[0] = R[z+0] & mask];
                c[1] = sfb[l[z+1]][m[1] = R[z+1] & mask];
                c[2] = sfb[l[z+2]][m[2] = R[z+2] & mask];
                c[3] = sfb[l[z+3]][m[3] = R[z+3] & mask];

                R[z+0] = fb[l[z+0]][c[0]].f * (R[z+0]>>TF_SHIFT_O1);
                R[z+0] += m[0] - fb[l[z+0]][c[0]].b;

                R[z+1] = fb[l[z+1]][c[1]].f * (R[z+1]>>TF_SHIFT_O1);
                R[z+1] += m[1] - fb[l[z+1]][c[1]].b;

                R[z+2] = fb[l[z+2]][c[2]].f * (R[z+2]>>TF_SHIFT_O1);
                R[z+2] += m[2] - fb[l[z+2]][c[2]].b;

                R[z+3] = fb[l[z+3]][c[3]].f * (R[z+3]>>TF_SHIFT_O1);
                R[z+3] += m[3] - fb[l[z+3]][c[3]].b;

                out[i4[z+0]++] = l[z+0] = c[0];
                out[i4[z+1]++] = l[z+1] = c[1];
                out[i4[z+2]++] = l[z+2] = c[2];
                out[i4[z+3]++] = l[z+3] = c[3];

                if (!low_ent && likely(ptr < ptr_end)) {
                    RansDecRenorm(&R[z+0], &ptr);
                    RansDecRenorm(&R[z+1], &ptr);
                    RansDecRenorm(&R[z+2], &ptr);
                    RansDecRenorm(&R[z+3], &ptr);
                } else {
                    RansDecRenormSafe(&R[z+0], &ptr, ptr_end+2*NX);
                    RansDecRenormSafe(&R[z+1], &ptr, ptr_end+2*NX);
                    RansDecRenormSafe(&R[z+2], &ptr, ptr_end+2*NX);
                    RansDecRenormSafe(&R[z+3], &ptr, ptr_end+2*NX);
                }
            }
        }

        // Remainder
        for (; i4[NX-1] < out_sz; i4[NX-1]++) {
            uint32_t m = R[NX-1] & ((1u<<TF_SHIFT_O1)-1);
            unsigned char c = sfb[l[NX-1]][m];
            out[i4[NX-1]] = c;
            R[NX-1] = fb[l[NX-1]][c].f * (R[NX-1]>>TF_SHIFT_O1) +
                m - fb[l[NX-1]][c].b;
            RansDecRenormSafe(&R[NX-1], &ptr, ptr_end + 2*NX);
            l[NX-1] = c;
        }
    } else {
        // TF_SHIFT_O1 = 10
        const uint32_t mask = ((1u << TF_SHIFT_O1_FAST)-1);
        for (; likely(i4[0] < isz4);) {
            for (z = 0; z < NX; z+=4) {
                // Merged sfb and fb into single s3 lookup.
                // The m[4] array completely vanishes in this method.
                uint32_t S[4] = {
                    s3[l[z+0]][R[z+0] & mask],
                    s3[l[z+1]][R[z+1] & mask],
                    s3[l[z+2]][R[z+2] & mask],
                    s3[l[z+3]][R[z+3] & mask],
                };

                l[z+0] = out[i4[z+0]++] = S[0];
                l[z+1] = out[i4[z+1]++] = S[1];
                l[z+2] = out[i4[z+2]++] = S[2];
                l[z+3] = out[i4[z+3]++] = S[3];

                uint32_t F[4] = {
                    S[0]>>(TF_SHIFT_O1_FAST+8),
                    S[1]>>(TF_SHIFT_O1_FAST+8),
                    S[2]>>(TF_SHIFT_O1_FAST+8),
                    S[3]>>(TF_SHIFT_O1_FAST+8),
                };
                uint32_t B[4] = {
                    (S[0]>>8) & mask,
                    (S[1]>>8) & mask,
                    (S[2]>>8) & mask,
                    (S[3]>>8) & mask,
                };

                R[z+0] = F[0] * (R[z+0]>>TF_SHIFT_O1_FAST) + B[0];
                R[z+1] = F[1] * (R[z+1]>>TF_SHIFT_O1_FAST) + B[1];
                R[z+2] = F[2] * (R[z+2]>>TF_SHIFT_O1_FAST) + B[2];
                R[z+3] = F[3] * (R[z+3]>>TF_SHIFT_O1_FAST) + B[3];

                if (!low_ent && (ptr < ptr_end)) {
                    // branchless & asm
                    RansDecRenorm(&R[z+0], &ptr);
                    RansDecRenorm(&R[z+1], &ptr);
                    RansDecRenorm(&R[z+2], &ptr);
                    RansDecRenorm(&R[z+3], &ptr);
                } else {
                    // branched, but better when predictable
                    RansDecRenormSafe(&R[z+0], &ptr, ptr_end+2*NX);
                    RansDecRenormSafe(&R[z+1], &ptr, ptr_end+2*NX);
                    RansDecRenormSafe(&R[z+2], &ptr, ptr_end+2*NX);
                    RansDecRenormSafe(&R[z+3], &ptr, ptr_end+2*NX);
                }
            }
        }

        // Remainder
        for (; i4[NX-1] < out_sz; i4[NX-1]++) {
            uint32_t S = s3[l[NX-1]][R[NX-1] & ((1u<<TF_SHIFT_O1_FAST)-1)];
            out[i4[NX-1]] = l[NX-1] = S&0xff;
            R[NX-1] = (S>>(TF_SHIFT_O1_FAST+8)) * (R[NX-1]>>TF_SHIFT_O1_FAST)
                + ((S>>8) & ((1u<<TF_SHIFT_O1_FAST)-1));
            RansDecRenormSafe(&R[NX-1], &ptr, ptr_end + 2*NX);
        }
    }
    //fprintf(stderr, "    1 Decoded %d bytes\n", (int)(ptr-in)); //c-size

    htscodecs_tls_free(sfb_);
    return out;

 err:
    htscodecs_tls_free(sfb_);
    free(out_free);
    free(c_freq);

    return NULL;
}
