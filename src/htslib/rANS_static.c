/*
 * Copyright (c) 2014-2022 Genome Research Ltd.
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

// Use 11 for order-1?
#define TF_SHIFT 12
#define TOTFREQ (1<<TF_SHIFT)

#include "htscodecs/rANS_byte.h"
#include "htscodecs/utils.h"

/*-------------------------------------------------------------------------- */
/*
 * Example wrapper to use the rans_byte.h functions included above.
 *
 * This demonstrates how to use, and unroll, an order-0 and order-1 frequency
 * model.
 */

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <assert.h>
#include <string.h>
#include <limits.h>
#include <sys/time.h>
#ifndef NO_THREADS
#include <pthread.h>
#endif

#include "htscodecs/rANS_static.h"

#define ABS(a) ((a)>0?(a):-(a))

typedef struct {
    unsigned char R[TOTFREQ];
} ari_decoder;

static
unsigned char *rans_uncompress_O0(unsigned char *in, unsigned int in_size,
                                  unsigned int *out_size) {
    /* Load in the static tables */
    unsigned char *cp = in + 9;
    unsigned char *cp_end = in + in_size;
    const uint32_t mask = (1u << TF_SHIFT)-1;
    int i, j, rle;
    unsigned int x, y;
    unsigned int out_sz, in_sz;
    char *out_buf;
    RansState R[4];
    RansState m[4];
    uint16_t sfreq[TOTFREQ+32];
    uint16_t ssym [TOTFREQ+32]; // faster, but only needs uint8_t
    uint32_t sbase[TOTFREQ+16]; // faster, but only needs uint16_t

    if (in_size < 26) // Need at least this many bytes just to start
        return NULL;

    if (*in++ != 0) // Order-0 check
        return NULL;
    
    in_sz  = ((in[0])<<0) | ((in[1])<<8) | ((in[2])<<16) | (((uint32_t)in[3])<<24);
    out_sz = ((in[4])<<0) | ((in[5])<<8) | ((in[6])<<16) | (((uint32_t)in[7])<<24);
    if (in_sz != in_size-9)
        return NULL;

    if (out_sz >= INT_MAX)
        return NULL; // protect against some overflow cases

    // For speeding up the fuzzer only.
    // Small input can lead to large uncompressed data.
    // We reject this as it just slows things up instead of testing more code
    // paths (once we've verified a few times for large data).
#ifdef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
    if (out_sz > 100000)
        return NULL;
#endif

    out_buf = malloc(out_sz);
    if (!out_buf)
        return NULL;

    //fprintf(stderr, "out_sz=%d\n", out_sz);

    // Precompute reverse lookup of frequency.
    rle = x = y = 0;
    j = *cp++;
    do {
        int F, C;
        if (cp > cp_end - 16) goto cleanup; // Not enough input bytes left
        if ((F = *cp++) >= 128) {
            F &= ~128;
            F = ((F & 127) << 8) | *cp++;
        }
        C = x;

        if (x + F > TOTFREQ)
            goto cleanup;

        for (y = 0; y < F; y++) {
            ssym [y + C] = j;
            sfreq[y + C] = F;
            sbase[y + C] = y;
        }
        x += F;

        if (!rle && j+1 == *cp) {
            j = *cp++;
            rle = *cp++;
        } else if (rle) {
            rle--;
            j++;
            if (j > 255)
                goto cleanup;
        } else {
            j = *cp++;
        }
    } while(j);

    if (x < TOTFREQ-1 || x > TOTFREQ)
        goto cleanup;
    if (x != TOTFREQ) {
        // Protection against accessing uninitialised memory in the case
        // where SUM(freqs) == 4095 and not 4096.
        ssym [x] = ssym [x-1];
        sfreq[x] = sfreq[x-1];
        sbase[x] = sbase[x-1]+1;
    }

    // 16 bytes of cp here. Also why cp - 16 in above loop.
    if (cp > cp_end - 16) goto cleanup; // Not enough input bytes left

    RansDecInit(&R[0], &cp); if (R[0] < RANS_BYTE_L) goto cleanup;
    RansDecInit(&R[1], &cp); if (R[1] < RANS_BYTE_L) goto cleanup;
    RansDecInit(&R[2], &cp); if (R[2] < RANS_BYTE_L) goto cleanup;
    RansDecInit(&R[3], &cp); if (R[3] < RANS_BYTE_L) goto cleanup;

    int out_end = (out_sz&~3);
    cp_end -= 8; // within 8 for simplicity of loop below
    // 2 x likely() here harms gcc 7.5 by about 8% rate drop, but only in O2
    for (i=0; likely(i < out_end); i+=4) {
        //                              /curr code
        // gcc7  O2 513/497   562/556++ 556/547 ok
        // gcc7  O3 566/552   569/553   581/563+
        // gcc10 O2 544/538   563/547   541/537-?
        // gcc10 O3 531/519   546/530   575/546+
        // gcc11 O2 512/490   588/540   540/535 mid
        // gcc11 O3 482/471   553/541   549/535
        // gcc12 O2 533/526   544/534   539/535
        // gcc12 O3 548/533   502/497-- 553/527 ok
        // clang10  555/542   564/549   560/541
        // clang13  560/553   572/559   556/559
        m[0] = R[0] & mask;
        R[0] = sfreq[m[0]] * (R[0] >> TF_SHIFT) + sbase[m[0]];

        m[1] = R[1] & mask;
        R[1] = sfreq[m[1]] * (R[1] >> TF_SHIFT) + sbase[m[1]];

        m[2] = R[2] & mask;
        R[2] = sfreq[m[2]] * (R[2] >> TF_SHIFT) + sbase[m[2]];

        m[3] = R[3] & mask;
        R[3] = sfreq[m[3]] * (R[3] >> TF_SHIFT) + sbase[m[3]];

        // likely() here harms gcc12 -O3
        if (cp<cp_end) {
            RansDecRenorm2(&R[0], &R[1], &cp);
            RansDecRenorm2(&R[2], &R[3], &cp);
        } else {
            RansDecRenormSafe(&R[0], &cp, cp_end+8);
            RansDecRenormSafe(&R[1], &cp, cp_end+8);
            RansDecRenormSafe(&R[2], &cp, cp_end+8);
            RansDecRenormSafe(&R[3], &cp, cp_end+8);
        }

        out_buf[i+0] = ssym[m[0]];
        out_buf[i+1] = ssym[m[1]];
        out_buf[i+2] = ssym[m[2]];
        out_buf[i+3] = ssym[m[3]];
    }


    switch(out_sz&3) {
    case 3:
        out_buf[out_end + 2] = ssym[R[2] & mask];
    case 2:
        out_buf[out_end + 1] = ssym[R[1] & mask];
    case 1:
        out_buf[out_end] = ssym[R[0] & mask];
    default:
        break;
    }
    
    *out_size = out_sz;
    return (unsigned char *)out_buf;

 cleanup:
    free(out_buf);
    return NULL;
}

static
unsigned char *rans_uncompress_O1(unsigned char *in, unsigned int in_size,
                                  unsigned int *out_size) {
    /* Load in the static tables */
    unsigned char *cp = in + 9;
    unsigned char *ptr_end = in + in_size;
    int i, j = -999, rle_i, rle_j;
    unsigned int x;
    unsigned int out_sz, in_sz;
    char *out_buf = NULL;

    // Sanity checking
    if (in_size < 27) // Need at least this many bytes to start
        return NULL;

    if (*in++ != 1) // Order-1 check
        return NULL;

    in_sz  = ((in[0])<<0) | ((in[1])<<8) | ((in[2])<<16) | (((uint32_t)in[3])<<24);
    out_sz = ((in[4])<<0) | ((in[5])<<8) | ((in[6])<<16) | (((uint32_t)in[7])<<24);
    if (in_sz != in_size-9)
        return NULL;

    if (out_sz >= INT_MAX)
        return NULL; // protect against some overflow cases

    // For speeding up the fuzzer only.
    // Small input can lead to large uncompressed data.
    // We reject this as it just slows things up instead of testing more code
    // paths (once we've verified a few times for large data).
#ifdef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
    if (out_sz > 100000)
        return NULL;
#endif

    // Allocate decoding lookup tables
    RansDecSymbol32 (*syms)[256];
    uint8_t *mem = htscodecs_tls_calloc(256, sizeof(ari_decoder)
                                        + sizeof(*syms));
    if (!mem)
        return NULL;
    ari_decoder *const D = (ari_decoder *)mem;
    syms = (RansDecSymbol32 (*)[256])(mem + 256*sizeof(ari_decoder));
    int16_t map[256], map_i = 0;
    
    memset(map, -1, 256*sizeof(*map));

    if (!D) goto cleanup;
    /* These memsets prevent illegal memory access in syms due to
       broken compressed data.  As D is calloc'd, all illegal transitions
       will end up in either row or column 0 of syms. */
    memset(&syms[0], 0, sizeof(syms[0]));
    for (i = 0; i < 256; i++)
        memset(&syms[i][0], 0, sizeof(syms[0][0]));

    //fprintf(stderr, "out_sz=%d\n", out_sz);

    //i = *cp++;
    rle_i = 0;
    i = *cp++;
    do {
        // Map arbitrary a,b,c to 0,1,2 to improve cache locality.
        if (map[i] == -1)
            map[i] = map_i++;
        int m_i = map[i];

        rle_j = x = 0;
        j = *cp++;
        do {
            if (map[j] == -1)
                map[j] = map_i++;

            int F, C;
            if (cp > ptr_end - 16) goto cleanup; // Not enough input bytes left
            if ((F = *cp++) >= 128) {
                F &= ~128;
                F = ((F & 127) << 8) | *cp++;
            }
            C = x;

            //fprintf(stderr, "i=%d j=%d F=%d C=%d\n", i, j, F, C);

            if (unlikely(!F))
                F = TOTFREQ;

            RansDecSymbolInit32(&syms[m_i][j], C, F);

            /* Build reverse lookup table */
            //if (!D[i].R) D[i].R = (unsigned char *)malloc(TOTFREQ);
            if (x + F > TOTFREQ)
                goto cleanup;

            memset(&D[m_i].R[x], j, F);
            x += F;

            if (!rle_j && j+1 == *cp) {
                j = *cp++;
                rle_j = *cp++;
            } else if (rle_j) {
                rle_j--;
                j++;
                if (j > 255)
                    goto cleanup;
            } else {
                j = *cp++;
            }
        } while(j);

        if (x < TOTFREQ-1 || x > TOTFREQ)
            goto cleanup;
        if (x < TOTFREQ) // historically we fill 4095, not 4096
            D[i].R[x] = D[i].R[x-1];

        if (!rle_i && i+1 == *cp) {
            i = *cp++;
            rle_i = *cp++;
        } else if (rle_i) {
            rle_i--;
            i++;
            if (i > 255)
                goto cleanup;
        } else {
            i = *cp++;
        }
    } while (i);
    for (i = 0; i < 256; i++)
        if (map[i] == -1)
            map[i] = 0;

    RansState rans0, rans1, rans2, rans3;
    uint8_t *ptr = cp;
    if (cp > ptr_end - 16) goto cleanup; // Not enough input bytes left
    RansDecInit(&rans0, &ptr); if (rans0 < RANS_BYTE_L) goto cleanup;
    RansDecInit(&rans1, &ptr); if (rans1 < RANS_BYTE_L) goto cleanup;
    RansDecInit(&rans2, &ptr); if (rans2 < RANS_BYTE_L) goto cleanup;
    RansDecInit(&rans3, &ptr); if (rans3 < RANS_BYTE_L) goto cleanup;

    RansState R[4];
    R[0] = rans0;
    R[1] = rans1;
    R[2] = rans2;
    R[3] = rans3;

    unsigned int isz4 = out_sz>>2;
    uint32_t l0 = 0;
    uint32_t l1 = 0;
    uint32_t l2 = 0;
    uint32_t l3 = 0;
    
    unsigned int i4[] = {0*isz4, 1*isz4, 2*isz4, 3*isz4};

    /* Allocate output buffer */
    out_buf = malloc(out_sz);
    if (!out_buf) goto cleanup;

    uint8_t cc0 = D[map[l0]].R[R[0] & ((1u << TF_SHIFT)-1)];
    uint8_t cc1 = D[map[l1]].R[R[1] & ((1u << TF_SHIFT)-1)];
    uint8_t cc2 = D[map[l2]].R[R[2] & ((1u << TF_SHIFT)-1)];
    uint8_t cc3 = D[map[l3]].R[R[3] & ((1u << TF_SHIFT)-1)];

    ptr_end -= 8;
    for (; likely(i4[0] < isz4); i4[0]++, i4[1]++, i4[2]++, i4[3]++) {
        // seq4-head2: file q40b
        //          O3      O2
        // gcc7     296/291 290/260
        // gcc10    292/292 290/261
        // gcc11    293/293 290/265
        // gcc12    293/290 291/266
        // clang10  293/290 296/272
        // clang13  300/290 290/266
        out_buf[i4[0]] = cc0;
        out_buf[i4[1]] = cc1;
        out_buf[i4[2]] = cc2;
        out_buf[i4[3]] = cc3;

        RansDecSymbol32 s[4] = {
            syms[l0][cc0],
            syms[l1][cc1],
            syms[l2][cc2],
            syms[l3][cc3],
        };
        RansDecAdvanceStep(&R[0], s[0].start, s[0].freq, TF_SHIFT);
        RansDecAdvanceStep(&R[1], s[1].start, s[1].freq, TF_SHIFT);
        RansDecAdvanceStep(&R[2], s[2].start, s[2].freq, TF_SHIFT);
        RansDecAdvanceStep(&R[3], s[3].start, s[3].freq, TF_SHIFT);

        // Likely here helps speed of high-entropy data by 10-11%,
        // but harms low entropy-data speed by 3-4%.
        if ((ptr < ptr_end)) {
            RansDecRenorm2(&R[0], &R[1], &ptr);
            RansDecRenorm2(&R[2], &R[3], &ptr);
        } else {
            RansDecRenormSafe(&R[0], &ptr, ptr_end+8);
            RansDecRenormSafe(&R[1], &ptr, ptr_end+8);
            RansDecRenormSafe(&R[2], &ptr, ptr_end+8);
            RansDecRenormSafe(&R[3], &ptr, ptr_end+8);
        }

        l0 = map[cc0];
        l1 = map[cc1];
        l2 = map[cc2];
        l3 = map[cc3];

        cc0 = D[l0].R[R[0] & ((1u << TF_SHIFT)-1)];
        cc1 = D[l1].R[R[1] & ((1u << TF_SHIFT)-1)];
        cc2 = D[l2].R[R[2] & ((1u << TF_SHIFT)-1)];
        cc3 = D[l3].R[R[3] & ((1u << TF_SHIFT)-1)];
    }

    // Remainder
    for (; i4[3] < out_sz; i4[3]++) {
        unsigned char c3 = D[l3].R[RansDecGet(&R[3], TF_SHIFT)];
        out_buf[i4[3]] = c3;

        uint32_t m = R[3] & ((1u << TF_SHIFT)-1);
        R[3] = syms[l3][c3].freq * (R[3]>>TF_SHIFT) + m - syms[l3][c3].start;
        RansDecRenormSafe(&R[3], &ptr, ptr_end+8);
        l3 = map[c3];
    }
    
    *out_size = out_sz;

 cleanup:
    htscodecs_tls_free(D);

    return (unsigned char *)out_buf;
}


unsigned char *rans_uncompress(unsigned char *in, unsigned int in_size,
                               unsigned int *out_size) {
    /* Both rans_uncompress functions need to be able to read at least 9
       bytes. */
    if (in_size < 9)
        return NULL;
    return in[0]
        ? rans_uncompress_O1(in, in_size, out_size)
        : rans_uncompress_O0(in, in_size, out_size);
}
