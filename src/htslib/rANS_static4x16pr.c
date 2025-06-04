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

// FIXME Can we get decoder to return the compressed sized read, avoiding
// us needing to store it?  Yes we can.  See c-size comments.  If we added all these
// together we could get rans_uncompress_to_4x16 to return the number of bytes
// consumed, avoiding the calling code from needed to explicitly stored the size.
// However the effect on name tokeniser is to save 0.1 to 0.2% so not worth it.

/*-------------------------------------------------------------------------- */
/*
 * Example wrapper to use the rans_byte.h functions included above.
 *
 * This demonstrates how to use, and unroll, an order-0 and order-1 frequency
 * model.
 */

#include "config.h"

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <assert.h>
#include <string.h>
#include <sys/time.h>
#include <limits.h>
#include <math.h>

#ifndef NO_THREADS
#include <pthread.h>
#endif

#include "htscodecs/rANS_word.h"
#include "htscodecs/rANS_static4x16.h"
#include "htscodecs/rANS_static16_int.h"
#include "htscodecs/pack.h"
#include "htscodecs/rle.h"
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

unsigned char *rans_uncompress_O0_4x16(unsigned char *in, unsigned int in_size,
                                       unsigned char *out, unsigned int out_sz) {
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
    unsigned char *cp_end = in + in_size - 8; // within 8 => be extra safe
    int i, j;
    unsigned int x, y;
    uint16_t sfreq[TOTFREQ+32];
    uint16_t sbase[TOTFREQ+32]; // faster to use 32-bit on clang
    uint8_t  ssym [TOTFREQ+64]; // faster to use 16-bit on clang

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
    for (j = x = 0; j < 256; j++) {
        if (F[j]) {
            if (F[j] > TOTFREQ - x)
                goto err;
            for (y = 0; y < F[j]; y++) {
                ssym [y + x] = j;
                sfreq[y + x] = F[j];
                sbase[y + x] = y;
            }
            x += F[j];
        }
    }

    if (x != TOTFREQ)
        goto err;

    if (cp+16 > cp_end+8)
        goto err;

    RansState R[4];
    RansDecInit(&R[0], &cp); if (R[0] < RANS_BYTE_L) goto err;
    RansDecInit(&R[1], &cp); if (R[1] < RANS_BYTE_L) goto err;
    RansDecInit(&R[2], &cp); if (R[2] < RANS_BYTE_L) goto err;
    RansDecInit(&R[3], &cp); if (R[3] < RANS_BYTE_L) goto err;

// Simple version is comparable to below, but only with -O3
//
//    for (i = 0; cp < cp_end-8 && i < (out_sz&~7); i+=8) {
//        for(j=0; j<8;j++) {
//          RansState m = RansDecGet(&R[j%4], TF_SHIFT);
//          R[j%4] = sfreq[m] * (R[j%4] >> TF_SHIFT) + sbase[m];
//          out[i+j] = ssym[m];
//          RansDecRenorm(&R[j%4], &cp);
//        }
//    }

    for (i = 0; cp < cp_end-8 && i < (out_sz&~7); i+=8) {
        for (j = 0; j < 8; j+=4) {
            RansState m0 = RansDecGet(&R[0], TF_SHIFT);
            RansState m1 = RansDecGet(&R[1], TF_SHIFT);
            out[i+j+0] = ssym[m0];
            out[i+j+1] = ssym[m1];

            R[0] = sfreq[m0] * (R[0] >> TF_SHIFT) + sbase[m0];
            R[1] = sfreq[m1] * (R[1] >> TF_SHIFT) + sbase[m1];

            RansState m2 = RansDecGet(&R[2], TF_SHIFT);
            RansState m3 = RansDecGet(&R[3], TF_SHIFT);

            RansDecRenorm(&R[0], &cp);
            RansDecRenorm(&R[1], &cp);

            R[2] = sfreq[m2] * (R[2] >> TF_SHIFT) + sbase[m2];
            R[3] = sfreq[m3] * (R[3] >> TF_SHIFT) + sbase[m3];

            RansDecRenorm(&R[2], &cp);
            RansDecRenorm(&R[3], &cp);

            out[i+j+2] = ssym[m2];
            out[i+j+3] = ssym[m3];
        }
    }

    // remainder
    for (; i < out_sz; i++) {
        RansState m = RansDecGet(&R[i%4], TF_SHIFT);
        R[i%4] = sfreq[m] * (R[i%4] >> TF_SHIFT) + sbase[m];
        out[i] = ssym[m];
        RansDecRenormSafe(&R[i%4], &cp, cp_end+8);
    }

    //fprintf(stderr, "    0 Decoded %d bytes\n", (int)(cp-in)); //c-size

    return out;

 err:
    free(out_free);
    return NULL;
}

//#define MAGIC2 111
#define MAGIC2 179
//#define MAGIC2 0

static
unsigned char *rans_uncompress_O1_4x16(unsigned char *in, unsigned int in_size,
                                       unsigned char *out, unsigned int out_sz) {
    if (in_size < 16) // 4-states at least
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
    int i, j = -999;
    unsigned int x;

    uint8_t *sfb_ = htscodecs_tls_alloc(256*(TOTFREQ_O1+MAGIC2)*sizeof(*sfb_));
    uint32_t (*s3)[TOTFREQ_O1_FAST] = (uint32_t (*)[TOTFREQ_O1_FAST])sfb_;
    // reuse the same memory for the fast mode lookup, but this only works
    // if we're on e.g. 12-bit freqs vs 10-bit freqs as needs 4x larger array.
    //uint32_t s3[256][TOTFREQ_O1_FAST];

    if (!sfb_)
        return NULL;
    fb_t (*fb)[256] = htscodecs_tls_alloc(256 * sizeof(*fb));
    if (!fb)
        goto err;
    uint8_t *sfb[256];
    if ((*cp >> 4) == TF_SHIFT_O1) {
        for (i = 0; i < 256; i++)
            sfb[i]=  sfb_ + i*(TOTFREQ_O1+MAGIC2);
    } else {
        for (i = 0; i < 256; i++)
            sfb[i]=  sfb_ + i*(TOTFREQ_O1_FAST+MAGIC2);
    }

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
        if (!(c_freq = rans_uncompress_O0_4x16(cp, c_freq_sz, NULL, u_freq_sz)))
            goto err;
        cp = c_freq;
        c_freq_end = c_freq + u_freq_sz;
    }

    // Decode order-0 symbol list; avoids needing in order-1 tables
    uint32_t F0[256] = {0};
    int fsz = decode_alphabet(cp, c_freq_end, F0);
    if (!fsz)
        goto err;
    cp += fsz;

    if (cp >= c_freq_end)
        goto err;

    const int s3_fast_on = in_size >= 100000;

    for (i = 0; i < 256; i++) {
        if (F0[i] == 0)
            continue;

        uint32_t F[256] = {0}, T = 0;
        fsz = decode_freq_d(cp, c_freq_end, F0, F, &T);
        if (!fsz)
            goto err;
        cp += fsz;

        if (!T) {
            //fprintf(stderr, "No freq for F_%d\n", i);
            continue;
        }

        normalise_freq_shift(F, T, 1<<shift);

        // Build symbols; fixme, do as part of decode, see the _d variant
        for (j = x = 0; j < 256; j++) {
            if (F[j]) {
                if (F[j] > (1<<shift) - x)
                    goto err;

                if (shift == TF_SHIFT_O1_FAST && s3_fast_on) {
                    int y;
                    for (y = 0; y < F[j]; y++)
                        s3[i][y+x] = (((uint32_t)F[j])<<(shift+8)) |(y<<8) |j;
                } else {
                    memset(&sfb[i][x], j, F[j]);
                    fb[i][j].f = F[j];
                    fb[i][j].b = x;
                }
                x += F[j];
            }
        }
        if (x != (1<<shift))
            goto err;
    }

    if (tab_end)
        cp = tab_end;
    free(c_freq);
    c_freq = NULL;

    if (cp+16 > cp_end)
        goto err;

    RansState rans0, rans1, rans2, rans3;
    uint8_t *ptr = cp, *ptr_end = in + in_size - 8;
    RansDecInit(&rans0, &ptr); if (rans0 < RANS_BYTE_L) goto err;
    RansDecInit(&rans1, &ptr); if (rans1 < RANS_BYTE_L) goto err;
    RansDecInit(&rans2, &ptr); if (rans2 < RANS_BYTE_L) goto err;
    RansDecInit(&rans3, &ptr); if (rans3 < RANS_BYTE_L) goto err;

    unsigned int isz4 = out_sz>>2;
    int l0 = 0, l1 = 0, l2 = 0, l3 = 0;
    unsigned int i4[] = {0*isz4, 1*isz4, 2*isz4, 3*isz4};

    RansState R[4];
    R[0] = rans0;
    R[1] = rans1;
    R[2] = rans2;
    R[3] = rans3;

    // Around 15% faster to specialise for 10/12 than to have one
    // loop with shift as a variable.
    if (shift == TF_SHIFT_O1) {
        // TF_SHIFT_O1 = 12

        const uint32_t mask = ((1u << TF_SHIFT_O1)-1);
        for (; i4[0] < isz4; i4[0]++, i4[1]++, i4[2]++, i4[3]++) {
            uint16_t m, c;
            c = sfb[l0][m = R[0] & mask];
            R[0] = fb[l0][c].f * (R[0]>>TF_SHIFT_O1) + m - fb[l0][c].b;
            out[i4[0]] = l0 = c;

            c = sfb[l1][m = R[1] & mask];
            R[1] = fb[l1][c].f * (R[1]>>TF_SHIFT_O1) + m - fb[l1][c].b;
            out[i4[1]] = l1 = c;

            c = sfb[l2][m = R[2] & mask];
            R[2] = fb[l2][c].f * (R[2]>>TF_SHIFT_O1) + m - fb[l2][c].b;
            out[i4[2]] = l2 = c;

            c = sfb[l3][m = R[3] & mask];
            R[3] = fb[l3][c].f * (R[3]>>TF_SHIFT_O1) + m - fb[l3][c].b;
            out[i4[3]] = l3 = c;

            if (ptr < ptr_end) {
                RansDecRenorm(&R[0], &ptr);
                RansDecRenorm(&R[1], &ptr);
                RansDecRenorm(&R[2], &ptr);
                RansDecRenorm(&R[3], &ptr);
            } else {
                RansDecRenormSafe(&R[0], &ptr, ptr_end+8);
                RansDecRenormSafe(&R[1], &ptr, ptr_end+8);
                RansDecRenormSafe(&R[2], &ptr, ptr_end+8);
                RansDecRenormSafe(&R[3], &ptr, ptr_end+8);
            }
        }

        // Remainder
        for (; i4[3] < out_sz; i4[3]++) {
            uint32_t m3 = R[3] & ((1u<<TF_SHIFT_O1)-1);
            unsigned char c3 = sfb[l3][m3];
            out[i4[3]] = c3;
            R[3] = fb[l3][c3].f * (R[3]>>TF_SHIFT_O1) + m3 - fb[l3][c3].b;
            RansDecRenormSafe(&R[3], &ptr, ptr_end + 8);
            l3 = c3;
        }
    } else if (!s3_fast_on) {
        // TF_SHIFT_O1 = 10 with sfb[256][1024] & fb[256]256] array lookup
        // Slightly faster for -o193 on q4 (high comp), but also less
        // initialisation cost for smaller data
        const uint32_t mask = ((1u << TF_SHIFT_O1_FAST)-1);
        for (; i4[0] < isz4; i4[0]++, i4[1]++, i4[2]++, i4[3]++) {
            uint16_t m, c;
            c = sfb[l0][m = R[0] & mask];
            R[0] = fb[l0][c].f * (R[0]>>TF_SHIFT_O1_FAST) + m - fb[l0][c].b;
            out[i4[0]] = l0 = c;

            c = sfb[l1][m = R[1] & mask];
            R[1] = fb[l1][c].f * (R[1]>>TF_SHIFT_O1_FAST) + m - fb[l1][c].b;
            out[i4[1]] = l1 = c;

            c = sfb[l2][m = R[2] & mask];
            R[2] = fb[l2][c].f * (R[2]>>TF_SHIFT_O1_FAST) + m - fb[l2][c].b;
            out[i4[2]] = l2 = c;

            c = sfb[l3][m = R[3] & mask];
            R[3] = fb[l3][c].f * (R[3]>>TF_SHIFT_O1_FAST) + m - fb[l3][c].b;
            out[i4[3]] = l3 = c;

            if (ptr < ptr_end) {
                RansDecRenorm(&R[0], &ptr);
                RansDecRenorm(&R[1], &ptr);
                RansDecRenorm(&R[2], &ptr);
                RansDecRenorm(&R[3], &ptr);
            } else {
                RansDecRenormSafe(&R[0], &ptr, ptr_end+8);
                RansDecRenormSafe(&R[1], &ptr, ptr_end+8);
                RansDecRenormSafe(&R[2], &ptr, ptr_end+8);
                RansDecRenormSafe(&R[3], &ptr, ptr_end+8);
            }
        }

        // Remainder
        for (; i4[3] < out_sz; i4[3]++) {
            uint32_t m3 = R[3] & ((1u<<TF_SHIFT_O1_FAST)-1);
            unsigned char c3 = sfb[l3][m3];
            out[i4[3]] = c3;
            R[3] = fb[l3][c3].f * (R[3]>>TF_SHIFT_O1_FAST) + m3 - fb[l3][c3].b;
            RansDecRenormSafe(&R[3], &ptr, ptr_end + 8);
            l3 = c3;
        }
    } else {
        // TF_SHIFT_O1_FAST.
        // Significantly faster for -o1 on q40 (low comp).
        // Higher initialisation cost, so only use if big blocks.
        const uint32_t mask = ((1u << TF_SHIFT_O1_FAST)-1);
        for (; i4[0] < isz4; i4[0]++, i4[1]++, i4[2]++, i4[3]++) {
            uint32_t S0 = s3[l0][R[0] & mask];
            uint32_t S1 = s3[l1][R[1] & mask];
            l0 = out[i4[0]] = S0;
            l1 = out[i4[1]] = S1;
            uint16_t F0 = S0>>(TF_SHIFT_O1_FAST+8);
            uint16_t F1 = S1>>(TF_SHIFT_O1_FAST+8);
            uint16_t B0 = (S0>>8) & mask;
            uint16_t B1 = (S1>>8) & mask;

            R[0] = F0 * (R[0]>>TF_SHIFT_O1_FAST) + B0;
            R[1] = F1 * (R[1]>>TF_SHIFT_O1_FAST) + B1;

            uint32_t S2 = s3[l2][R[2] & mask];
            uint32_t S3 = s3[l3][R[3] & mask];
            l2 = out[i4[2]] = S2;
            l3 = out[i4[3]] = S3;
            uint16_t F2 = S2>>(TF_SHIFT_O1_FAST+8);
            uint16_t F3 = S3>>(TF_SHIFT_O1_FAST+8);
            uint16_t B2 = (S2>>8) & mask;
            uint16_t B3 = (S3>>8) & mask;

            R[2] = F2 * (R[2]>>TF_SHIFT_O1_FAST) + B2;
            R[3] = F3 * (R[3]>>TF_SHIFT_O1_FAST) + B3;

            if (ptr < ptr_end) {
                RansDecRenorm(&R[0], &ptr);
                RansDecRenorm(&R[1], &ptr);
                RansDecRenorm(&R[2], &ptr);
                RansDecRenorm(&R[3], &ptr);
            } else {
                RansDecRenormSafe(&R[0], &ptr, ptr_end+8);
                RansDecRenormSafe(&R[1], &ptr, ptr_end+8);
                RansDecRenormSafe(&R[2], &ptr, ptr_end+8);
                RansDecRenormSafe(&R[3], &ptr, ptr_end+8);
            }
        }

        // Remainder
        for (; i4[3] < out_sz; i4[3]++) {
            uint32_t S = s3[l3][R[3] & ((1u<<TF_SHIFT_O1_FAST)-1)];
            l3 = out[i4[3]] = S;
            R[3] = (S>>(TF_SHIFT_O1_FAST+8)) * (R[3]>>TF_SHIFT_O1_FAST)
                + ((S>>8) & ((1u<<TF_SHIFT_O1_FAST)-1));
            RansDecRenormSafe(&R[3], &ptr, ptr_end + 8);
        }
    }
    //fprintf(stderr, "    1 Decoded %d bytes\n", (int)(ptr-in)); //c-size

    htscodecs_tls_free(fb);
    htscodecs_tls_free(sfb_);
    return out;

 err:
    htscodecs_tls_free(fb);
    htscodecs_tls_free(sfb_);
    free(out_free);
    free(c_freq);

    return NULL;
}

/*-----------------------------------------------------------------------------
 * r32x16 implementation, included here for now for simplicity
 */
#include "htscodecs/rANS_static32x16pr.h"

// Test interface for restricting the auto-detection methods so we
// can forcibly compare different implementations on the same machine.
// See RANS_CPU_ defines in rANS_static4x16.h
static int rans_cpu = 0xFFFF; // all
void rans_set_cpu(int opts) {
    rans_cpu = opts;
}

#if (defined(__GNUC__) || defined(__clang__)) && defined(__x86_64__)
// Icc and Clang both also set __GNUC__ on linux, but not on Windows.
#include <cpuid.h>

#if defined(__clang__) && defined(__has_attribute)
#  if __has_attribute(unused)
#    define UNUSED __attribute__((unused))
#  else
#    define UNUSED
#  endif
#elif defined(__GNUC__) && __GNUC__ >= 3
#  define UNUSED __attribute__((unused))
#else
#  define UNUSED
#endif

// CPU detection is performed once.  NB this has an assumption that we're
// not migrating between processes with different instruction stes, but
// to date the only systems I know of that support this don't have different
// capabilities (that we use) per core.
#ifndef NO_THREADS
static pthread_once_t rans_cpu_once = PTHREAD_ONCE_INIT;
#endif

static int have_ssse3   UNUSED = 0;
static int have_sse4_1  UNUSED = 0;
static int have_popcnt  UNUSED = 0;
static int have_avx2    UNUSED = 0;
static int have_avx512f UNUSED = 0;
static int is_amd       UNUSED = 0;

static void htscodecs_tls_cpu_init(void) {
    unsigned int eax = 0, ebx = 0, ecx = 0, edx = 0;
    // These may be unused, depending on HAVE_* config.h macros

    int level = __get_cpuid_max(0, NULL);
    __cpuid_count(0, 0, eax, ebx, ecx, edx);
    is_amd = (ecx == 0x444d4163);
    if (level >= 1) {
        __cpuid_count(1, 0, eax, ebx, ecx, edx);
#if defined(bit_SSSE3)
        have_ssse3 = ecx & bit_SSSE3;
#endif
#if defined(bit_POPCNT)
        have_popcnt = ecx & bit_POPCNT;
#endif
#if defined(bit_SSE4_1)
        have_sse4_1 = ecx & bit_SSE4_1;
#endif
    }
    if (level >= 7) {
        __cpuid_count(7, 0, eax, ebx, ecx, edx);
#if defined(bit_AVX2)
        have_avx2 = ebx & bit_AVX2;
#endif
#if defined(bit_AVX512F)
        have_avx512f = ebx & bit_AVX512F;
#endif
    }

    if (!have_popcnt) have_avx512f = have_avx2 = have_sse4_1 = 0;
    if (!have_ssse3)  have_sse4_1 = 0;

    if (!(rans_cpu & RANS_CPU_ENC_AVX512)) have_avx512f = 0;
    if (!(rans_cpu & RANS_CPU_ENC_AVX2))   have_avx2 = 0;
    if (!(rans_cpu & RANS_CPU_ENC_SSE4))   have_sse4_1 = 0;
}

static inline
unsigned char *(*rans_dec_func(int do_simd, int order))
    (unsigned char *in,
     unsigned int in_size,
     unsigned char *out,
     unsigned int out_size) {

    if (!do_simd) { // SIMD disabled
        return order & 1
            ? rans_uncompress_O1_4x16
            : rans_uncompress_O0_4x16;
    }

#ifdef NO_THREADS
    htscodecs_tls_cpu_init();
#else
    int err = pthread_once(&rans_cpu_once, htscodecs_tls_cpu_init);
    if (err != 0) {
        // commented to comply with cran's warning
        //fprintf(stderr, "Initialising TLS data failed: pthread_once: %s\n",
                //strerror(err));
        //fprintf(stderr, "Using scalar code only\n");
    }
#endif

    if (order & 1) {
#if defined(HAVE_AVX512)
        if (have_avx512f)
            return rans_uncompress_O1_32x16_avx512;
#endif
#if defined(HAVE_AVX2)
        if (have_avx2)
            return rans_uncompress_O1_32x16_avx2;
#endif
#if defined(HAVE_SSE4_1) && defined(HAVE_SSSE3) && defined(HAVE_POPCNT)
        if (have_sse4_1)
            return rans_uncompress_O1_32x16_sse4;
#endif
        return rans_uncompress_O1_32x16;
    } else {
#if defined(HAVE_AVX512)
        if (have_avx512f && (!is_amd || !have_avx2))
            return rans_uncompress_O0_32x16_avx512;
#endif
#if defined(HAVE_AVX2)
        if (have_avx2)
            return rans_uncompress_O0_32x16_avx2;
#endif
#if defined(HAVE_SSE4_1) && defined(HAVE_SSSE3) && defined(HAVE_POPCNT)
        if (have_sse4_1)
            return rans_uncompress_O0_32x16_sse4;
#endif
        return rans_uncompress_O0_32x16;
    }
}

#elif defined(__ARM_NEON) && defined(__aarch64__)

#if defined(__linux__) || defined(__FreeBSD__)
#include <sys/auxv.h>
#elif defined(_WIN32)
#include <processthreadsapi.h>
#endif

static inline int have_neon() {
#if defined(__linux__) && defined(__arm__)
    return (getauxval(AT_HWCAP) & HWCAP_NEON) != 0;
#elif defined(__linux__) && defined(__aarch64__)
    return (getauxval(AT_HWCAP) & HWCAP_ASIMD) != 0;
#elif defined(__APPLE__)
    return 1;
#elif defined(__FreeBSD__) && defined(__arm__)
    u_long cap;
    if (elf_aux_info(AT_HWCAP, &cap, sizeof cap) != 0) return 0;
    return (cap & HWCAP_NEON) != 0;
#elif defined(__FreeBSD__) && defined(__aarch64__)
    u_long cap;
    if (elf_aux_info(AT_HWCAP, &cap, sizeof cap) != 0) return 0;
    return (cap & HWCAP_ASIMD) != 0;
#elif defined(_WIN32)
    return IsProcessorFeaturePresent(PF_ARM_V8_INSTRUCTIONS_AVAILABLE) != 0;
#else
    return 0;
#endif
}

static inline
unsigned char *(*rans_dec_func(int do_simd, int order))
    (unsigned char *in,
     unsigned int in_size,
     unsigned char *out,
     unsigned int out_size) {

    if (do_simd) {
        if ((rans_cpu & RANS_CPU_DEC_NEON) && have_neon())
            return order & 1
                ? rans_uncompress_O1_32x16_neon
                : rans_uncompress_O0_32x16_neon;
        else
            return order & 1
                ? rans_uncompress_O1_32x16
                : rans_uncompress_O0_32x16;
    } else {
        return order & 1
            ? rans_uncompress_O1_4x16
            : rans_uncompress_O0_4x16;
    }
}

#else // !(defined(__GNUC__) && defined(__x86_64__)) && !defined(__ARM_NEON)

static inline
unsigned char *(*rans_dec_func(int do_simd, int order))
    (unsigned char *in,
     unsigned int in_size,
     unsigned char *out,
     unsigned int out_size) {

    if (do_simd) {
        return order & 1
            ? rans_uncompress_O1_32x16
            : rans_uncompress_O0_32x16;
    } else {
        return order & 1
            ? rans_uncompress_O1_4x16
            : rans_uncompress_O0_4x16;
    }
}

#endif

unsigned char *rans_uncompress_to_4x16(unsigned char *in,  unsigned int in_size,
                                       unsigned char *out, unsigned int *out_size) {
    unsigned char *in_end = in + in_size;
    unsigned char *out_free = NULL, *tmp_free = NULL, *meta_free = NULL;

    if (in_size == 0)
        return NULL;

    if (*in & RANS_ORDER_STRIPE) {
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
#ifdef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
            if (ulen > 100000)
                return NULL;
#endif
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
            if (!rans_uncompress_to_4x16(in+c_meta_len, in_size-c_meta_len, outN + idxN[i], &olen)
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
    int do_pack = order & RANS_ORDER_PACK;
    int do_rle  = order & RANS_ORDER_RLE;
    int do_cat  = order & RANS_ORDER_CAT;
    int no_size = order & RANS_ORDER_NOSZ;
    int do_simd = order & RANS_ORDER_X32;
    order &= 1;

    int sz = 0;
    unsigned int osz;
    if (!no_size) {
        sz = var_get_u32(in, in_end, &osz);
    } else
        sz = 0, osz = *out_size;
    in += sz;
    in_size -= sz;

#ifdef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
    if (osz > 100000)
        return NULL;
#endif

    if (no_size && !out)
        goto err; // Need one or the other

    if (!out) {
        *out_size = osz;
        if (!(out = out_free = malloc(*out_size)))
            return NULL;
    } else {
        if (*out_size < osz)
        goto err;
        *out_size = osz;
    }

//    if (do_pack || do_rle) {
//      in += sz; // size field not needed when pure rANS
//      in_size -= sz;
//    }

    uint32_t c_meta_size = 0;
    unsigned int tmp1_size = *out_size;
    unsigned int tmp2_size = *out_size;
    unsigned int tmp3_size = *out_size;
    unsigned char *tmp1 = NULL, *tmp2 = NULL, *tmp3 = NULL, *tmp = NULL;

    // Need In, Out and Tmp buffers with temporary buffer of the same size
    // as output.  All use rANS, but with optional transforms (none, RLE,
    // Pack, or both).
    //
    //                    rans   unrle  unpack
    // If none:     in -> out
    // If RLE:      in -> tmp -> out
    // If Pack:     in -> tmp        -> out
    // If RLE+Pack: in -> out -> tmp -> out
    //                    tmp1   tmp2   tmp3
    //
    // So rans is in   -> tmp1
    // RLE     is tmp1 -> tmp2
    // Unpack  is tmp2 -> tmp3

    // Format is meta data (Pack and RLE in that order if present),
    // followed by rANS compressed data.

    if (do_pack || do_rle) {
        if (!(tmp = tmp_free = malloc(*out_size)))
            goto err;
        if (do_pack && do_rle) {
            tmp1 = out;
            tmp2 = tmp;
            tmp3 = out;
        } else if (do_pack) {
            tmp1 = tmp;
            tmp2 = tmp1;
            tmp3 = out;
        } else if (do_rle) {
            tmp1 = tmp;
            tmp2 = out;
            tmp3 = out;
        }
    } else {
        // neither
        tmp  = NULL;
        tmp1 = out;
        tmp2 = out;
        tmp3 = out;
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

    uint8_t *meta = NULL;
    uint32_t u_meta_size = 0;
    if (do_rle) {
        // Uncompress meta data
        uint32_t c_meta_size, rle_len, sz;
        sz  = var_get_u32(in,    in_end, &u_meta_size);
        sz += var_get_u32(in+sz, in_end, &rle_len);
        if (rle_len > tmp1_size) // should never grow
            goto err;
        if (u_meta_size & 1) {
            meta = in + sz;
            u_meta_size = u_meta_size/2 > (in_end-meta) ? (in_end-meta) : u_meta_size/2;
            c_meta_size = u_meta_size;
        } else {
            sz += var_get_u32(in+sz, in_end, &c_meta_size);
            u_meta_size /= 2;

            meta_free = meta = rans_dec_func(do_simd, 0)(in+sz, in_size-sz, NULL, u_meta_size);
            if (!meta)
                goto err;
        }
        if (c_meta_size+sz > in_size)
            goto err;
        in      += c_meta_size+sz;
        in_size -= c_meta_size+sz;
        tmp1_size = rle_len;
    }
    //fprintf(stderr, "    meta_size %d bytes\n", (int)(in - orig_in)); //c-size

    // uncompress RLE data.  in -> tmp1
    if (in_size) {
        if (do_cat) {
            //fprintf(stderr, "    CAT %d\n", tmp1_size); //c-size
            if (tmp1_size > in_size)
                goto err;
            if (tmp1_size > *out_size)
                goto err;
            memcpy(tmp1, in, tmp1_size);
        } else {
            tmp1 = rans_dec_func(do_simd, order)(in, in_size, tmp1, tmp1_size);
            if (!tmp1)
                goto err;
        }
    } else {
        tmp1_size = 0;
    }
    tmp2_size = tmp3_size = tmp1_size;

    if (do_rle) {
        // Unpack RLE.  tmp1 -> tmp2.
        if (u_meta_size == 0)
            goto err;
        uint64_t unrle_size = *out_size;
        int rle_nsyms = *meta ? *meta : 256;
        if (u_meta_size < 1+rle_nsyms)
            goto err;
        if (!hts_rle_decode(tmp1, tmp1_size,
                            meta+1+rle_nsyms, u_meta_size-(1+rle_nsyms),
                            meta+1, rle_nsyms, tmp2, &unrle_size))
            goto err;
        tmp3_size = tmp2_size = unrle_size;
        free(meta_free);
        meta_free = NULL;
    }
    if (do_pack) {
        // Unpack bits via pack-map.  tmp2 -> tmp3
        if (npacked_sym == 1)
            unpacked_sz = tmp2_size;
        //uint8_t *porig = unpack(tmp2, tmp2_size, unpacked_sz, npacked_sym, map);
        //memcpy(tmp3, porig, unpacked_sz);
        if (!hts_unpack(tmp2, tmp2_size, tmp3, unpacked_sz, npacked_sym, map))
            goto err;
        tmp3_size = unpacked_sz;
    }

    if (tmp)
        free(tmp);

    *out_size = tmp3_size;
    return tmp3;

 err:
    free(meta_free);
    free(out_free);
    free(tmp_free);
    return NULL;
}

unsigned char *rans_uncompress_4x16(unsigned char *in, unsigned int in_size,
                                    unsigned int *out_size) {
    return rans_uncompress_to_4x16(in, in_size, NULL, out_size);
}
