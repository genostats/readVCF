/*
 * Copyright (c) 2019-2021 Genome Research Ltd.
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
#include <string.h>
#include <stdio.h>

#include "htscodecs/varint.h"
#include "htscodecs/rle.h"


// On input *out_len holds the allocated size of out[].
// On output it holds the used size of out[].
uint8_t *hts_rle_decode(uint8_t *lit, uint64_t lit_len,
                        uint8_t *run, uint64_t run_len,
                        uint8_t *rle_syms, int rle_nsyms,
                        uint8_t *out, uint64_t *out_len) {
    uint64_t j;
    uint8_t *run_end = run + run_len;

#ifdef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
    if (*out_len > 100000)
        return NULL;
#endif

    int saved[256] = {0};
    for (j = 0; j < rle_nsyms; j++)
        saved[rle_syms[j]] = 1;

    uint8_t *lit_end = lit + lit_len;
    uint8_t *out_end = out + *out_len;
    uint8_t *outp = out;

    while (lit < lit_end) {
        if (outp >= out_end)
            goto err;

        uint8_t b = *lit;
        if (saved[b]) {
            uint32_t rlen;
            run += var_get_u32(run, run_end, &rlen);
            if (rlen) {
                if (outp + rlen >= out_end)
                    goto err;
                memset(outp, b, rlen+1);
                outp += rlen+1;
            } else {
                *outp++ = b;
            }
        } else {
            *outp++ = b;
        }
        lit++;
    }

    *out_len = outp-out;
    return out;

 err:
    return NULL;
}