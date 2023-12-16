/*
Copyright (c) 2012-2023 Genome Research Ltd.
Author: James Bonfield <jkb@sanger.ac.uk>

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

   3. Neither the names Genome Research Ltd and Wellcome Trust Sanger
Institute nor the names of its contributors may be used to endorse or promote
products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH LTD OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/*
 * CRAM I/O primitives.
 *
 * - ITF8 encoding and decoding.
 * - Block based I/O
 * - Zlib inflating and deflating (memory)
 * - CRAM basic data structure reading and writing
 * - File opening / closing
 * - Reference sequence handling
 */

#define HTS_BUILDING_LIBRARY // Enables HTSLIB_EXPORT, see htslib/hts_defs.h
#include <config.h>

#include <stdio.h>
#include <errno.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>
#include <bzlib.h>
#ifdef HAVE_LZMA_H
#include <lzma.h>
#else
#include "../os/lzma_stub.h"
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>
#include <stdint.h>

#ifdef HAVE_LIBDEFLATE
#include <libdeflate.h>
#define crc32(a,b,c) libdeflate_crc32((a),(b),(c))
#endif

#ifdef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
#include "../fuzz_settings.h"
#endif

#include "cram/cram.h"
#include "cram/os.h"
#include "htslib/hts.h"
#include "cram/open_trace_file.h"

#include "htscodecs/rANS_static.h"
#include "htscodecs/rANS_static4x16.h"
#include "htscodecs/arith_dynamic.h"
#include "htscodecs/varint.h"

//#define REF_DEBUG

#ifdef REF_DEBUG
#include <sys/syscall.h>
#define gettid() (int)syscall(SYS_gettid)

#define RP(...) fprintf (stderr, __VA_ARGS__)
#else
#define RP(...)
#endif

#include "htslib/hfile.h"
#include "htslib/bgzf.h"
#include "htslib/faidx.h"
#include "hts_internal.h"

#ifndef PATH_MAX
#define PATH_MAX FILENAME_MAX
#endif

#define TRIAL_SPAN 70
#define NTRIALS 3

#define CRAM_DEFAULT_LEVEL 5

/* ----------------------------------------------------------------------
 * ITF8 encoding and decoding.
 *
 * Also see the itf8_get and itf8_put macros in cram_io.h
 */

/*
 * LEGACY: consider using itf8_decode_crc.
 *
 * Reads an integer in ITF-8 encoding from 'cp' and stores it in
 * *val.
 *
 * Returns the number of bytes read on success
 *        -1 on failure
 */
int itf8_decode(cram_fd *fd, int32_t *val_p) {
    static int nbytes[16] = {
        0,0,0,0, 0,0,0,0,                               // 0000xxxx - 0111xxxx
        1,1,1,1,                                        // 1000xxxx - 1011xxxx
        2,2,                                            // 1100xxxx - 1101xxxx
        3,                                              // 1110xxxx
        4,                                              // 1111xxxx
    };

    static int nbits[16] = {
        0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f, // 0000xxxx - 0111xxxx
        0x3f, 0x3f, 0x3f, 0x3f,                         // 1000xxxx - 1011xxxx
        0x1f, 0x1f,                                     // 1100xxxx - 1101xxxx
        0x0f,                                           // 1110xxxx
        0x0f,                                           // 1111xxxx
    };

    int32_t val = hgetc(fd->fp);
    if (val == -1)
        return -1;

    int i = nbytes[val>>4];
    val &= nbits[val>>4];

    switch(i) {
    case 0:
        *val_p = val;
        return 1;

    case 1:
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        *val_p = val;
        return 2;

    case 2:
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        *val_p = val;
        return 3;

    case 3:
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        *val_p = val;
        return 4;

    case 4: // really 3.5 more, why make it different?
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<4) | (((unsigned char)hgetc(fd->fp)) & 0x0f);
        *val_p = val;
    }

    return 5;
}

int itf8_decode_crc(cram_fd *fd, int32_t *val_p, uint32_t *crc) {
    static int nbytes[16] = {
        0,0,0,0, 0,0,0,0,                               // 0000xxxx - 0111xxxx
        1,1,1,1,                                        // 1000xxxx - 1011xxxx
        2,2,                                            // 1100xxxx - 1101xxxx
        3,                                              // 1110xxxx
        4,                                              // 1111xxxx
    };

    static int nbits[16] = {
        0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f, // 0000xxxx - 0111xxxx
        0x3f, 0x3f, 0x3f, 0x3f,                         // 1000xxxx - 1011xxxx
        0x1f, 0x1f,                                     // 1100xxxx - 1101xxxx
        0x0f,                                           // 1110xxxx
        0x0f,                                           // 1111xxxx
    };
    unsigned char c[5];

    int32_t val = hgetc(fd->fp);
    if (val == -1)
        return -1;

    c[0]=val;

    int i = nbytes[val>>4];
    val &= nbits[val>>4];

    if (i > 0) {
        if (hread(fd->fp, &c[1], i) < i)
            return -1;
    }

    switch(i) {
    case 0:
        *val_p = val;
        *crc = crc32(*crc, c, 1);
        return 1;

    case 1:
        val = (val<<8) | c[1];
        *val_p = val;
        *crc = crc32(*crc, c, 2);
        return 2;

    case 2:
        val = (val<<8) | c[1];
        val = (val<<8) | c[2];
        *val_p = val;
        *crc = crc32(*crc, c, 3);
        return 3;

    case 3:
        val = (val<<8) | c[1];
        val = (val<<8) | c[2];
        val = (val<<8) | c[3];
        *val_p = val;
        *crc = crc32(*crc, c, 4);
        return 4;

    case 4: // really 3.5 more, why make it different?
        {
            uint32_t uv = val;
            uv = (uv<<8) |  c[1];
            uv = (uv<<8) |  c[2];
            uv = (uv<<8) |  c[3];
            uv = (uv<<4) | (c[4] & 0x0f);
            // Avoid implementation-defined behaviour on negative values
            *val_p = uv < 0x80000000UL ? (int32_t) uv : -((int32_t) (0xffffffffUL - uv)) - 1;
            *crc = crc32(*crc, c, 5);
        }
    }

    return 5;
}

/*
 * Stores a value to memory in ITF-8 format.
 *
 * Returns the number of bytes required to store the number.
 * This is a maximum of 5 bytes.
 */
static inline int itf8_put(char *cp, int32_t val) {
    unsigned char *up = (unsigned char *)cp;
    if        (!(val & ~0x00000007f)) { // 1 byte
        *up = val;
        return 1;
    } else if (!(val & ~0x00003fff)) { // 2 byte
        *up++ = (val >> 8 ) | 0x80;
        *up   = val & 0xff;
        return 2;
    } else if (!(val & ~0x01fffff)) { // 3 byte
        *up++ = (val >> 16) | 0xc0;
        *up++ = (val >> 8 ) & 0xff;
        *up   = val & 0xff;
        return 3;
    } else if (!(val & ~0x0fffffff)) { // 4 byte
        *up++ = (val >> 24) | 0xe0;
        *up++ = (val >> 16) & 0xff;
        *up++ = (val >> 8 ) & 0xff;
        *up   = val & 0xff;
        return 4;
    } else {                           // 5 byte
        *up++ = 0xf0 | ((val>>28) & 0xff);
        *up++ = (val >> 20) & 0xff;
        *up++ = (val >> 12) & 0xff;
        *up++ = (val >> 4 ) & 0xff;
        *up = val & 0x0f;
        return 5;
    }
}


/* 64-bit itf8 variant */
static inline int ltf8_put(char *cp, int64_t val) {
    unsigned char *up = (unsigned char *)cp;
    if        (!(val & ~((1LL<<7)-1))) {
        *up = val;
        return 1;
    } else if (!(val & ~((1LL<<(6+8))-1))) {
        *up++ = (val >> 8 ) | 0x80;
        *up   = val & 0xff;
        return 2;
    } else if (!(val & ~((1LL<<(5+2*8))-1))) {
        *up++ = (val >> 16) | 0xc0;
        *up++ = (val >> 8 ) & 0xff;
        *up   = val & 0xff;
        return 3;
    } else if (!(val & ~((1LL<<(4+3*8))-1))) {
        *up++ = (val >> 24) | 0xe0;
        *up++ = (val >> 16) & 0xff;
        *up++ = (val >> 8 ) & 0xff;
        *up   = val & 0xff;
        return 4;
    } else if (!(val & ~((1LL<<(3+4*8))-1))) {
        *up++ = (val >> 32) | 0xf0;
        *up++ = (val >> 24) & 0xff;
        *up++ = (val >> 16) & 0xff;
        *up++ = (val >> 8 ) & 0xff;
        *up   = val & 0xff;
        return 5;
    } else if (!(val & ~((1LL<<(2+5*8))-1))) {
        *up++ = (val >> 40) | 0xf8;
        *up++ = (val >> 32) & 0xff;
        *up++ = (val >> 24) & 0xff;
        *up++ = (val >> 16) & 0xff;
        *up++ = (val >> 8 ) & 0xff;
        *up   = val & 0xff;
        return 6;
    } else if (!(val & ~((1LL<<(1+6*8))-1))) {
        *up++ = (val >> 48) | 0xfc;
        *up++ = (val >> 40) & 0xff;
        *up++ = (val >> 32) & 0xff;
        *up++ = (val >> 24) & 0xff;
        *up++ = (val >> 16) & 0xff;
        *up++ = (val >> 8 ) & 0xff;
        *up   = val & 0xff;
        return 7;
    } else if (!(val & ~((1LL<<(7*8))-1))) {
        *up++ = (val >> 56) | 0xfe;
        *up++ = (val >> 48) & 0xff;
        *up++ = (val >> 40) & 0xff;
        *up++ = (val >> 32) & 0xff;
        *up++ = (val >> 24) & 0xff;
        *up++ = (val >> 16) & 0xff;
        *up++ = (val >> 8 ) & 0xff;
        *up   = val & 0xff;
        return 8;
    } else {
        *up++ = 0xff;
        *up++ = (val >> 56) & 0xff;
        *up++ = (val >> 48) & 0xff;
        *up++ = (val >> 40) & 0xff;
        *up++ = (val >> 32) & 0xff;
        *up++ = (val >> 24) & 0xff;
        *up++ = (val >> 16) & 0xff;
        *up++ = (val >> 8 ) & 0xff;
        *up   = val & 0xff;
        return 9;
    }
}

/*
 * Encodes and writes a single integer in ITF-8 format.
 * Returns 0 on success
 *        -1 on failure
 */
int itf8_encode(cram_fd *fd, int32_t val) {
    char buf[5];
    int len = itf8_put(buf, val);
    return hwrite(fd->fp, buf, len) == len ? 0 : -1;
}

const int itf8_bytes[16] = {
    1, 1, 1, 1,  1, 1, 1, 1,
    2, 2, 2, 2,  3, 3, 4, 5
};

const int ltf8_bytes[256] = {
    1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,
    1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,
    1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,
    1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,

    1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,
    1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,
    1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,
    1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,

    2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,
    2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,
    2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,
    2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,

    3, 3, 3, 3,  3, 3, 3, 3,  3, 3, 3, 3,  3, 3, 3, 3,
    3, 3, 3, 3,  3, 3, 3, 3,  3, 3, 3, 3,  3, 3, 3, 3,

    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,

    5, 5, 5, 5,  5, 5, 5, 5,  6, 6, 6, 6,  7, 7, 8, 9
};

/*
 * LEGACY: consider using ltf8_decode_crc.
 */
int ltf8_decode(cram_fd *fd, int64_t *val_p) {
    int c = hgetc(fd->fp);
    int64_t val = (unsigned char)c;
    if (c == -1)
        return -1;

    if (val < 0x80) {
        *val_p =   val;
        return 1;

    } else if (val < 0xc0) {
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        *val_p = val & (((1LL<<(6+8)))-1);
        return 2;

    } else if (val < 0xe0) {
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        *val_p = val & ((1LL<<(5+2*8))-1);
        return 3;

    } else if (val < 0xf0) {
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        *val_p = val & ((1LL<<(4+3*8))-1);
        return 4;

    } else if (val < 0xf8) {
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        *val_p = val & ((1LL<<(3+4*8))-1);
        return 5;

    } else if (val < 0xfc) {
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        *val_p = val & ((1LL<<(2+5*8))-1);
        return 6;

    } else if (val < 0xfe) {
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        *val_p = val & ((1LL<<(1+6*8))-1);
        return 7;

    } else if (val < 0xff) {
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        *val_p = val & ((1LL<<(7*8))-1);
        return 8;

    } else {
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        *val_p = val;
    }

    return 9;
}

int ltf8_decode_crc(cram_fd *fd, int64_t *val_p, uint32_t *crc) {
    unsigned char c[9];
    int64_t val = hgetc(fd->fp);
    if (val < 0)
        return -1;

    c[0] = val;

    if (val < 0x80) {
        *val_p =   val;
        *crc = crc32(*crc, c, 1);
        return 1;

    } else if (val < 0xc0) {
        int v = hgetc(fd->fp);
        if (v < 0)
            return -1;
        val = (val<<8) | (c[1]=v);
        *val_p = val & (((1LL<<(6+8)))-1);
        *crc = crc32(*crc, c, 2);
        return 2;

    } else if (val < 0xe0) {
        if (hread(fd->fp, &c[1], 2) < 2)
            return -1;
        val = (val<<8) | c[1];
        val = (val<<8) | c[2];
        *val_p = val & ((1LL<<(5+2*8))-1);
        *crc = crc32(*crc, c, 3);
        return 3;

    } else if (val < 0xf0) {
        if (hread(fd->fp, &c[1], 3) < 3)
            return -1;
        val = (val<<8) | c[1];
        val = (val<<8) | c[2];
        val = (val<<8) | c[3];
        *val_p = val & ((1LL<<(4+3*8))-1);
        *crc = crc32(*crc, c, 4);
        return 4;

    } else if (val < 0xf8) {
        if (hread(fd->fp, &c[1], 4) < 4)
            return -1;
        val = (val<<8) | c[1];
        val = (val<<8) | c[2];
        val = (val<<8) | c[3];
        val = (val<<8) | c[4];
        *val_p = val & ((1LL<<(3+4*8))-1);
        *crc = crc32(*crc, c, 5);
        return 5;

    } else if (val < 0xfc) {
        if (hread(fd->fp, &c[1], 5) < 5)
            return -1;
        val = (val<<8) | c[1];
        val = (val<<8) | c[2];
        val = (val<<8) | c[3];
        val = (val<<8) | c[4];
        val = (val<<8) | c[5];
        *val_p = val & ((1LL<<(2+5*8))-1);
        *crc = crc32(*crc, c, 6);
        return 6;

    } else if (val < 0xfe) {
        if (hread(fd->fp, &c[1], 6) < 6)
            return -1;
        val = (val<<8) | c[1];
        val = (val<<8) | c[2];
        val = (val<<8) | c[3];
        val = (val<<8) | c[4];
        val = (val<<8) | c[5];
        val = (val<<8) | c[6];
        *val_p = val & ((1LL<<(1+6*8))-1);
        *crc = crc32(*crc, c, 7);
        return 7;

    } else if (val < 0xff) {
        uint64_t uval = val;
        if (hread(fd->fp, &c[1], 7) < 7)
            return -1;
        uval = (uval<<8) | c[1];
        uval = (uval<<8) | c[2];
        uval = (uval<<8) | c[3];
        uval = (uval<<8) | c[4];
        uval = (uval<<8) | c[5];
        uval = (uval<<8) | c[6];
        uval = (uval<<8) | c[7];
        *val_p = uval & ((1ULL<<(7*8))-1);
        *crc = crc32(*crc, c, 8);
        return 8;

    } else {
        uint64_t uval;
        if (hread(fd->fp, &c[1], 8) < 8)
            return -1;
        uval =             c[1];
        uval = (uval<<8) | c[2];
        uval = (uval<<8) | c[3];
        uval = (uval<<8) | c[4];
        uval = (uval<<8) | c[5];
        uval = (uval<<8) | c[6];
        uval = (uval<<8) | c[7];
        uval = (uval<<8) | c[8];
        *crc = crc32(*crc, c, 9);
        // Avoid implementation-defined behaviour on negative values
        *val_p = c[1] < 0x80 ? (int64_t) uval : -((int64_t) (0xffffffffffffffffULL - uval)) - 1;
    }

    return 9;
}

/*
 * Pushes a value in ITF8 format onto the end of a block.
 * This shouldn't be used for high-volume data as it is not the fastest
 * method.
 *
 * Returns the number of bytes written
 */
int itf8_put_blk(cram_block *blk, int32_t val) {
    char buf[5];
    int sz;

    sz = itf8_put(buf, val);
    BLOCK_APPEND(blk, buf, sz);
    return sz;

 block_err:
    return -1;
}

int ltf8_put_blk(cram_block *blk, int64_t val) {
    char buf[9];
    int sz;

    sz = ltf8_put(buf, val);
    BLOCK_APPEND(blk, buf, sz);
    return sz;

 block_err:
    return -1;
}

static int64_t safe_itf8_get(char **cp, const char *endp, int *err) {
    const unsigned char *up = (unsigned char *)*cp;

    if (endp && endp - *cp < 5 &&
        (*cp >= endp || endp - *cp < itf8_bytes[up[0]>>4])) {
        if (err) *err = 1;
        return 0;
    }

    if (up[0] < 0x80) {
        (*cp)++;
        return up[0];
    } else if (up[0] < 0xc0) {
        (*cp)+=2;
        return ((up[0] <<8) |  up[1])                           & 0x3fff;
    } else if (up[0] < 0xe0) {
        (*cp)+=3;
        return ((up[0]<<16) | (up[1]<< 8) |  up[2])             & 0x1fffff;
    } else if (up[0] < 0xf0) {
        (*cp)+=4;
        uint32_t uv = (((uint32_t)up[0]<<24) | (up[1]<<16) | (up[2]<<8) | up[3]) & 0x0fffffff;
        return (int32_t)uv;
    } else {
        (*cp)+=5;
        uint32_t uv = (((uint32_t)up[0] & 0x0f)<<28) | (up[1]<<20) | (up[2]<<12) | (up[3]<<4) | (up[4] & 0x0f);
        return (int32_t)uv;
    }
}

static int64_t safe_ltf8_get(char **cp, const char *endp, int *err) {
    unsigned char *up = (unsigned char *)*cp;

    if (endp && endp - *cp < 9 &&
        (*cp >= endp || endp - *cp < ltf8_bytes[up[0]])) {
        if (err) *err = 1;
        return 0;
    }

    if (up[0] < 0x80) {
        (*cp)++;
        return up[0];
    } else if (up[0] < 0xc0) {
        (*cp)+=2;
        return (((uint64_t)up[0]<< 8) |
                 (uint64_t)up[1]) & (((1LL<<(6+8)))-1);
    } else if (up[0] < 0xe0) {
        (*cp)+=3;
        return (((uint64_t)up[0]<<16) |
                ((uint64_t)up[1]<< 8) |
                 (uint64_t)up[2]) & ((1LL<<(5+2*8))-1);
    } else if (up[0] < 0xf0) {
        (*cp)+=4;
        return (((uint64_t)up[0]<<24) |
                ((uint64_t)up[1]<<16) |
                ((uint64_t)up[2]<< 8) |
                 (uint64_t)up[3]) & ((1LL<<(4+3*8))-1);
    } else if (up[0] < 0xf8) {
        (*cp)+=5;
        return (((uint64_t)up[0]<<32) |
                ((uint64_t)up[1]<<24) |
                ((uint64_t)up[2]<<16) |
                ((uint64_t)up[3]<< 8) |
                 (uint64_t)up[4]) & ((1LL<<(3+4*8))-1);
    } else if (up[0] < 0xfc) {
        (*cp)+=6;
        return (((uint64_t)up[0]<<40) |
                ((uint64_t)up[1]<<32) |
                ((uint64_t)up[2]<<24) |
                ((uint64_t)up[3]<<16) |
                ((uint64_t)up[4]<< 8) |
                 (uint64_t)up[5]) & ((1LL<<(2+5*8))-1);
    } else if (up[0] < 0xfe) {
        (*cp)+=7;
        return (((uint64_t)up[0]<<48) |
                ((uint64_t)up[1]<<40) |
                ((uint64_t)up[2]<<32) |
                ((uint64_t)up[3]<<24) |
                ((uint64_t)up[4]<<16) |
                ((uint64_t)up[5]<< 8) |
                 (uint64_t)up[6]) & ((1LL<<(1+6*8))-1);
    } else if (up[0] < 0xff) {
        (*cp)+=8;
        return (((uint64_t)up[1]<<48) |
                ((uint64_t)up[2]<<40) |
                ((uint64_t)up[3]<<32) |
                ((uint64_t)up[4]<<24) |
                ((uint64_t)up[5]<<16) |
                ((uint64_t)up[6]<< 8) |
                 (uint64_t)up[7]) & ((1LL<<(7*8))-1);
    } else {
        (*cp)+=9;
        return (((uint64_t)up[1]<<56) |
                ((uint64_t)up[2]<<48) |
                ((uint64_t)up[3]<<40) |
                ((uint64_t)up[4]<<32) |
                ((uint64_t)up[5]<<24) |
                ((uint64_t)up[6]<<16) |
                ((uint64_t)up[7]<< 8) |
                 (uint64_t)up[8]);
    }
}

// Wrapper for now
static int safe_itf8_put(char *cp, char *cp_end, int32_t val) {
    return itf8_put(cp, val);
}

static int safe_ltf8_put(char *cp, char *cp_end, int64_t val) {
    return ltf8_put(cp, val);
}

static int itf8_size(int64_t v) {
    return ((!((v)&~0x7f))?1:(!((v)&~0x3fff))?2:(!((v)&~0x1fffff))?3:(!((v)&~0xfffffff))?4:5);
}

//-----------------------------------------------------------------------------

// CRAM v4.0 onwards uses a different variable sized integer encoding
// that is size agnostic.

// Local interface to varint.h inline version, so we can use in func ptr.
// Note a lot of these use the unsigned interface but take signed int64_t.
// This is because the old CRAM ITF8 inteface had signed -1 as unsigned
// 0xffffffff.
static int uint7_size(int64_t v) {
    return var_size_u64(v);
}

static int64_t uint7_get_32(char **cp, const char *endp, int *err) {
    uint32_t val;
    int nb = var_get_u32((uint8_t *)(*cp), (const uint8_t *)endp, &val);
    (*cp) += nb;
    if (!nb && err) *err = 1;
    return val;
}

static int64_t sint7_get_32(char **cp, const char *endp, int *err) {
    int32_t val;
    int nb = var_get_s32((uint8_t *)(*cp), (const uint8_t *)endp, &val);
    (*cp) += nb;
    if (!nb && err) *err = 1;
    return val;
}

static int64_t uint7_get_64(char **cp, const char *endp, int *err) {
    uint64_t val;
    int nb = var_get_u64((uint8_t *)(*cp), (const uint8_t *)endp, &val);
    (*cp) += nb;
    if (!nb && err) *err = 1;
    return val;
}

static int64_t sint7_get_64(char **cp, const char *endp, int *err) {
    int64_t val;
    int nb = var_get_s64((uint8_t *)(*cp), (const uint8_t *)endp, &val);
    (*cp) += nb;
    if (!nb && err) *err = 1;
    return val;
}

static int uint7_put_32(char *cp, char *endp, int32_t val) {
    return var_put_u32((uint8_t *)cp, (uint8_t *)endp, val);
}

static int sint7_put_32(char *cp, char *endp, int32_t val) {
    return var_put_s32((uint8_t *)cp, (uint8_t *)endp, val);
}

static int uint7_put_64(char *cp, char *endp, int64_t val) {
    return var_put_u64((uint8_t *)cp, (uint8_t *)endp, val);
}

static int sint7_put_64(char *cp, char *endp, int64_t val) {
    return var_put_s64((uint8_t *)cp, (uint8_t *)endp, val);
}

// Put direct to to cram_block
static int uint7_put_blk_32(cram_block *blk, int32_t v) {
    uint8_t buf[10];
    int sz = var_put_u32(buf, buf+10, v);
    BLOCK_APPEND(blk, buf, sz);
    return sz;

 block_err:
    return -1;
}

static int sint7_put_blk_32(cram_block *blk, int32_t v) {
    uint8_t buf[10];
    int sz = var_put_s32(buf, buf+10, v);
    BLOCK_APPEND(blk, buf, sz);
    return sz;

 block_err:
    return -1;
}

static int uint7_put_blk_64(cram_block *blk, int64_t v) {
    uint8_t buf[10];
    int sz = var_put_u64(buf, buf+10, v);
    BLOCK_APPEND(blk, buf, sz);
    return sz;

 block_err:
    return -1;
}

static int sint7_put_blk_64(cram_block *blk, int64_t v) {
    uint8_t buf[10];
    int sz = var_put_s64(buf, buf+10, v);
    BLOCK_APPEND(blk, buf, sz);
    return sz;

 block_err:
    return -1;
}

// Decode 32-bits with CRC update from cram_fd
static int uint7_decode_crc32(cram_fd *fd, int32_t *val_p, uint32_t *crc) {
    uint8_t b[5], i = 0;
    int c;
    uint32_t v = 0;

#ifdef VARINT2
    b[0] = hgetc(fd->fp);
    if (b[0] < 177) {
    } else if (b[0] < 241) {
        b[1] = hgetc(fd->fp);
    } else if (b[0] < 249) {
        b[1] = hgetc(fd->fp);
        b[2] = hgetc(fd->fp);
    } else {
        int n = b[0]+2, z = 1;
        while (n-- >= 249)
            b[z++] = hgetc(fd->fp);
    }
    i = var_get_u32(b, NULL, &v);
#else
//    // Little endian
//    int s = 0;
//    do {
//        b[i++] = c = hgetc(fd->fp);
//        if (c < 0)
//            return -1;
//        v |= (c & 0x7f) << s;
//      s += 7;
//    } while (i < 5 && (c & 0x80));

    // Big endian, see also htscodecs/varint.h
    do {
        b[i++] = c = hgetc(fd->fp);
        if (c < 0)
            return -1;
        v = (v<<7) | (c & 0x7f);
    } while (i < 5 && (c & 0x80));
#endif
    *crc = crc32(*crc, b, i);

    *val_p = v;
    return i;
}

// Decode 32-bits with CRC update from cram_fd
static int sint7_decode_crc32(cram_fd *fd, int32_t *val_p, uint32_t *crc) {
    uint8_t b[5], i = 0;
    int c;
    uint32_t v = 0;

#ifdef VARINT2
    b[0] = hgetc(fd->fp);
    if (b[0] < 177) {
    } else if (b[0] < 241) {
        b[1] = hgetc(fd->fp);
    } else if (b[0] < 249) {
        b[1] = hgetc(fd->fp);
        b[2] = hgetc(fd->fp);
    } else {
        int n = b[0]+2, z = 1;
        while (n-- >= 249)
            b[z++] = hgetc(fd->fp);
    }
    i = var_get_u32(b, NULL, &v);
#else
//    // Little endian
//    int s = 0;
//    do {
//        b[i++] = c = hgetc(fd->fp);
//        if (c < 0)
//            return -1;
//        v |= (c & 0x7f) << s;
//      s += 7;
//    } while (i < 5 && (c & 0x80));

    // Big endian, see also htscodecs/varint.h
    do {
        b[i++] = c = hgetc(fd->fp);
        if (c < 0)
            return -1;
        v = (v<<7) | (c & 0x7f);
    } while (i < 5 && (c & 0x80));
#endif
    *crc = crc32(*crc, b, i);

    *val_p = (v>>1) ^ -(v&1);
    return i;
}


// Decode 64-bits with CRC update from cram_fd
static int uint7_decode_crc64(cram_fd *fd, int64_t *val_p, uint32_t *crc) {
    uint8_t b[10], i = 0;
    int c;
    uint64_t v = 0;

#ifdef VARINT2
    b[0] = hgetc(fd->fp);
    if (b[0] < 177) {
    } else if (b[0] < 241) {
        b[1] = hgetc(fd->fp);
    } else if (b[0] < 249) {
        b[1] = hgetc(fd->fp);
        b[2] = hgetc(fd->fp);
    } else {
        int n = b[0]+2, z = 1;
        while (n-- >= 249)
            b[z++] = hgetc(fd->fp);
    }
    i = var_get_u64(b, NULL, &v);
#else
//    // Little endian
//    int s = 0;
//    do {
//        b[i++] = c = hgetc(fd->fp);
//        if (c < 0)
//            return -1;
//        v |= (c & 0x7f) << s;
//      s += 7;
//    } while (i < 10 && (c & 0x80));

    // Big endian, see also htscodecs/varint.h
    do {
        b[i++] = c = hgetc(fd->fp);
        if (c < 0)
            return -1;
        v = (v<<7) | (c & 0x7f);
    } while (i < 5 && (c & 0x80));
#endif
    *crc = crc32(*crc, b, i);

    *val_p = v;
    return i;
}

//-----------------------------------------------------------------------------

/*
 * Decodes a 32-bit little endian value from fd and stores in val.
 *
 * Returns the number of bytes read on success
 *         -1 on failure
 */
static int int32_decode(cram_fd *fd, int32_t *val) {
    int32_t i;
    if (4 != hread(fd->fp, &i, 4))
        return -1;

    *val = le_int4(i);
    return 4;
}

/*
 * Encodes a 32-bit little endian value 'val' and writes to fd.
 *
 * Returns the number of bytes written on success
 *         -1 on failure
 */
static int int32_encode(cram_fd *fd, int32_t val) {
    uint32_t v = le_int4(val);
    if (4 != hwrite(fd->fp, &v, 4))
        return -1;

    return 4;
}

/* As int32_decoded/encode, but from/to blocks instead of cram_fd */
int int32_get_blk(cram_block *b, int32_t *val) {
    if (b->uncomp_size - BLOCK_SIZE(b) < 4)
        return -1;

    uint32_t v =
         ((uint32_t) b->data[b->byte  ])        |
        (((uint32_t) b->data[b->byte+1]) <<  8) |
        (((uint32_t) b->data[b->byte+2]) << 16) |
        (((uint32_t) b->data[b->byte+3]) << 24);
    // Avoid implementation-defined behaviour on negative values
    *val = v < 0x80000000U ? (int32_t) v : -((int32_t) (0xffffffffU - v)) - 1;
    BLOCK_SIZE(b) += 4;
    return 4;
}


#ifdef HAVE_LIBDEFLATE
/* ----------------------------------------------------------------------
 * libdeflate compression code, with interface to match
 * zlib_mem_{in,de}flate for simplicity elsewhere.
 */

// Named the same as the version that uses zlib as we always use libdeflate for
// decompression when available.
char *zlib_mem_inflate(char *cdata, size_t csize, size_t *size) {
    struct libdeflate_decompressor *z = libdeflate_alloc_decompressor();
    if (!z) {
        hts_log_error("Call to libdeflate_alloc_decompressor failed");
        return NULL;
    }

    uint8_t *data = NULL, *new_data;
    if (!*size)
        *size = csize*2;
    for(;;) {
        new_data = realloc(data, *size);
        if (!new_data) {
            hts_log_error("Memory allocation failure");
            goto fail;
        }
        data = new_data;

        int ret = libdeflate_gzip_decompress(z, cdata, csize, data, *size, size);

        // Auto grow output buffer size if needed and try again.
        // Fortunately for all bar one call of this we know the size already.
        if (ret == LIBDEFLATE_INSUFFICIENT_SPACE) {
            (*size) *= 1.5;
            continue;
        }

        if (ret != LIBDEFLATE_SUCCESS) {
            hts_log_error("Inflate operation failed: %d", ret);
            goto fail;
        } else {
            break;
        }
    }

    libdeflate_free_decompressor(z);
    return (char *)data;

 fail:
    libdeflate_free_decompressor(z);
    free(data);
    return NULL;
}

// Named differently as we use both zlib/libdeflate for compression.
static char *libdeflate_deflate(char *data, size_t size, size_t *cdata_size,
                                int level, int strat) {
    level = level > 0 ? level : 6; // libdeflate doesn't honour -1 as default
    level *= 1.23;     // NB levels go up to 12 here; 5 onwards is +1
    level += level>=8; // 5,6,7->6,7,8  8->10  9->12
    if (level > 12) level = 12;

    if (strat == Z_RLE) // not supported by libdeflate
        level = 1;

    struct libdeflate_compressor *z = libdeflate_alloc_compressor(level);
    if (!z) {
        hts_log_error("Call to libdeflate_alloc_compressor failed");
        return NULL;
    }

    unsigned char *cdata = NULL; /* Compressed output */
    size_t cdata_alloc;
    cdata = malloc(cdata_alloc = size*1.05+100);
    if (!cdata) {
        hts_log_error("Memory allocation failure");
        libdeflate_free_compressor(z);
        return NULL;
    }

    *cdata_size = libdeflate_gzip_compress(z, data, size, cdata, cdata_alloc);
    libdeflate_free_compressor(z);

    if (*cdata_size == 0) {
        hts_log_error("Call to libdeflate_gzip_compress failed");
        free(cdata);
        return NULL;
    }

    return (char *)cdata;
}

#else

/* ----------------------------------------------------------------------
 * zlib compression code - from Gap5's tg_iface_g.c
 * They're static here as they're only used within the cram_compress_block
 * and cram_uncompress_block functions, which are the external interface.
 */
char *zlib_mem_inflate(char *cdata, size_t csize, size_t *size) {
    z_stream s;
    unsigned char *data = NULL; /* Uncompressed output */
    int data_alloc = 0;
    int err;

    /* Starting point at uncompressed size, and scale after that */
    data = malloc(data_alloc = csize*1.2+100);
    if (!data)
        return NULL;

    /* Initialise zlib stream */
    s.zalloc = Z_NULL; /* use default allocation functions */
    s.zfree  = Z_NULL;
    s.opaque = Z_NULL;
    s.next_in  = (unsigned char *)cdata;
    s.avail_in = csize;
    s.total_in = 0;
    s.next_out  = data;
    s.avail_out = data_alloc;
    s.total_out = 0;

    //err = inflateInit(&s);
    err = inflateInit2(&s, 15 + 32);
    if (err != Z_OK) {
        hts_log_error("Call to zlib inflateInit failed: %s", s.msg);
        free(data);
        return NULL;
    }

    /* Decode to 'data' array */
    for (;s.avail_in;) {
        unsigned char *data_tmp;
        int alloc_inc;

        s.next_out = &data[s.total_out];
        err = inflate(&s, Z_NO_FLUSH);
        if (err == Z_STREAM_END)
            break;

        if (err != Z_OK) {
            hts_log_error("Call to zlib inflate failed: %s", s.msg);
            free(data);
            inflateEnd(&s);
            return NULL;
        }

        /* More to come, so realloc based on growth so far */
        alloc_inc = (double)s.avail_in/s.total_in * s.total_out + 100;
        data = realloc((data_tmp = data), data_alloc += alloc_inc);
        if (!data) {
            free(data_tmp);
            inflateEnd(&s);
            return NULL;
        }
        s.avail_out += alloc_inc;
    }
    inflateEnd(&s);

    *size = s.total_out;
    return (char *)data;
}
#endif

/* ----------------------------------------------------------------------
 * CRAM blocks - the dynamically growable data block. We have code to
 * create, update, (un)compress and read/write.
 *
 * These are derived from the deflate_interlaced.c blocks, but with the
 * CRAM extension of content types and IDs.
 */

/*
 * Allocates a new cram_block structure with a specified content_type and
 * id.
 *
 * Returns block pointer on success
 *         NULL on failure
 */
cram_block *cram_new_block(enum cram_content_type content_type,
                           int content_id) {
    cram_block *b = malloc(sizeof(*b));
    if (!b)
        return NULL;
    b->method = b->orig_method = RAW;
    b->content_type = content_type;
    b->content_id = content_id;
    b->comp_size = 0;
    b->uncomp_size = 0;
    b->data = NULL;
    b->alloc = 0;
    b->byte = 0;
    b->bit = 7; // MSB
    b->crc32 = 0;
    b->idx = 0;
    b->m = NULL;

    return b;
}

/*
 * Reads a block from a cram file.
 * Returns cram_block pointer on success.
 *         NULL on failure
 */
cram_block *cram_read_block(cram_fd *fd) {
    cram_block *b = malloc(sizeof(*b));
    unsigned char c;
    uint32_t crc = 0;
    if (!b)
        return NULL;

    //fprintf(stderr, "Block at %d\n", (int)ftell(fd->fp));

    if (-1 == (b->method      = hgetc(fd->fp))) { free(b); return NULL; }
    c = b->method; crc = crc32(crc, &c, 1);
    if (-1 == (b->content_type= hgetc(fd->fp))) { free(b); return NULL; }
    c = b->content_type; crc = crc32(crc, &c, 1);
    if (-1 == fd->vv.varint_decode32_crc(fd, &b->content_id, &crc))  { free(b); return NULL; }
    if (-1 == fd->vv.varint_decode32_crc(fd, &b->comp_size, &crc))   { free(b); return NULL; }
    if (-1 == fd->vv.varint_decode32_crc(fd, &b->uncomp_size, &crc)) { free(b); return NULL; }

    //fprintf(stderr, "  method %d, ctype %d, cid %d, csize %d, ucsize %d\n",
    //      b->method, b->content_type, b->content_id, b->comp_size, b->uncomp_size);

    if (b->method == RAW) {
        if (b->uncomp_size < 0 || b->comp_size != b->uncomp_size) {
            free(b);
            return NULL;
        }
        b->alloc = b->uncomp_size;
        if (!(b->data = malloc(b->uncomp_size))){ free(b); return NULL; }
        if (b->uncomp_size != hread(fd->fp, b->data, b->uncomp_size)) {
            free(b->data);
            free(b);
            return NULL;
        }
    } else {
        if (b->comp_size < 0 || b->uncomp_size < 0) {
            free(b);
            return NULL;
        }
        b->alloc = b->comp_size;
        if (!(b->data = malloc(b->comp_size)))  { free(b); return NULL; }
        if (b->comp_size != hread(fd->fp, b->data, b->comp_size)) {
            free(b->data);
            free(b);
            return NULL;
        }
    }

    if (CRAM_MAJOR_VERS(fd->version) >= 3) {
        if (-1 == int32_decode(fd, (int32_t *)&b->crc32)) {
            free(b->data);
            free(b);
            return NULL;
        }

        b->crc32_checked = fd->ignore_md5;
        b->crc_part = crc;
    } else {
        b->crc32_checked = 1; // CRC not present
    }

    b->orig_method = b->method;
    b->idx = 0;
    b->byte = 0;
    b->bit = 7; // MSB

    return b;
}


/*
 * Writes a CRAM block.
 * Returns 0 on success
 *        -1 on failure
 */
int cram_write_block(cram_fd *fd, cram_block *b) {
    char vardata[100];
    int vardata_o = 0;

    assert(b->method != RAW || (b->comp_size == b->uncomp_size));

    if (hputc(b->method,       fd->fp)  == EOF) return -1;
    if (hputc(b->content_type, fd->fp)  == EOF) return -1;
    vardata_o += fd->vv.varint_put32(vardata          , vardata+100, b->content_id);
    vardata_o += fd->vv.varint_put32(vardata+vardata_o, vardata+100, b->comp_size);
    vardata_o += fd->vv.varint_put32(vardata+vardata_o, vardata+100, b->uncomp_size);
    if (vardata_o != hwrite(fd->fp, vardata, vardata_o))
        return -1;

    if (b->data) {
        if (b->method == RAW) {
            if (b->uncomp_size != hwrite(fd->fp, b->data, b->uncomp_size))
                return -1;
        } else {
            if (b->comp_size != hwrite(fd->fp, b->data, b->comp_size))
                return -1;
        }
    } else {
        // Absent blocks should be size 0
        assert(b->method == RAW && b->uncomp_size == 0);
    }

    if (CRAM_MAJOR_VERS(fd->version) >= 3) {
        char dat[100], *cp = (char *)dat;
        uint32_t crc;

        *cp++ = b->method;
        *cp++ = b->content_type;
        cp += fd->vv.varint_put32(cp, dat+100, b->content_id);
        cp += fd->vv.varint_put32(cp, dat+100, b->comp_size);
        cp += fd->vv.varint_put32(cp, dat+100, b->uncomp_size);
        crc = crc32(0L, (uc *)dat, cp-dat);

        if (b->method == RAW) {
            b->crc32 = crc32(crc, b->data ? b->data : (uc*)"", b->uncomp_size);
        } else {
            b->crc32 = crc32(crc, b->data ? b->data : (uc*)"", b->comp_size);
        }

        if (-1 == int32_encode(fd, b->crc32))
            return -1;
    }

    return 0;
}

/*
 * Frees a CRAM block, deallocating internal data too.
 */
void cram_free_block(cram_block *b) {
    if (!b)
        return;
    if (b->data)
        free(b->data);
    free(b);
}

cram_metrics *cram_new_metrics(void) {
    cram_metrics *m = calloc(1, sizeof(*m));
    if (!m)
        return NULL;
    m->trial = NTRIALS-1;
    m->next_trial = TRIAL_SPAN/2; // learn quicker at start
    m->method = RAW;
    m->strat = 0;
    m->revised_method = 0;
    m->unpackable = 0;

    return m;
}

char *cram_content_type2str(enum cram_content_type t) {
    switch (t) {
    case FILE_HEADER:         return "FILE_HEADER";
    case COMPRESSION_HEADER:  return "COMPRESSION_HEADER";
    case MAPPED_SLICE:        return "MAPPED_SLICE";
    case UNMAPPED_SLICE:      return "UNMAPPED_SLICE";
    case EXTERNAL:            return "EXTERNAL";
    case CORE:                return "CORE";
    case CT_ERROR:            break;
    }
    return "?";
}

/* ----------------------------------------------------------------------
 * Reference sequence handling
 *
 * These revolve around the refs_t structure, which may potentially be
 * shared between multiple cram_fd.
 *
 * We start with refs_create() to allocate an empty refs_t and then
 * populate it with @SQ line data using refs_from_header(). This is done on
 * cram_open().  Also at start up we can call cram_load_reference() which
 * is used with "scramble -r foo.fa". This replaces the fd->refs with the
 * new one specified. In either case refs2id() is then called which
 * maps ref_entry names to @SQ ids (refs_t->ref_id[]).
 *
 * Later, possibly within a thread, we will want to know the actual ref
 * seq itself, obtained by calling cram_get_ref().  This may use the
 * UR: or M5: fields or the filename specified in the original
 * cram_load_reference() call.
 *
 * Given the potential for multi-threaded reference usage, we have
 * reference counting (sorry for the confusing double use of "ref") to
 * track the number of callers interested in any specific reference.
 */

/*
 * Frees/unmaps a reference sequence and associated file handles.
 */
static void ref_entry_free_seq(ref_entry *e) {
    if (e->mf)
        mfclose(e->mf);
    if (e->seq && !e->mf)
        free(e->seq);

    e->seq = NULL;
    e->mf = NULL;
}

void refs_free(refs_t *r) {
    RP("refs_free()\n");

    if (--r->count > 0)
        return;

    if (!r)
        return;

    if (r->pool)
        string_pool_destroy(r->pool);

    if (r->h_meta) {
        khint_t k;

        for (k = kh_begin(r->h_meta); k != kh_end(r->h_meta); k++) {
            ref_entry *e;

            if (!kh_exist(r->h_meta, k))
                continue;
            if (!(e = kh_val(r->h_meta, k)))
                continue;
            ref_entry_free_seq(e);
            free(e);
        }

        kh_destroy(refs, r->h_meta);
    }

    if (r->ref_id)
        free(r->ref_id);

    if (r->fp)
        bgzf_close(r->fp);

    pthread_mutex_destroy(&r->lock);

    free(r);
}

static refs_t *refs_create(void) {
    refs_t *r = calloc(1, sizeof(*r));

    RP("refs_create()\n");

    if (!r)
        return NULL;

    if (!(r->pool = string_pool_create(8192)))
        goto err;

    r->ref_id = NULL; // see refs2id() to populate.
    r->count = 1;
    r->last = NULL;
    r->last_id = -1;

    if (!(r->h_meta = kh_init(refs)))
        goto err;

    pthread_mutex_init(&r->lock, NULL);

    return r;

 err:
    refs_free(r);
    return NULL;
}

/*
 * Opens a reference fasta file as a BGZF stream, allowing for
 * compressed files.  It automatically builds a .fai file if
 * required and if compressed a .gzi bgzf index too.
 *
 * Returns a BGZF handle on success;
 *         NULL on failure.
 */
static BGZF *bgzf_open_ref(char *fn, char *mode, int is_md5) {
    BGZF *fp;

    if (!is_md5 && !hisremote(fn)) {
        char fai_file[PATH_MAX];

        snprintf(fai_file, PATH_MAX, "%s.fai", fn);
        if (access(fai_file, R_OK) != 0)
            if (fai_build(fn) != 0)
                return NULL;
    }

    if (!(fp = bgzf_open(fn, mode))) {
        perror(fn);
        return NULL;
    }

    if (fp->is_compressed == 1 && bgzf_index_load(fp, fn, ".gzi") < 0) {
        hts_log_error("Unable to load .gzi index '%s.gzi'", fn);
        bgzf_close(fp);
        return NULL;
    }

    return fp;
}

/*
 * Loads a FAI file for a reference.fasta.
 * "is_err" indicates whether failure to load is worthy of emitting an
 * error message. In some cases (eg with embedded references) we
 * speculatively load, just in case, and silently ignore errors.
 *
 * Returns the refs_t struct on success (maybe newly allocated);
 *         NULL on failure
 */
static refs_t *refs_load_fai(refs_t *r_orig, const char *fn, int is_err) {
    hFILE *fp = NULL;
    char fai_fn[PATH_MAX];
    char line[8192];
    refs_t *r = r_orig;
    size_t fn_l = strlen(fn);
    int id = 0, id_alloc = 0;

    RP("refs_load_fai %s\n", fn);

    if (!r)
        if (!(r = refs_create()))
            goto err;

    if (r->fp)
        if (bgzf_close(r->fp) != 0)
            goto err;
    r->fp = NULL;

    /* Look for a FASTA##idx##FAI format */
    char *fn_delim = strstr(fn, HTS_IDX_DELIM);
    if (fn_delim) {
        if (!(r->fn = string_ndup(r->pool, fn, fn_delim - fn)))
            goto err;
        fn_delim += strlen(HTS_IDX_DELIM);
        snprintf(fai_fn, PATH_MAX, "%s", fn_delim);
    } else {
        /* An index file was provided, instead of the actual reference file */
        if (fn_l > 4 && strcmp(&fn[fn_l-4], ".fai") == 0) {
            if (!r->fn) {
                if (!(r->fn = string_ndup(r->pool, fn, fn_l-4)))
                    goto err;
            }
            snprintf(fai_fn, PATH_MAX, "%s", fn);
        } else {
        /* Only the reference file provided. Get the index file name from it */
            if (!(r->fn = string_dup(r->pool, fn)))
                goto err;
            snprintf(fai_fn, PATH_MAX, "%.*s.fai", PATH_MAX-5, fn);
        }
    }

    if (!(r->fp = bgzf_open_ref(r->fn, "r", 0))) {
        hts_log_error("Failed to open reference file '%s'", r->fn);
        goto err;
    }

    if (!(fp = hopen(fai_fn, "r"))) {
        hts_log_error("Failed to open index file '%s'", fai_fn);
        if (is_err)
            perror(fai_fn);
        goto err;
    }
    while (hgets(line, 8192, fp) != NULL) {
        ref_entry *e = malloc(sizeof(*e));
        char *cp;
        int n;
        khint_t k;

        if (!e)
            return NULL;

        // id
        for (cp = line; *cp && !isspace_c(*cp); cp++)
            ;
        *cp++ = 0;
        e->name = string_dup(r->pool, line);

        // length
        while (*cp && isspace_c(*cp))
            cp++;
        e->length = strtoll(cp, &cp, 10);

        // offset
        while (*cp && isspace_c(*cp))
            cp++;
        e->offset = strtoll(cp, &cp, 10);

        // bases per line
        while (*cp && isspace_c(*cp))
            cp++;
        e->bases_per_line = strtol(cp, &cp, 10);

        // line length
        while (*cp && isspace_c(*cp))
            cp++;
        e->line_length = strtol(cp, &cp, 10);

        // filename
        e->fn = r->fn;

        e->count = 0;
        e->seq = NULL;
        e->mf = NULL;
        e->is_md5 = 0;
        e->validated_md5 = 0;

        k = kh_put(refs, r->h_meta, e->name, &n);
        if (-1 == n)  {
            free(e);
            return NULL;
        }

        if (n) {
            kh_val(r->h_meta, k) = e;
        } else {
            ref_entry *re = kh_val(r->h_meta, k);
            if (re && (re->count != 0 || re->length != 0)) {
                /* Keep old */
                free(e);
            } else {
                /* Replace old */
                if (re)
                    free(re);
                kh_val(r->h_meta, k) = e;
            }
        }

        if (id >= id_alloc) {
            ref_entry **new_refs;
            int x;

            id_alloc = id_alloc ?id_alloc*2 : 16;
            new_refs = realloc(r->ref_id, id_alloc * sizeof(*r->ref_id));
            if (!new_refs)
                goto err;
            r->ref_id = new_refs;

            for (x = id; x < id_alloc; x++)
                r->ref_id[x] = NULL;
        }
        r->ref_id[id] = e;
        r->nref = ++id;
    }

    if(hclose(fp) < 0)
        goto err;
    return r;

 err:
    if (fp)
        hclose_abruptly(fp);

    if (!r_orig)
        refs_free(r);

    return NULL;
}

/*
 * Verifies that the CRAM @SQ lines and .fai files match.
 */
static void sanitise_SQ_lines(cram_fd *fd) {
    int i;

    if (!fd->header || !fd->header->hrecs)
        return;

    if (!fd->refs || !fd->refs->h_meta)
        return;

    for (i = 0; i < fd->header->hrecs->nref; i++) {
        const char *name = fd->header->hrecs->ref[i].name;
        khint_t k = kh_get(refs, fd->refs->h_meta, name);
        ref_entry *r;

        // We may have @SQ lines which have no known .fai, but do not
        // in themselves pose a problem because they are unused in the file.
        if (k == kh_end(fd->refs->h_meta))
            continue;

        if (!(r = (ref_entry *)kh_val(fd->refs->h_meta, k)))
            continue;

        if (r->length && r->length != fd->header->hrecs->ref[i].len) {
            assert(strcmp(r->name, fd->header->hrecs->ref[i].name) == 0);

            // Should we also check MD5sums here to ensure the correct
            // reference was given?
            hts_log_warning("Header @SQ length mismatch for ref %s, %"PRIhts_pos" vs %d",
                            r->name, fd->header->hrecs->ref[i].len, (int)r->length);

            // Fixing the parsed @SQ header will make MD:Z: strings work
            // and also stop it producing N for the sequence.
            fd->header->hrecs->ref[i].len = r->length;
        }
    }
}

/*
 * Indexes references by the order they appear in a BAM file. This may not
 * necessarily be the same order they appear in the fasta reference file.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int refs2id(refs_t *r, sam_hdr_t *hdr) {
    int i;
    sam_hrecs_t *h = hdr->hrecs;

    if (r->ref_id)
        free(r->ref_id);
    if (r->last)
        r->last = NULL;

    r->ref_id = calloc(h->nref, sizeof(*r->ref_id));
    if (!r->ref_id)
        return -1;

    r->nref = h->nref;
    for (i = 0; i < h->nref; i++) {
        khint_t k = kh_get(refs, r->h_meta, h->ref[i].name);
        if (k != kh_end(r->h_meta)) {
            r->ref_id[i] = kh_val(r->h_meta, k);
        } else {
            hts_log_warning("Unable to find ref name '%s'", h->ref[i].name);
        }
    }

    return 0;
}

/*
 * Generates refs_t entries based on @SQ lines in the header.
 * Returns 0 on success
 *         -1 on failure
 */
static int refs_from_header(cram_fd *fd) {
    if (!fd)
        return -1;

    refs_t *r = fd->refs;
    if (!r)
        return -1;

    sam_hdr_t *h = fd->header;
    if (!h)
        return 0;

    if (!h->hrecs) {
        if (-1 == sam_hdr_fill_hrecs(h))
            return -1;
    }

    if (h->hrecs->nref == 0)
        return 0;

    //fprintf(stderr, "refs_from_header for %p mode %c\n", fd, fd->mode);

    /* Existing refs are fine, as long as they're compatible with the hdr. */
    ref_entry **new_ref_id = realloc(r->ref_id, (r->nref + h->hrecs->nref) * sizeof(*r->ref_id));
    if (!new_ref_id)
        return -1;
    r->ref_id = new_ref_id;

    int i, j;
    /* Copy info from h->ref[i] over to r */
    for (i = 0, j = r->nref; i < h->hrecs->nref; i++) {
        sam_hrec_type_t *ty;
        sam_hrec_tag_t *tag;
        khint_t k;
        int n;

        k = kh_get(refs, r->h_meta, h->hrecs->ref[i].name);
        if (k != kh_end(r->h_meta))
            // Ref already known about
            continue;

        if (!(r->ref_id[j] = calloc(1, sizeof(ref_entry))))
            return -1;

        if (!h->hrecs->ref[i].name)
            return -1;

        r->ref_id[j]->name = string_dup(r->pool, h->hrecs->ref[i].name);
        if (!r->ref_id[j]->name) return -1;
        r->ref_id[j]->length = 0; // marker for not yet loaded

        /* Initialise likely filename if known */
        if ((ty = sam_hrecs_find_type_id(h->hrecs, "SQ", "SN", h->hrecs->ref[i].name))) {
            if ((tag = sam_hrecs_find_key(ty, "M5", NULL))) {
                r->ref_id[j]->fn = string_dup(r->pool, tag->str+3);
                //fprintf(stderr, "Tagging @SQ %s / %s\n", r->ref_id[h]->name, r->ref_id[h]->fn);
            }
        }

        k = kh_put(refs, r->h_meta, r->ref_id[j]->name, &n);
        if (n <= 0) // already exists or error
            return -1;
        kh_val(r->h_meta, k) = r->ref_id[j];

        j++;
    }
    r->nref = j;

    return 0;
}

/*
 * Attaches a header to a cram_fd.
 *
 * This should be used when creating a new cram_fd for writing where
 * we have a header already constructed (eg from a file we've read
 * in).
 */
int cram_set_header2(cram_fd *fd, const sam_hdr_t *hdr) {
    if (!fd || !hdr )
        return -1;

    if (fd->header != hdr) {
        if (fd->header)
            sam_hdr_destroy(fd->header);
        fd->header = sam_hdr_dup(hdr);
        if (!fd->header)
            return -1;
    }
    return refs_from_header(fd);
}

int cram_set_header(cram_fd *fd, sam_hdr_t *hdr) {
    return cram_set_header2(fd, hdr);
}

/*
 * Returns whether the path refers to a directory.
 */
static int is_directory(char *fn) {
    struct stat buf;
    if ( stat(fn,&buf) ) return 0;
    return S_ISDIR(buf.st_mode);
}

/*
 * Converts a directory and a filename into an expanded path, replacing %s
 * in directory with the filename and %[0-9]+s with portions of the filename
 * Any remaining parts of filename are added to the end with /%s.
 */
static int expand_cache_path(char *path, char *dir, const char *fn) {
    char *cp, *start = path;
    size_t len;
    size_t sz = PATH_MAX;

    while ((cp = strchr(dir, '%'))) {
        if (cp-dir >= sz) return -1;
        strncpy(path, dir, cp-dir);
        path += cp-dir;
        sz -= cp-dir;

        if (*++cp == 's') {
            len = strlen(fn);
            if (len >= sz) return -1;
            strcpy(path, fn);
            path += len;
            sz -= len;
            fn += len;
            cp++;
        } else if (*cp >= '0' && *cp <= '9') {
            char *endp;
            long l;

            l = strtol(cp, &endp, 10);
            l = MIN(l, strlen(fn));
            if (*endp == 's') {
                if (l >= sz) return -1;
                strncpy(path, fn, l);
                path += l;
                fn += l;
                sz -= l;
                *path = 0;
                cp = endp+1;
            } else {
                if (sz < 3) return -1;
                *path++ = '%';
                *path++ = *cp++;
            }
        } else {
            if (sz < 3) return -1;
            *path++ = '%';
            *path++ = *cp++;
        }
        dir = cp;
    }

    len = strlen(dir);
    if (len >= sz) return -1;
    strcpy(path, dir);
    path += len;
    sz -= len;

    len = strlen(fn) + ((*fn && path > start && path[-1] != '/') ? 1 : 0);
    if (len >= sz) return -1;
    if (*fn && path > start && path[-1] != '/')
        *path++ = '/';
    strcpy(path, fn);
    return 0;
}

/*
 * Make the directory containing path and any prefix directories.
 */
static void mkdir_prefix(char *path, int mode) {
    char *cp = strrchr(path, '/');
    if (!cp)
        return;

    *cp = 0;
    if (is_directory(path)) {
        *cp = '/';
        return;
    }

    if (mkdir(path, mode) == 0) {
        chmod(path, mode);
        *cp = '/';
        return;
    }

    mkdir_prefix(path, mode);
    mkdir(path, mode);
    chmod(path, mode);
    *cp = '/';
}

/*
 * Return the cache directory to use, based on the first of these
 * environment variables to be set to a non-empty value.
 */
static const char *get_cache_basedir(const char **extra) {
    char *base;

    *extra = "";

    base = getenv("XDG_CACHE_HOME");
    if (base && *base) return base;

    base = getenv("HOME");
    if (base && *base) { *extra = "/.cache"; return base; }

    base = getenv("TMPDIR");
    if (base && *base) return base;

    base = getenv("TEMP");
    if (base && *base) return base;

    return "/tmp";
}

/*
 * Queries the M5 string from the header and attempts to populate the
 * reference from this using the REF_PATH environment.
 *
 * Returns 0 on success
 *        -1 on failure
 */
static int cram_populate_ref(cram_fd *fd, int id, ref_entry *r) {
    char *ref_path = getenv("REF_PATH");
    sam_hrec_type_t *ty;
    sam_hrec_tag_t *tag;
    char path[PATH_MAX];
    kstring_t path_tmp = KS_INITIALIZE;
    char cache[PATH_MAX], cache_root[PATH_MAX];
    char *local_cache = getenv("REF_CACHE");
    mFILE *mf;
    int local_path = 0;

    hts_log_info("Running cram_populate_ref on fd %p, id %d", (void *)fd, id);

    cache_root[0] = '\0';

    if (!ref_path || *ref_path == '\0') {
        /*
         * If we have no ref path, we use the EBI server.
         * However to avoid spamming it we require a local ref cache too.
         */
        ref_path = "https://www.ebi.ac.uk/ena/cram/md5/%s";
        if (!local_cache || *local_cache == '\0') {
            const char *extra;
            const char *base = get_cache_basedir(&extra);
            snprintf(cache_root, PATH_MAX, "%s%s/hts-ref", base, extra);
            snprintf(cache,PATH_MAX, "%s%s/hts-ref/%%2s/%%2s/%%s", base, extra);
            local_cache = cache;
            hts_log_info("Populating local cache: %s", local_cache);
        }
    }

    if (!r->name)
        return -1;

    if (!(ty = sam_hrecs_find_type_id(fd->header->hrecs, "SQ", "SN", r->name)))
        return -1;

    if (!(tag = sam_hrecs_find_key(ty, "M5", NULL)))
        goto no_M5;

    hts_log_info("Querying ref %s", tag->str+3);

    /* Use cache if available */
    if (local_cache && *local_cache) {
        if (expand_cache_path(path, local_cache, tag->str+3) == 0)
            local_path = 1;
    }

#ifndef HAVE_MMAP
    char *path2;
    /* Search local files in REF_PATH; we can open them and return as above */
    if (!local_path && (path2 = find_path(tag->str+3, ref_path))) {
        int len = snprintf(path, PATH_MAX, "%s", path2);
        free(path2);
        if (len > 0 && len < PATH_MAX) // in case it's too long
            local_path = 1;
    }
#endif

    /* Found via REF_CACHE or local REF_PATH file */
    if (local_path) {
        struct stat sb;
        BGZF *fp;

        if (0 == stat(path, &sb)
            && S_ISREG(sb.st_mode)
            && (fp = bgzf_open(path, "r"))) {
            r->length = sb.st_size;
            r->offset = r->line_length = r->bases_per_line = 0;

            r->fn = string_dup(fd->refs->pool, path);

            if (fd->refs->fp)
                if (bgzf_close(fd->refs->fp) != 0)
                    return -1;
            fd->refs->fp = fp;
            fd->refs->fn = r->fn;
            r->is_md5 = 1;
            r->validated_md5 = 1;

            // Fall back to cram_get_ref() where it'll do the actual
            // reading of the file.
            return 0;
        }
    }


    /* Otherwise search full REF_PATH; slower as loads entire file */
    if ((mf = open_path_mfile(tag->str+3, ref_path, NULL))) {
        size_t sz;
        r->seq = mfsteal(mf, &sz);
        if (r->seq) {
            r->mf = NULL;
        } else {
            // keep mf around as we couldn't detach
            r->seq = mf->data;
            r->mf = mf;
        }
        r->length = sz;
        r->is_md5 = 1;
        r->validated_md5 = 1;
    } else {
        refs_t *refs;
        const char *fn;

    no_M5:
        /* Failed to find in search path or M5 cache, see if @SQ UR: tag? */
        if (!(tag = sam_hrecs_find_key(ty, "UR", NULL)))
            return -1;

        fn = (strncmp(tag->str+3, "file:", 5) == 0)
            ? tag->str+8
            : tag->str+3;

        if (fd->refs->fp) {
            if (bgzf_close(fd->refs->fp) != 0)
                return -1;
            fd->refs->fp = NULL;
        }
        if (!(refs = refs_load_fai(fd->refs, fn, 0)))
            return -1;
        sanitise_SQ_lines(fd);

        fd->refs = refs;
        if (fd->refs->fp) {
            if (bgzf_close(fd->refs->fp) != 0)
                return -1;
            fd->refs->fp = NULL;
        }

        if (!fd->refs->fn)
            return -1;

        if (-1 == refs2id(fd->refs, fd->header))
            return -1;
        if (!fd->refs->ref_id || !fd->refs->ref_id[id])
            return -1;

        // Local copy already, so fall back to cram_get_ref().
        return 0;
    }

    /* Populate the local disk cache if required */
    if (local_cache && *local_cache) {
        hFILE *fp;

        if (*cache_root && !is_directory(cache_root)) {
            hts_log_warning("Creating reference cache directory %s\n"
                            "This may become large; see the samtools(1) manual page REF_CACHE discussion",
                            cache_root);
        }

        if (expand_cache_path(path, local_cache, tag->str+3) < 0) {
            return 0; // Not fatal - we have the data already so keep going.
        }
        hts_log_info("Writing cache file '%s'", path);
        mkdir_prefix(path, 01777);

        fp = hts_open_tmpfile(path, "wx", &path_tmp);
        if (!fp) {
            perror(path_tmp.s);
            free(path_tmp.s);

            // Not fatal - we have the data already so keep going.
            return 0;
        }

        // Check md5sum
        hts_md5_context *md5;
        char unsigned md5_buf1[16];
        char md5_buf2[33];

        if (!(md5 = hts_md5_init())) {
            hclose_abruptly(fp);
            unlink(path_tmp.s);
            free(path_tmp.s);
            return -1;
        }
        hts_md5_update(md5, r->seq, r->length);
        hts_md5_final(md5_buf1, md5);
        hts_md5_destroy(md5);
        hts_md5_hex(md5_buf2, md5_buf1);

        if (strncmp(tag->str+3, md5_buf2, 32) != 0) {
            hts_log_error("Mismatching md5sum for downloaded reference");
            hclose_abruptly(fp);
            unlink(path_tmp.s);
            free(path_tmp.s);
            return -1;
        }

        ssize_t length_written = hwrite(fp, r->seq, r->length);
        if (hclose(fp) < 0 || length_written != r->length ||
            chmod(path_tmp.s, 0444) < 0 ||
            rename(path_tmp.s, path) < 0) {
            hts_log_error("Creating reference at %s failed: %s",
                          path, strerror(errno));
            unlink(path_tmp.s);
        }
    }

    free(path_tmp.s);
    return 0;
}

static void cram_ref_incr_locked(refs_t *r, int id) {
    RP("%d INC REF %d, %d %p\n", gettid(), id,
       (int)(id>=0 && r->ref_id[id]?r->ref_id[id]->count+1:-999),
       id>=0 && r->ref_id[id]?r->ref_id[id]->seq:(char *)1);

    if (id < 0 || !r->ref_id[id] || !r->ref_id[id]->seq)
        return;

    if (r->last_id == id)
        r->last_id = -1;

    ++r->ref_id[id]->count;
}

void cram_ref_incr(refs_t *r, int id) {
    pthread_mutex_lock(&r->lock);
    cram_ref_incr_locked(r, id);
    pthread_mutex_unlock(&r->lock);
}

static void cram_ref_decr_locked(refs_t *r, int id) {
    RP("%d DEC REF %d, %d %p\n", gettid(), id,
       (int)(id>=0 && r->ref_id[id]?r->ref_id[id]->count-1:-999),
       id>=0 && r->ref_id[id]?r->ref_id[id]->seq:(char *)1);

    if (id < 0 || !r->ref_id[id] || !r->ref_id[id]->seq) {
        return;
    }

    if (--r->ref_id[id]->count <= 0) {
        assert(r->ref_id[id]->count == 0);
        if (r->last_id >= 0) {
            if (r->ref_id[r->last_id]->count <= 0 &&
                r->ref_id[r->last_id]->seq) {
                RP("%d FREE REF %d (%p)\n", gettid(),
                   r->last_id, r->ref_id[r->last_id]->seq);
                ref_entry_free_seq(r->ref_id[r->last_id]);
                if (r->ref_id[r->last_id]->is_md5) r->ref_id[r->last_id]->length = 0;
            }
        }
        r->last_id = id;
    }
}

void cram_ref_decr(refs_t *r, int id) {
    pthread_mutex_lock(&r->lock);
    cram_ref_decr_locked(r, id);
    pthread_mutex_unlock(&r->lock);
}

/*
 * Used by cram_ref_load and cram_get_ref. The file handle will have
 * already been opened, so we can catch it. The ref_entry *e informs us
 * of whether this is a multi-line fasta file or a raw MD5 style file.
 * Either way we create a single contiguous sequence.
 *
 * Returns all or part of a reference sequence on success (malloced);
 *         NULL on failure.
 */
static char *load_ref_portion(BGZF *fp, ref_entry *e, int start, int end) {
    off_t offset, len;
    char *seq;

    if (end < start)
        end = start;

    /*
     * Compute locations in file. This is trivial for the MD5 files, but
     * is still necessary for the fasta variants.
     *
     * Note the offset here, as with faidx, has the assumption that white-
     * space (the diff between line_length and bases_per_line) only occurs
     * at the end of a line of text.
     */
    offset = e->line_length
        ? e->offset + (start-1)/e->bases_per_line * e->line_length +
          (start-1) % e->bases_per_line
        : start-1;

    len = (e->line_length
           ? e->offset + (end-1)/e->bases_per_line * e->line_length +
             (end-1) % e->bases_per_line
           : end-1) - offset + 1;

    if (bgzf_useek(fp, offset, SEEK_SET) < 0) {
        perror("bgzf_useek() on reference file");
        return NULL;
    }

    if (len == 0 || !(seq = malloc(len))) {
        return NULL;
    }

    if (len != bgzf_read(fp, seq, len)) {
        perror("bgzf_read() on reference file");
        free(seq);
        return NULL;
    }

    /* Strip white-space if required. */
    if (len != end-start+1) {
        hts_pos_t i, j;
        char *cp = seq;
        char *cp_to;

        // Copy up to the first white-space, and then repeatedly just copy
        // bases_per_line verbatim, and use the slow method to end again.
        //
        // This may seem excessive, but this code can be a significant
        // portion of total CRAM decode CPU time for shallow data sets.
        for (i = j = 0; i < len; i++) {
            if (!isspace_c(cp[i]))
                cp[j++] = cp[i] & ~0x20;
            else
                break;
        }
        while (i < len && isspace_c(cp[i]))
            i++;
        while (i < len - e->line_length) {
            hts_pos_t j_end = j + e->bases_per_line;
            while (j < j_end)
                cp[j++] = cp[i++] & ~0x20; // toupper equiv
            i += e->line_length - e->bases_per_line;
        }
        for (; i < len; i++) {
            if (!isspace_c(cp[i]))
                cp[j++] = cp[i] & ~0x20;
        }

        cp_to = cp+j;

        if (cp_to - seq != end-start+1) {
            hts_log_error("Malformed reference file");
            free(seq);
            return NULL;
        }
    } else {
        int i;
        for (i = 0; i < len; i++) {
            seq[i] = toupper_c(seq[i]);
        }
    }

    return seq;
}

/*
 * Load the entire reference 'id'.
 * This also increments the reference count by 1.
 *
 * Returns ref_entry on success;
 *         NULL on failure
 */
ref_entry *cram_ref_load(refs_t *r, int id, int is_md5) {
    ref_entry *e = r->ref_id[id];
    int start = 1, end = e->length;
    char *seq;

    if (e->seq) {
        return e;
    }

    assert(e->count == 0);

    if (r->last) {
#ifdef REF_DEBUG
        int idx = 0;
        for (idx = 0; idx < r->nref; idx++)
            if (r->last == r->ref_id[idx])
                break;
        RP("%d cram_ref_load DECR %d\n", gettid(), idx);
#endif
        assert(r->last->count > 0);
        if (--r->last->count <= 0) {
            RP("%d FREE REF %d (%p)\n", gettid(), id, r->ref_id[id]->seq);
            if (r->last->seq)
                ref_entry_free_seq(r->last);
        }
    }

    if (!r->fn)
        return NULL;

    /* Open file if it's not already the current open reference */
    if (strcmp(r->fn, e->fn) || r->fp == NULL) {
        if (r->fp)
            if (bgzf_close(r->fp) != 0)
                return NULL;
        r->fn = e->fn;
        if (!(r->fp = bgzf_open_ref(r->fn, "r", is_md5)))
            return NULL;
    }

    RP("%d Loading ref %d (%d..%d)\n", gettid(), id, start, end);

    if (!(seq = load_ref_portion(r->fp, e, start, end))) {
        return NULL;
    }

    RP("%d Loaded ref %d (%d..%d) = %p\n", gettid(), id, start, end, seq);

    RP("%d INC REF %d, %"PRId64"\n", gettid(), id, (e->count+1));
    e->seq = seq;
    e->mf = NULL;
    e->count++;

    /*
     * Also keep track of last used ref so incr/decr loops on the same
     * sequence don't cause load/free loops.
     */
    RP("%d cram_ref_load INCR %d => %"PRId64"\n", gettid(), id, e->count+1);
    r->last = e;
    e->count++;

    return e;
}

/*
 * Returns a portion of a reference sequence from start to end inclusive.
 * The returned pointer is owned by either the cram_file fd or by the
 * internal refs_t structure and should not be freed  by the caller.
 *
 * The difference is whether or not this refs_t is in use by just the one
 * cram_fd or by multiples, or whether we have multiple threads accessing
 * references. In either case fd->shared will be true and we start using
 * reference counting to track the number of users of a specific reference
 * sequence.
 *
 * Otherwise the ref seq returned is allocated as part of cram_fd itself
 * and will be freed up on the next call to cram_get_ref or cram_close.
 *
 * To return the entire reference sequence, specify start as 1 and end
 * as 0.
 *
 * To cease using a reference, call cram_ref_decr().
 *
 * Returns reference on success,
 *         NULL on failure
 */
char *cram_get_ref(cram_fd *fd, int id, int start, int end) {
    ref_entry *r;
    char *seq;
    int ostart = start;

    if (id == -1 || start < 1)
        return NULL;

    /* FIXME: axiomatic query of r->seq being true?
     * Or shortcut for unsorted data where we load once and never free?
     */

    //fd->shared_ref = 1; // hard code for now to simplify things

    pthread_mutex_lock(&fd->ref_lock);

    RP("%d cram_get_ref on fd %p, id %d, range %d..%d\n", gettid(), fd, id, start, end);

    /*
     * Unsorted data implies we want to fetch an entire reference at a time.
     * We just deal with this at the moment by claiming we're sharing
     * references instead, which has the same requirement.
     */
    if (fd->unsorted)
        fd->shared_ref = 1;


    /* Sanity checking: does this ID exist? */
    if (id >= fd->refs->nref) {
        hts_log_error("No reference found for id %d", id);
        pthread_mutex_unlock(&fd->ref_lock);
        return NULL;
    }

    if (!fd->refs || !fd->refs->ref_id[id]) {
        hts_log_error("No reference found for id %d", id);
        pthread_mutex_unlock(&fd->ref_lock);
        return NULL;
    }

    if (!(r = fd->refs->ref_id[id])) {
        hts_log_error("No reference found for id %d", id);
        pthread_mutex_unlock(&fd->ref_lock);
        return NULL;
    }


    /*
     * It has an entry, but may not have been populated yet.
     * Any manually loaded .fai files have their lengths known.
     * A ref entry computed from @SQ lines (M5 or UR field) will have
     * r->length == 0 unless it's been loaded once and verified that we have
     * an on-disk filename for it.
     *
     * 19 Sep 2013: Moved the lock here as the cram_populate_ref code calls
     * open_path_mfile and libcurl, which isn't multi-thread safe unless I
     * rewrite my code to have one curl handle per thread.
     */
    pthread_mutex_lock(&fd->refs->lock);
    if (r->length == 0) {
        if (fd->ref_fn)
            hts_log_warning("Reference file given, but ref '%s' not present",
                            r->name);
        if (cram_populate_ref(fd, id, r) == -1) {
            hts_log_warning("Failed to populate reference for id %d", id);
            pthread_mutex_unlock(&fd->refs->lock);
            pthread_mutex_unlock(&fd->ref_lock);
            return NULL;
        }
        r = fd->refs->ref_id[id];
        if (fd->unsorted)
            cram_ref_incr_locked(fd->refs, id);
    }


    /*
     * We now know that we the filename containing the reference, so check
     * for limits. If it's over half the reference we'll load all of it in
     * memory as this will speed up subsequent calls.
     */
    if (end < 1)
        end = r->length;
    if (end >= r->length)
        end  = r->length;

    if (end - start >= 0.5*r->length || fd->shared_ref) {
        start = 1;
        end = r->length;
    }

    /*
     * Maybe we have it cached already? If so use it.
     *
     * Alternatively if we don't have the sequence but we're sharing
     * references and/or are asking for the entire length of it, then
     * load the full reference into the refs structure and return
     * a pointer to that one instead.
     */
    if (fd->shared_ref || r->seq || (start == 1 && end == r->length)) {
        char *cp;

        if (id >= 0) {
            if (r->seq) {
                cram_ref_incr_locked(fd->refs, id);
            } else {
                ref_entry *e;
                if (!(e = cram_ref_load(fd->refs, id, r->is_md5))) {
                    pthread_mutex_unlock(&fd->refs->lock);
                    pthread_mutex_unlock(&fd->ref_lock);
                    return NULL;
                }

                /* unsorted data implies cache ref indefinitely, to avoid
                 * continually loading and unloading.
                 */
                if (fd->unsorted)
                    cram_ref_incr_locked(fd->refs, id);
            }

            fd->ref = NULL; /* We never access it directly */
            fd->ref_start = 1;
            fd->ref_end   = r->length;
            fd->ref_id    = id;

            cp = fd->refs->ref_id[id]->seq + ostart-1;
        } else {
            fd->ref = NULL;
            cp = NULL;
        }

        RP("%d cram_get_ref returning for id %d, count %d\n", gettid(), id, (int)r->count);

        pthread_mutex_unlock(&fd->refs->lock);
        pthread_mutex_unlock(&fd->ref_lock);
        return cp;
    }

    /*
     * Otherwise we're not sharing, we don't have a copy of it already and
     * we're only asking for a small portion of it.
     *
     * In this case load up just that segment ourselves, freeing any old
     * small segments in the process.
     */

    /* Unmapped ref ID */
    if (id < 0 || !fd->refs->fn) {
        if (fd->ref_free) {
            free(fd->ref_free);
            fd->ref_free = NULL;
        }
        fd->ref = NULL;
        fd->ref_id = id;
        pthread_mutex_unlock(&fd->refs->lock);
        pthread_mutex_unlock(&fd->ref_lock);
        return NULL;
    }

    /* Open file if it's not already the current open reference */
    if (strcmp(fd->refs->fn, r->fn) || fd->refs->fp == NULL) {
        if (fd->refs->fp)
            if (bgzf_close(fd->refs->fp) != 0)
                return NULL;
        fd->refs->fn = r->fn;
        if (!(fd->refs->fp = bgzf_open_ref(fd->refs->fn, "r", r->is_md5))) {
            pthread_mutex_unlock(&fd->refs->lock);
            pthread_mutex_unlock(&fd->ref_lock);
            return NULL;
        }
    }

    if (!(fd->ref = load_ref_portion(fd->refs->fp, r, start, end))) {
        pthread_mutex_unlock(&fd->refs->lock);
        pthread_mutex_unlock(&fd->ref_lock);
        return NULL;
    }

    if (fd->ref_free)
        free(fd->ref_free);

    fd->ref_id    = id;
    fd->ref_start = start;
    fd->ref_end   = end;
    fd->ref_free = fd->ref;
    seq = fd->ref;

    pthread_mutex_unlock(&fd->refs->lock);
    pthread_mutex_unlock(&fd->ref_lock);

    return seq ? seq + ostart - start : NULL;
}

/*
 * If fd has been opened for reading, it may be permitted to specify 'fn'
 * as NULL and let the code auto-detect the reference by parsing the
 * SAM header @SQ lines.
 */
int cram_load_reference(cram_fd *fd, char *fn) {
    int ret = 0;

    if (fn) {
        fd->refs = refs_load_fai(fd->refs, fn,
                                 !(fd->embed_ref>0 && fd->mode == 'r'));
        fn = fd->refs ? fd->refs->fn : NULL;
        if (!fn)
            ret = -1;
        sanitise_SQ_lines(fd);
    }
    fd->ref_fn = fn;

    if ((!fd->refs || (fd->refs->nref == 0 && !fn)) && fd->header) {
        if (fd->refs)
            refs_free(fd->refs);
        if (!(fd->refs = refs_create()))
            return -1;
        if (-1 == refs_from_header(fd))
            return -1;
    }

    if (fd->header)
        if (-1 == refs2id(fd->refs, fd->header))
            return -1;

    return ret;
}

void cram_free_container(cram_container *c) {
    enum cram_DS_ID id;
    int i;

    if (!c)
        return;

    if (c->refs_used)
        free(c->refs_used);

    if (c->landmark)
        free(c->landmark);

    if (c->comp_hdr)
        cram_free_compression_header(c->comp_hdr);

    if (c->comp_hdr_block)
        cram_free_block(c->comp_hdr_block);

    // Free the slices; filled out by encoder only
    if (c->slices) {
        for (i = 0; i < c->max_slice; i++) {
            if (c->slices[i])
                cram_free_slice(c->slices[i]);
            if (c->slices[i] == c->slice)
                c->slice = NULL;
        }
        free(c->slices);
    }

    // Free the current slice; set by both encoder & decoder
    if (c->slice) {
        cram_free_slice(c->slice);
        c->slice = NULL;
    }

    for (id = DS_RN; id < DS_TN; id++)
        if (c->stats[id]) cram_stats_free(c->stats[id]);

    //if (c->aux_B_stats) cram_stats_free(c->aux_B_stats);

    if (c->tags_used) {
        khint_t k;

        for (k = kh_begin(c->tags_used); k != kh_end(c->tags_used); k++) {
            if (!kh_exist(c->tags_used, k))
                continue;

            cram_tag_map *tm = (cram_tag_map *)kh_val(c->tags_used, k);
            if (tm) {
                cram_codec *c = tm->codec;

                if (c) c->free(c);
                free(tm);
            }
        }

        kh_destroy(m_tagmap, c->tags_used);
    }

    if (c->ref_free)
        free(c->ref);

    free(c);
}

/*
 * Reads a container header.
 *
 * Returns cram_container on success
 *         NULL on failure or no container left (fd->err == 0).
 */
cram_container *cram_read_container(cram_fd *fd) {
    cram_container c2, *c;
    int i, s;
    size_t rd = 0;
    uint32_t crc = 0;

    fd->err = 0;
    fd->eof = 0;

    memset(&c2, 0, sizeof(c2));
    if (CRAM_MAJOR_VERS(fd->version) == 1) {
        if ((s = fd->vv.varint_decode32_crc(fd, &c2.length, &crc)) == -1) {
            fd->eof = fd->empty_container ? 1 : 2;
            return NULL;
        } else {
            rd+=s;
        }
    } else if (CRAM_MAJOR_VERS(fd->version) < 4) {
        uint32_t len;
        if ((s = int32_decode(fd, &c2.length)) == -1) {
            if (CRAM_MAJOR_VERS(fd->version) == 2 &&
                CRAM_MINOR_VERS(fd->version) == 0)
                fd->eof = 1; // EOF blocks arrived in v2.1
            else
                fd->eof = fd->empty_container ? 1 : 2;
            return NULL;
        } else {
            rd+=s;
        }
        len = le_int4(c2.length);
        crc = crc32(0L, (unsigned char *)&len, 4);
    } else {
        if ((s = fd->vv.varint_decode32_crc(fd, &c2.length, &crc))   == -1) {
            fd->eof = fd->empty_container ? 1 : 2;
            return NULL;
        } else {
            rd+=s;
        }
    }
    if ((s = fd->vv.varint_decode32s_crc(fd, &c2.ref_seq_id, &crc))   == -1) return NULL; else rd+=s;
    if (CRAM_MAJOR_VERS(fd->version) >= 4) {
        int64_t i64;
        if ((s = fd->vv.varint_decode64_crc(fd, &i64, &crc))== -1) return NULL; else rd+=s;
        c2.ref_seq_start = i64;
        if ((s = fd->vv.varint_decode64_crc(fd, &i64, &crc)) == -1) return NULL; else rd+=s;
        c2.ref_seq_span = i64;
    } else {
        int32_t i32;
        if ((s = fd->vv.varint_decode32_crc(fd, &i32, &crc))== -1) return NULL; else rd+=s;
        c2.ref_seq_start = i32;
        if ((s = fd->vv.varint_decode32_crc(fd, &i32, &crc)) == -1) return NULL; else rd+=s;
        c2.ref_seq_span = i32;
    }
    if ((s = fd->vv.varint_decode32_crc(fd, &c2.num_records, &crc))  == -1) return NULL; else rd+=s;

    if (CRAM_MAJOR_VERS(fd->version) == 1) {
        c2.record_counter = 0;
        c2.num_bases = 0;
    } else {
        if (CRAM_MAJOR_VERS(fd->version) >= 3) {
            if ((s = fd->vv.varint_decode64_crc(fd, &c2.record_counter, &crc)) == -1)
                return NULL;
            else
                rd += s;
        } else {
            int32_t i32;
            if ((s = fd->vv.varint_decode32_crc(fd, &i32, &crc)) == -1)
                return NULL;
            else
                rd += s;
            c2.record_counter = i32;
        }

        if ((s = fd->vv.varint_decode64_crc(fd, &c2.num_bases, &crc))== -1)
            return NULL;
        else
            rd += s;
    }
    if ((s = fd->vv.varint_decode32_crc(fd, &c2.num_blocks, &crc))   == -1)
        return NULL;
    else
        rd+=s;
    if ((s = fd->vv.varint_decode32_crc(fd, &c2.num_landmarks, &crc))== -1)
        return NULL;
    else
        rd+=s;

    if (c2.num_landmarks < 0 || c2.num_landmarks >= SIZE_MAX / sizeof(int32_t))
        return NULL;

    if (!(c = calloc(1, sizeof(*c))))
        return NULL;

    *c = c2;
#ifdef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
    if (c->num_landmarks > FUZZ_ALLOC_LIMIT/sizeof(int32_t)) {
        fd->err = errno = ENOMEM;
        cram_free_container(c);
        return NULL;
    }
#endif
    if (c->num_landmarks && !(c->landmark = malloc(c->num_landmarks * sizeof(int32_t)))) {
        fd->err = errno;
        cram_free_container(c);
        return NULL;
    }
    for (i = 0; i < c->num_landmarks; i++) {
        if ((s = fd->vv.varint_decode32_crc(fd, &c->landmark[i], &crc)) == -1) {
            cram_free_container(c);
            return NULL;
        } else {
            rd += s;
        }
    }

    if (CRAM_MAJOR_VERS(fd->version) >= 3) {
        if (-1 == int32_decode(fd, (int32_t *)&c->crc32)) {
            cram_free_container(c);
            return NULL;
        } else {
            rd+=4;
        }

#ifdef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
        // Pretend the CRC was OK so the fuzzer doesn't have to get it right
        crc = c->crc32;
#endif

        if (crc != c->crc32) {
            hts_log_error("Container header CRC32 failure");
            cram_free_container(c);
            return NULL;
        }
    }

    c->offset = rd;
    c->slices = NULL;
    c->slice = NULL;
    c->curr_slice = 0;
    c->max_slice = c->num_landmarks;
    c->slice_rec = 0;
    c->curr_rec = 0;
    c->max_rec = 0;

    if (c->ref_seq_id == -2) {
        c->multi_seq = 1;
        fd->multi_seq = 1;
    }

    fd->empty_container =
        (c->num_records == 0 &&
         c->ref_seq_id == -1 &&
         c->ref_seq_start == 0x454f46 /* EOF */) ? 1 : 0;

    return c;
}


/* MAXIMUM storage size needed for the container. */
int cram_container_size(cram_container *c) {
    return 55 + 5*c->num_landmarks;
}


/*
 * Stores the container structure in dat and returns *size as the
 * number of bytes written to dat[].  The input size of dat is also
 * held in *size and should be initialised to cram_container_size(c).
 *
 * Returns 0 on success;
 *        -1 on failure
 */
int cram_store_container(cram_fd *fd, cram_container *c, char *dat, int *size)
{
    char *cp = (char *)dat;
    int i;

    // Check the input buffer is large enough according to our stated
    // requirements. (NOTE: it may actually take less.)
    if (cram_container_size(c) > *size)
        return -1;

    if (CRAM_MAJOR_VERS(fd->version) == 1) {
        cp += itf8_put(cp, c->length);
    } else {
        *(int32_t *)cp = le_int4(c->length);
        cp += 4;
    }
    if (c->multi_seq) {
        cp += fd->vv.varint_put32(cp, NULL, -2);
        cp += fd->vv.varint_put32(cp, NULL, 0);
        cp += fd->vv.varint_put32(cp, NULL, 0);
    } else {
        cp += fd->vv.varint_put32s(cp, NULL, c->ref_seq_id);
        if (CRAM_MAJOR_VERS(fd->version) >= 4) {
            cp += fd->vv.varint_put64(cp, NULL, c->ref_seq_start);
            cp += fd->vv.varint_put64(cp, NULL, c->ref_seq_span);
        } else {
            cp += fd->vv.varint_put32(cp, NULL, c->ref_seq_start);
            cp += fd->vv.varint_put32(cp, NULL, c->ref_seq_span);
        }
    }
    cp += fd->vv.varint_put32(cp, NULL, c->num_records);
    if (CRAM_MAJOR_VERS(fd->version) == 2) {
        cp += fd->vv.varint_put64(cp, NULL, c->record_counter);
    } else if (CRAM_MAJOR_VERS(fd->version) >= 3) {
        cp += fd->vv.varint_put32(cp, NULL, c->record_counter);
    }
    cp += fd->vv.varint_put64(cp, NULL, c->num_bases);
    cp += fd->vv.varint_put32(cp, NULL, c->num_blocks);
    cp += fd->vv.varint_put32(cp, NULL, c->num_landmarks);
    for (i = 0; i < c->num_landmarks; i++)
        cp += fd->vv.varint_put32(cp, NULL, c->landmark[i]);

    if (CRAM_MAJOR_VERS(fd->version) >= 3) {
        c->crc32 = crc32(0L, (uc *)dat, cp-dat);
        cp[0] =  c->crc32        & 0xff;
        cp[1] = (c->crc32 >>  8) & 0xff;
        cp[2] = (c->crc32 >> 16) & 0xff;
        cp[3] = (c->crc32 >> 24) & 0xff;
        cp += 4;
    }

    *size = cp-dat; // actual used size

    return 0;
}


/*
 * Writes a container structure.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_write_container(cram_fd *fd, cram_container *c) {
    char buf_a[1024], *buf = buf_a, *cp;
    int i;

    if (61 + c->num_landmarks * 10 >= 1024) {
        buf = malloc(61 + c->num_landmarks * 10);
        if (!buf)
            return -1;
    }
    cp = buf;

    if (CRAM_MAJOR_VERS(fd->version) == 1) {
        cp += itf8_put(cp, c->length);
    } else if (CRAM_MAJOR_VERS(fd->version) <= 3) {
        *(int32_t *)cp = le_int4(c->length);
        cp += 4;
    } else {
        cp += fd->vv.varint_put32(cp, NULL, c->length);
    }
    if (c->multi_seq) {
        cp += fd->vv.varint_put32(cp, NULL, (uint32_t)-2);
        cp += fd->vv.varint_put32(cp, NULL, 0);
        cp += fd->vv.varint_put32(cp, NULL, 0);
    } else {
        cp += fd->vv.varint_put32s(cp, NULL, c->ref_seq_id);
        if (CRAM_MAJOR_VERS(fd->version) >= 4) {
            cp += fd->vv.varint_put64(cp, NULL, c->ref_seq_start);
            cp += fd->vv.varint_put64(cp, NULL, c->ref_seq_span);
        } else {
            cp += fd->vv.varint_put32(cp, NULL, c->ref_seq_start);
            cp += fd->vv.varint_put32(cp, NULL, c->ref_seq_span);
        }
    }
    cp += fd->vv.varint_put32(cp, NULL, c->num_records);
    if (CRAM_MAJOR_VERS(fd->version) >= 3)
        cp += fd->vv.varint_put64(cp, NULL, c->record_counter);
    else
        cp += fd->vv.varint_put32(cp, NULL, c->record_counter);
    cp += fd->vv.varint_put64(cp, NULL, c->num_bases);
    cp += fd->vv.varint_put32(cp, NULL, c->num_blocks);
    cp += fd->vv.varint_put32(cp, NULL, c->num_landmarks);
    for (i = 0; i < c->num_landmarks; i++)
        cp += fd->vv.varint_put32(cp, NULL, c->landmark[i]);

    if (CRAM_MAJOR_VERS(fd->version) >= 3) {
        c->crc32 = crc32(0L, (uc *)buf, cp-buf);
        cp[0] =  c->crc32        & 0xff;
        cp[1] = (c->crc32 >>  8) & 0xff;
        cp[2] = (c->crc32 >> 16) & 0xff;
        cp[3] = (c->crc32 >> 24) & 0xff;
        cp += 4;
    }

    if (cp-buf != hwrite(fd->fp, buf, cp-buf)) {
        if (buf != buf_a)
            free(buf);
        return -1;
    }

    if (buf != buf_a)
        free(buf);

    return 0;
}

typedef struct {
    cram_fd *fd;
    cram_container *c;
} cram_job;

void *cram_flush_thread(void *arg) {
    // cram_job *j = (cram_job *)arg;

    // /* Encode the container blocks and generate compression header */
    // if (0 != cram_encode_container(j->fd, j->c)) {
    //     hts_log_error("Call to cram_encode_container failed");
    //     return NULL;
    // }

    return 0;
}

// Note: called while metrics_lock is held.
// Will be left in this state too, but may temporarily unlock.
void reset_metrics(cram_fd *fd) {
    int i;

    if (fd->pool) {
        // If multi-threaded we have multiple blocks being
        // compressed already and several on the to-do list
        // (fd->rqueue->pending).  It's tricky to reset the
        // metrics exactly the correct point, so instead we
        // just flush the pool, reset, and then continue again.

        // Don't bother starting a new trial before then though.
        for (i = 0; i < DS_END; i++) {
            cram_metrics *m = fd->m[i];
            if (!m)
                continue;
            m->next_trial = 999;
        }

        pthread_mutex_unlock(&fd->metrics_lock);
        hts_tpool_process_flush(fd->rqueue);
        pthread_mutex_lock(&fd->metrics_lock);
    }

    for (i = 0; i < DS_END; i++) {
        cram_metrics *m = fd->m[i];
        if (!m)
            continue;

        m->trial = NTRIALS;
        m->next_trial = TRIAL_SPAN;
        m->revised_method = 0;
        m->unpackable = 0;

        memset(m->sz, 0, sizeof(m->sz));
    }
}

void cram_free_compression_header(cram_block_compression_hdr *hdr) {
    int i;

    if (hdr->landmark)
        free(hdr->landmark);

    if (hdr->preservation_map)
        kh_destroy(map, hdr->preservation_map);

    for (i = 0; i < CRAM_MAP_HASH; i++) {
        cram_map *m, *m2;
        for (m = hdr->rec_encoding_map[i]; m; m = m2) {
            m2 = m->next;
            if (m->codec)
                m->codec->free(m->codec);
            free(m);
        }
    }

    for (i = 0; i < CRAM_MAP_HASH; i++) {
        cram_map *m, *m2;
        for (m = hdr->tag_encoding_map[i]; m; m = m2) {
            m2 = m->next;
            if (m->codec)
                m->codec->free(m->codec);
            free(m);
        }
    }

    for (i = 0; i < DS_END; i++) {
        if (hdr->codecs[i])
            hdr->codecs[i]->free(hdr->codecs[i]);
    }

    if (hdr->TL)
        free(hdr->TL);
    if (hdr->TD_blk)
        cram_free_block(hdr->TD_blk);
    if (hdr->TD_hash)
        kh_destroy(m_s2i, hdr->TD_hash);
    if (hdr->TD_keys)
        string_pool_destroy(hdr->TD_keys);

    free(hdr);
}


/* ----------------------------------------------------------------------
 * Slices and slice headers
 */

void cram_free_slice_header(cram_block_slice_hdr *hdr) {
    if (!hdr)
        return;

    if (hdr->block_content_ids)
        free(hdr->block_content_ids);

    free(hdr);

    return;
}

void cram_free_slice(cram_slice *s) {
    if (!s)
        return;

    if (s->hdr_block)
        cram_free_block(s->hdr_block);

    if (s->block) {
        int i;

        if (s->hdr) {
            for (i = 0; i < s->hdr->num_blocks; i++) {
                if (i > 0 && s->block[i] == s->block[0])
                    continue;
                cram_free_block(s->block[i]);
            }
        }
        free(s->block);
    }

    if (s->block_by_id)
        free(s->block_by_id);

    if (s->hdr)
        cram_free_slice_header(s->hdr);

    if (s->seqs_blk)
        cram_free_block(s->seqs_blk);

    if (s->qual_blk)
        cram_free_block(s->qual_blk);

    if (s->name_blk)
        cram_free_block(s->name_blk);

    if (s->aux_blk)
        cram_free_block(s->aux_blk);

    if (s->base_blk)
        cram_free_block(s->base_blk);

    if (s->soft_blk)
        cram_free_block(s->soft_blk);

    if (s->cigar)
        free(s->cigar);

    if (s->crecs)
        free(s->crecs);

    if (s->features)
        free(s->features);

    if (s->TN)
        free(s->TN);

    if (s->pair_keys)
        string_pool_destroy(s->pair_keys);

    if (s->pair[0])
        kh_destroy(m_s2i, s->pair[0]);
    if (s->pair[1])
        kh_destroy(m_s2i, s->pair[1]);

    if (s->aux_block)
        free(s->aux_block);

    free(s);
}

/*
 * Creates a new empty slice in memory, for subsequent writing to
 * disk.
 *
 * Returns cram_slice ptr on success
 *         NULL on failure
 */
cram_slice *cram_new_slice(enum cram_content_type type, int nrecs) {
    cram_slice *s = calloc(1, sizeof(*s));
    if (!s)
        return NULL;

    if (!(s->hdr = (cram_block_slice_hdr *)calloc(1, sizeof(*s->hdr))))
        goto err;
    s->hdr->content_type = type;

    s->hdr_block = NULL;
    s->block = NULL;
    s->block_by_id = NULL;
    s->last_apos = 0;
    if (!(s->crecs = malloc(nrecs * sizeof(cram_record))))  goto err;
    s->cigar_alloc = 1024;
    if (!(s->cigar = malloc(s->cigar_alloc * sizeof(*s->cigar)))) goto err;
    s->ncigar = 0;

    if (!(s->seqs_blk = cram_new_block(EXTERNAL, 0)))       goto err;
    if (!(s->qual_blk = cram_new_block(EXTERNAL, DS_QS)))   goto err;
    if (!(s->name_blk = cram_new_block(EXTERNAL, DS_RN)))   goto err;
    if (!(s->aux_blk  = cram_new_block(EXTERNAL, DS_aux)))  goto err;
    if (!(s->base_blk = cram_new_block(EXTERNAL, DS_IN)))   goto err;
    if (!(s->soft_blk = cram_new_block(EXTERNAL, DS_SC)))   goto err;

    s->features = NULL;
    s->nfeatures = s->afeatures = 0;

#ifndef TN_external
    s->TN = NULL;
    s->nTN = s->aTN = 0;
#endif

    // Volatile keys as we do realloc in dstring
    if (!(s->pair_keys = string_pool_create(8192))) goto err;
    if (!(s->pair[0] = kh_init(m_s2i)))             goto err;
    if (!(s->pair[1] = kh_init(m_s2i)))             goto err;

#ifdef BA_external
    s->BA_len = 0;
#endif

    return s;

 err:
    if (s)
        cram_free_slice(s);

    return NULL;
}


/* ----------------------------------------------------------------------
 * CRAM file definition (header)
 */

/*
 * Reads a CRAM file definition structure.
 * Returns file_def ptr on success
 *         NULL on failure
 */
cram_file_def *cram_read_file_def(cram_fd *fd) {
    cram_file_def *def = malloc(sizeof(*def));
    if (!def)
        return NULL;

    if (26 != hread(fd->fp, &def->magic[0], 26)) {
        free(def);
        return NULL;
    }

    if (memcmp(def->magic, "CRAM", 4) != 0) {
        free(def);
        return NULL;
    }

    if (def->major_version > 4) {
        hts_log_error("CRAM version number mismatch. Expected 1.x, 2.x, 3.x or 4.x, got %d.%d",
                      def->major_version, def->minor_version);
        free(def);
        return NULL;
    }

    fd->first_container += 26;
    fd->curr_position = fd->first_container;
    fd->last_slice = 0;

    return def;
}

/*
 * Writes a cram_file_def structure to cram_fd.
 * Returns 0 on success
 *        -1 on failure
 */
int cram_write_file_def(cram_fd *fd, cram_file_def *def) {
    return (hwrite(fd->fp, &def->magic[0], 26) == 26) ? 0 : -1;
}

void cram_free_file_def(cram_file_def *def) {
    if (def) free(def);
}

/* ----------------------------------------------------------------------
 * The top-level cram opening, closing and option handling
 */

/*
 * Sets CRAM variable sized integer decode function tables.
 * CRAM 1, 2, and 3.x all used ITF8 for uint32 and UTF8 for uint64.
 * CRAM 4.x uses the same encoding mechanism for 32-bit and 64-bit
 * (or anything inbetween), but also now supports signed values.
 *
 * Version is the CRAM major version number.
 * vv is the vector table (probably &cram_fd->vv)
 */
static void cram_init_varint(varint_vec *vv, int version) {
    if (version >= 4) {
        vv->varint_get32 = uint7_get_32; // FIXME: varint.h API should be size agnostic
        vv->varint_get32s = sint7_get_32;
        vv->varint_get64 = uint7_get_64;
        vv->varint_get64s = sint7_get_64;
        vv->varint_put32 = uint7_put_32;
        vv->varint_put32s = sint7_put_32;
        vv->varint_put64 = uint7_put_64;
        vv->varint_put64s = sint7_put_64;
        vv->varint_put32_blk = uint7_put_blk_32;
        vv->varint_put32s_blk = sint7_put_blk_32;
        vv->varint_put64_blk = uint7_put_blk_64;
        vv->varint_put64s_blk = sint7_put_blk_64;
        vv->varint_size = uint7_size;
        vv->varint_decode32_crc = uint7_decode_crc32;
        vv->varint_decode32s_crc = sint7_decode_crc32;
        vv->varint_decode64_crc = uint7_decode_crc64;
    } else {
        vv->varint_get32 = safe_itf8_get;
        vv->varint_get32s = safe_itf8_get;
        vv->varint_get64 = safe_ltf8_get;
        vv->varint_get64s = safe_ltf8_get;
        vv->varint_put32 = safe_itf8_put;
        vv->varint_put32s = safe_itf8_put;
        vv->varint_put64 = safe_ltf8_put;
        vv->varint_put64s = safe_ltf8_put;
        vv->varint_put32_blk = itf8_put_blk;
        vv->varint_put32s_blk = itf8_put_blk;
        vv->varint_put64_blk = ltf8_put_blk;
        vv->varint_put64s_blk = ltf8_put_blk;
        vv->varint_size = itf8_size;
        vv->varint_decode32_crc = itf8_decode_crc;
        vv->varint_decode32s_crc = itf8_decode_crc;
        vv->varint_decode64_crc = ltf8_decode_crc;
    }
}

/*
 * Initialises the lookup tables. These could be global statics, but they're
 * clumsy to setup in a multi-threaded environment unless we generate
 * verbatim code and include that.
 */
static void cram_init_tables(cram_fd *fd) {
    int i;

    memset(fd->L1, 4, 256);
    fd->L1['A'] = 0; fd->L1['a'] = 0;
    fd->L1['C'] = 1; fd->L1['c'] = 1;
    fd->L1['G'] = 2; fd->L1['g'] = 2;
    fd->L1['T'] = 3; fd->L1['t'] = 3;

    memset(fd->L2, 5, 256);
    fd->L2['A'] = 0; fd->L2['a'] = 0;
    fd->L2['C'] = 1; fd->L2['c'] = 1;
    fd->L2['G'] = 2; fd->L2['g'] = 2;
    fd->L2['T'] = 3; fd->L2['t'] = 3;
    fd->L2['N'] = 4; fd->L2['n'] = 4;

    if (CRAM_MAJOR_VERS(fd->version) == 1) {
        for (i = 0; i < 0x200; i++) {
            int f = 0;

            if (i & CRAM_FPAIRED)      f |= BAM_FPAIRED;
            if (i & CRAM_FPROPER_PAIR) f |= BAM_FPROPER_PAIR;
            if (i & CRAM_FUNMAP)       f |= BAM_FUNMAP;
            if (i & CRAM_FREVERSE)     f |= BAM_FREVERSE;
            if (i & CRAM_FREAD1)       f |= BAM_FREAD1;
            if (i & CRAM_FREAD2)       f |= BAM_FREAD2;
            if (i & CRAM_FSECONDARY)   f |= BAM_FSECONDARY;
            if (i & CRAM_FQCFAIL)      f |= BAM_FQCFAIL;
            if (i & CRAM_FDUP)         f |= BAM_FDUP;

            fd->bam_flag_swap[i]  = f;
        }

        for (i = 0; i < 0x1000; i++) {
            int g = 0;

            if (i & BAM_FPAIRED)           g |= CRAM_FPAIRED;
            if (i & BAM_FPROPER_PAIR)  g |= CRAM_FPROPER_PAIR;
            if (i & BAM_FUNMAP)        g |= CRAM_FUNMAP;
            if (i & BAM_FREVERSE)      g |= CRAM_FREVERSE;
            if (i & BAM_FREAD1)        g |= CRAM_FREAD1;
            if (i & BAM_FREAD2)        g |= CRAM_FREAD2;
            if (i & BAM_FSECONDARY)    g |= CRAM_FSECONDARY;
            if (i & BAM_FQCFAIL)       g |= CRAM_FQCFAIL;
            if (i & BAM_FDUP)          g |= CRAM_FDUP;

            fd->cram_flag_swap[i] = g;
        }
    } else {
        /* NOP */
        for (i = 0; i < 0x1000; i++)
            fd->bam_flag_swap[i] = i;
        for (i = 0; i < 0x1000; i++)
            fd->cram_flag_swap[i] = i;
    }

    memset(fd->cram_sub_matrix, 4, 32*32);
    for (i = 0; i < 32; i++) {
        fd->cram_sub_matrix[i]['A'&0x1f]=0;
        fd->cram_sub_matrix[i]['C'&0x1f]=1;
        fd->cram_sub_matrix[i]['G'&0x1f]=2;
        fd->cram_sub_matrix[i]['T'&0x1f]=3;
        fd->cram_sub_matrix[i]['N'&0x1f]=4;
    }
    for (i = 0; i < 20; i+=4) {
        int j;
        for (j = 0; j < 20; j++) {
            fd->cram_sub_matrix["ACGTN"[i>>2]&0x1f][j]=3;
            fd->cram_sub_matrix["ACGTN"[i>>2]&0x1f][j]=3;
            fd->cram_sub_matrix["ACGTN"[i>>2]&0x1f][j]=3;
            fd->cram_sub_matrix["ACGTN"[i>>2]&0x1f][j]=3;
        }
        fd->cram_sub_matrix["ACGTN"[i>>2]&0x1f][CRAM_SUBST_MATRIX[i+0]&0x1f]=0;
        fd->cram_sub_matrix["ACGTN"[i>>2]&0x1f][CRAM_SUBST_MATRIX[i+1]&0x1f]=1;
        fd->cram_sub_matrix["ACGTN"[i>>2]&0x1f][CRAM_SUBST_MATRIX[i+2]&0x1f]=2;
        fd->cram_sub_matrix["ACGTN"[i>>2]&0x1f][CRAM_SUBST_MATRIX[i+3]&0x1f]=3;
    }

    cram_init_varint(&fd->vv, CRAM_MAJOR_VERS(fd->version));
}

/*
 * Seek within a CRAM file.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_seek(cram_fd *fd, off_t offset, int whence) {
    char buf[65536];

    fd->ooc = 0;

    cram_drain_rqueue(fd);

    if (hseek(fd->fp, offset, whence) >= 0) {
        return 0;
    }

    if (!(whence == SEEK_CUR && offset >= 0))
        return -1;

    /* Couldn't fseek, but we're in SEEK_CUR mode so read instead */
    while (offset > 0) {
        int len = MIN(65536, offset);
        if (len != hread(fd->fp, buf, len))
            return -1;
        offset -= len;
    }

    return 0;
}

/*
 * Returns 1 if we hit an EOF while reading.
 */
int cram_eof(cram_fd *fd) {
    return fd->eof;
}


/*
 * Sets options on the cram_fd. See CRAM_OPT_* definitions in cram_structs.h.
 * Use this immediately after opening.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_set_option(cram_fd *fd, enum hts_fmt_option opt, ...) {
    int r;
    va_list args;

    va_start(args, opt);
    r = cram_set_voption(fd, opt, args);
    va_end(args);

    return r;
}

/*
 * Sets options on the cram_fd. See CRAM_OPT_* definitions in cram_structs.h.
 * Use this immediately after opening.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_set_voption(cram_fd *fd, enum hts_fmt_option opt, va_list args) {
    refs_t *refs;

    if (!fd) {
        errno = EBADF;
        return -1;
    }

    switch (opt) {
    case CRAM_OPT_DECODE_MD:
        fd->decode_md = va_arg(args, int);
        break;

    case CRAM_OPT_PREFIX:
        if (fd->prefix)
            free(fd->prefix);
        if (!(fd->prefix = strdup(va_arg(args, char *))))
            return -1;
        break;

    case CRAM_OPT_VERBOSITY:
        break;

    case CRAM_OPT_SEQS_PER_SLICE:
        fd->seqs_per_slice = va_arg(args, int);
        if (fd->bases_per_slice == BASES_PER_SLICE)
            fd->bases_per_slice = fd->seqs_per_slice * 500;
        break;

    case CRAM_OPT_BASES_PER_SLICE:
        fd->bases_per_slice = va_arg(args, int);
        break;

    case CRAM_OPT_SLICES_PER_CONTAINER:
        fd->slices_per_container = va_arg(args, int);
        break;

    case CRAM_OPT_EMBED_REF:
        fd->embed_ref = va_arg(args, int);
        break;

    case CRAM_OPT_NO_REF:
        fd->no_ref = va_arg(args, int);
        break;

    case CRAM_OPT_POS_DELTA:
        fd->ap_delta = va_arg(args, int);
        break;

    case CRAM_OPT_IGNORE_MD5:
        fd->ignore_md5 = va_arg(args, int);
        break;

    case CRAM_OPT_LOSSY_NAMES:
        fd->lossy_read_names = va_arg(args, int);
        // Currently lossy read names required paired (attached) reads.
        // TLEN 0 or being 1 out causes read pairs to be detached, breaking
        // the lossy read name compression, so we have extra options to
        // slacken the exact TLEN round-trip checks.
        fd->tlen_approx = fd->lossy_read_names;
        fd->tlen_zero = fd->lossy_read_names;
        break;

    case CRAM_OPT_USE_BZIP2:
        fd->use_bz2 = va_arg(args, int);
        break;

    case CRAM_OPT_USE_RANS:
        fd->use_rans = va_arg(args, int);
        break;

    case CRAM_OPT_USE_TOK:
        fd->use_tok = va_arg(args, int);
        break;

    case CRAM_OPT_USE_FQZ:
        fd->use_fqz = va_arg(args, int);
        break;

    case CRAM_OPT_USE_ARITH:
        fd->use_arith = va_arg(args, int);
        break;

    case CRAM_OPT_USE_LZMA:
        fd->use_lzma = va_arg(args, int);
        break;

    case CRAM_OPT_SHARED_REF:
        fd->shared_ref = 1;
        refs = va_arg(args, refs_t *);
        if (refs != fd->refs) {
            if (fd->refs)
                refs_free(fd->refs);
            fd->refs = refs;
            fd->refs->count++;
        }
        break;

    case CRAM_OPT_RANGE: {
        int r = cram_seek_to_refpos(fd, va_arg(args, cram_range *));
        pthread_mutex_lock(&fd->range_lock);
        if (fd->range.refid != -2)
            fd->required_fields |= SAM_POS;
        pthread_mutex_unlock(&fd->range_lock);
        return r;
    }

    case CRAM_OPT_RANGE_NOSEEK: {
        // As per CRAM_OPT_RANGE, but no seeking
        pthread_mutex_lock(&fd->range_lock);
        cram_range *r = va_arg(args, cram_range *);
        fd->range = *r;
        if (r->refid == HTS_IDX_NOCOOR) {
            fd->range.refid = -1;
            fd->range.start = 0;
        } else if (r->refid == HTS_IDX_START || r->refid == HTS_IDX_REST) {
            fd->range.refid = -2; // special case in cram_next_slice
        }
        if (fd->range.refid != -2)
            fd->required_fields |= SAM_POS;
        fd->ooc = 0;
        fd->eof = 0;
        pthread_mutex_unlock(&fd->range_lock);
        return 0;
    }

    case CRAM_OPT_REFERENCE:
        return cram_load_reference(fd, va_arg(args, char *));

    case CRAM_OPT_VERSION: {
        int major, minor;
        char *s = va_arg(args, char *);
        if (2 != sscanf(s, "%d.%d", &major, &minor)) {
            hts_log_error("Malformed version string %s", s);
            return -1;
        }
        if (!((major == 1 &&  minor == 0) ||
              (major == 2 && (minor == 0 || minor == 1)) ||
              (major == 3 && (minor == 0 || minor == 1)) ||
              (major == 4 &&  minor == 0))) {
            hts_log_error("Unknown version string; use 1.0, 2.0, 2.1, 3.0, 3.1 or 4.0");
            errno = EINVAL;
            return -1;
        }

        if (major > 3 || (major == 3 && minor > 1)) {
            hts_log_warning(
                "CRAM version %s is still a draft and subject to change.\n"
                "This is a technology demonstration that should not be "
                "used for archival data.", s);
        }

        fd->version = major*256 + minor;

        fd->use_rans = (CRAM_MAJOR_VERS(fd->version) >= 3) ? 1 : 0;

        fd->use_tok = ((CRAM_MAJOR_VERS(fd->version) == 3 &&
                        CRAM_MINOR_VERS(fd->version) >= 1) ||
                        CRAM_MAJOR_VERS(fd->version) >= 4) ? 1 : 0;
        cram_init_tables(fd);

        break;
    }

    case CRAM_OPT_MULTI_SEQ_PER_SLICE:
        fd->multi_seq_user = fd->multi_seq = va_arg(args, int);
        break;

    case CRAM_OPT_NTHREADS: {
        int nthreads =  va_arg(args, int);
        if (nthreads >= 1) {
            if (!(fd->pool = hts_tpool_init(nthreads)))
                return -1;

            fd->rqueue = hts_tpool_process_init(fd->pool, nthreads*2, 0);
            pthread_mutex_init(&fd->metrics_lock, NULL);
            pthread_mutex_init(&fd->ref_lock, NULL);
            pthread_mutex_init(&fd->range_lock, NULL);
            pthread_mutex_init(&fd->bam_list_lock, NULL);
            fd->shared_ref = 1;
            fd->own_pool = 1;
        }
        break;
    }

    case CRAM_OPT_THREAD_POOL: {
        htsThreadPool *p = va_arg(args, htsThreadPool *);
        fd->pool = p ? p->pool : NULL;
        if (fd->pool) {
            fd->rqueue = hts_tpool_process_init(fd->pool,
                                                p->qsize ? p->qsize : hts_tpool_size(fd->pool)*2,
                                                0);
            pthread_mutex_init(&fd->metrics_lock, NULL);
            pthread_mutex_init(&fd->ref_lock, NULL);
            pthread_mutex_init(&fd->range_lock, NULL);
            pthread_mutex_init(&fd->bam_list_lock, NULL);
        }
        fd->shared_ref = 1; // Needed to avoid clobbering ref between threads
        fd->own_pool = 0;

        //fd->qsize = 1;
        //fd->decoded = calloc(fd->qsize, sizeof(cram_container *));
        //hts_tpool_dispatch(fd->pool, cram_decoder_thread, fd);
        break;
    }

    case CRAM_OPT_REQUIRED_FIELDS:
        fd->required_fields = va_arg(args, int);
        if (fd->range.refid != -2)
            fd->required_fields |= SAM_POS;
        break;

    case CRAM_OPT_STORE_MD:
        fd->store_md = va_arg(args, int);
        break;

    case CRAM_OPT_STORE_NM:
        fd->store_nm = va_arg(args, int);
        break;

    case HTS_OPT_COMPRESSION_LEVEL:
        fd->level = va_arg(args, int);
        break;

    case HTS_OPT_PROFILE: {
        enum hts_profile_option prof = va_arg(args, int);
        switch (prof) {
        case HTS_PROFILE_FAST:
            if (fd->level == CRAM_DEFAULT_LEVEL) fd->level = 1;
            fd->use_tok = 0;
            fd->seqs_per_slice = 10000;
            break;

        case HTS_PROFILE_NORMAL:
            break;

        case HTS_PROFILE_SMALL:
            if (fd->level == CRAM_DEFAULT_LEVEL) fd->level = 6;
            fd->use_bz2 = 1;
            fd->use_fqz = 1;
            fd->seqs_per_slice = 25000;
            break;

        case HTS_PROFILE_ARCHIVE:
            if (fd->level == CRAM_DEFAULT_LEVEL) fd->level = 7;
            fd->use_bz2 = 1;
            fd->use_fqz = 1;
            fd->use_arith = 1;
            if (fd->level > 7)
                fd->use_lzma = 1;
            fd->seqs_per_slice = 100000;
            break;
        }

        if (fd->bases_per_slice == BASES_PER_SLICE)
            fd->bases_per_slice = fd->seqs_per_slice * 500;
        break;
    }

    default:
        hts_log_error("Unknown CRAM option code %d", opt);
        errno = EINVAL;
        return -1;
    }

    return 0;
}

int cram_check_EOF(cram_fd *fd)
{
    // Byte 9 in these templates is & with 0x0f to resolve differences
    // between ITF-8 interpretations between early Java and C
    // implementations of CRAM
    static const unsigned char TEMPLATE_2_1[30] = {
        0x0b, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0x0f, 0xe0,
        0x45, 0x4f, 0x46, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00,
        0x01, 0x00, 0x06, 0x06, 0x01, 0x00, 0x01, 0x00, 0x01, 0x00
    };
    static const unsigned char TEMPLATE_3[38] = {
        0x0f, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0x0f, 0xe0,
        0x45, 0x4f, 0x46, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x05,
        0xbd, 0xd9, 0x4f, 0x00, 0x01, 0x00, 0x06, 0x06, 0x01, 0x00,
        0x01, 0x00, 0x01, 0x00, 0xee, 0x63, 0x01, 0x4b
    };

    unsigned char buf[38]; // max(sizeof TEMPLATE_*)

    uint8_t major = CRAM_MAJOR_VERS(fd->version);
    uint8_t minor = CRAM_MINOR_VERS(fd->version);

    const unsigned char *template;
    ssize_t template_len;
    if ((major < 2) ||
        (major == 2 && minor == 0)) {
        return 3; // No EOF support in cram versions less than 2.1
    } else if (major == 2 && minor == 1) {
        template = TEMPLATE_2_1;
        template_len = sizeof TEMPLATE_2_1;
    } else {
        template = TEMPLATE_3;
        template_len = sizeof TEMPLATE_3;
    }

    off_t offset = htell(fd->fp);
    if (hseek(fd->fp, -template_len, SEEK_END) < 0) {
        if (errno == ESPIPE) {
            hclearerr(fd->fp);
            return 2;
        }
        else {
            return -1;
        }
    }
    if (hread(fd->fp, buf, template_len) != template_len) return -1;
    if (hseek(fd->fp, offset, SEEK_SET) < 0) return -1;
    buf[8] &= 0x0f;
    return (memcmp(template, buf, template_len) == 0)? 1 : 0;
}
