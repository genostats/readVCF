/*
Copyright (c) 2012-2020, 2022-2023 Genome Research Ltd.
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
 * - In-memory decoding of CRAM data structures.
 * - Iterator for reading CRAM record by record.
 */

#define HTS_BUILDING_LIBRARY // Enables HTSLIB_EXPORT, see htslib/hts_defs.h
#include <config.h>

#include <stdio.h>
#include <errno.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>
#include <stdint.h>
#include <inttypes.h>

#include "cram/cram.h"
#include "cram/os.h"
#include "htslib/hts.h"

//Whether CIGAR has just M or uses = and X to indicate match and mismatch
//#define USE_X

/* ----------------------------------------------------------------------
 * CRAM compression headers
 */

/*
 * Decodes the Tag Dictionary record in the preservation map
 * Updates the cram compression header.
 *
 * Returns number of bytes decoded on success
 *        -1 on failure
 */
int cram_decode_TD(cram_fd *fd, char *cp, const char *endp,
                   cram_block_compression_hdr *h) {
    char *op = cp;
    unsigned char *dat;
    cram_block *b;
    int32_t blk_size = 0;
    int nTL, i, sz, err = 0;

    if (!(b = cram_new_block(0, 0)))
        return -1;

    if (h->TD_blk || h->TL) {
        hts_log_warning("More than one TD block found in compression header");
        cram_free_block(h->TD_blk);
        free(h->TL);
        h->TD_blk = NULL;
        h->TL = NULL;
    }

    /* Decode */
    blk_size = fd->vv.varint_get32(&cp, endp, &err);
    if (!blk_size) {
        h->nTL = 0;
        cram_free_block(b);
        return cp - op;
    }

    if (err || blk_size < 0 || endp - cp < blk_size) {
        cram_free_block(b);
        return -1;
    }

    BLOCK_APPEND(b, cp, blk_size);
    cp += blk_size;
    sz = cp - op;
    // Force nul termination if missing
    if (BLOCK_DATA(b)[BLOCK_SIZE(b)-1])
        BLOCK_APPEND_CHAR(b, '\0');

    /* Set up TL lookup table */
    dat = BLOCK_DATA(b);

    // Count
    for (nTL = i = 0; i < BLOCK_SIZE(b); i++) {
        nTL++;
        while (dat[i])
            i++;
    }

    // Copy
    if (!(h->TL = calloc(nTL, sizeof(*h->TL)))) {
        cram_free_block(b);
        return -1;
    }
    for (nTL = i = 0; i < BLOCK_SIZE(b); i++) {
        h->TL[nTL++] = &dat[i];
        while (dat[i])
            i++;
    }
    h->TD_blk = b;
    h->nTL = nTL;

    return sz;

 block_err:
    cram_free_block(b);
    return -1;
}

/* ----------------------------------------------------------------------
 * Primary CRAM sequence decoder
 */

static inline int add_md_char(cram_slice *s, int decode_md, char c, int32_t *md_dist) {
    if (decode_md) {
        BLOCK_APPEND_UINT(s->aux_blk, *md_dist);
        BLOCK_APPEND_CHAR(s->aux_blk, c);
        *md_dist = 0;
    }
    return 0;

 block_err:
    return -1;
}

typedef struct {
    cram_fd *fd;
    cram_container *c;
    cram_slice *s;
    sam_hdr_t *h;
    int exit_code;
} cram_decode_job;

/*
 * Drains and frees the decode read-queue for a multi-threaded reader.
 */
void cram_drain_rqueue(cram_fd *fd) {
    cram_container *lc = NULL;

    if (!fd->pool || !fd->rqueue)
        return;

    // drain queue of any in-flight decode jobs
    while (!hts_tpool_process_empty(fd->rqueue)) {
        hts_tpool_result *r = hts_tpool_next_result_wait(fd->rqueue);
        if (!r)
            break;
        cram_decode_job *j = (cram_decode_job *)hts_tpool_result_data(r);
        if (j->c->slice == j->s)
            j->c->slice = NULL;
        if (j->c != lc) {
            if (lc) {
                if (fd->ctr == lc)
                    fd->ctr = NULL;
                if (fd->ctr_mt == lc)
                    fd->ctr_mt = NULL;
                cram_free_container(lc);
            }
            lc = j->c;
        }
        cram_free_slice(j->s);
        hts_tpool_delete_result(r, 1);
    }

    // Also tidy up any pending decode job that we didn't submit to the workers
    // due to the input queue being full.
    if (fd->job_pending) {
        cram_decode_job *j = (cram_decode_job *)fd->job_pending;
        if (j->c->slice == j->s)
            j->c->slice = NULL;
        if (j->c != lc) {
            if (lc) {
                if (fd->ctr == lc)
                    fd->ctr = NULL;
                if (fd->ctr_mt == lc)
                    fd->ctr_mt = NULL;
                cram_free_container(lc);
            }
            lc = j->c;
        }
        cram_free_slice(j->s);
        free(j);
        fd->job_pending = NULL;
    }

    if (lc) {
        if (fd->ctr == lc)
            fd->ctr = NULL;
        if (fd->ctr_mt == lc)
            fd->ctr_mt = NULL;
        cram_free_container(lc);
    }
}
