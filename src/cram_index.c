/*
Copyright (c) 2013-2020, 2023 Genome Research Ltd.
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
 * The index is a gzipped tab-delimited text file with one line per slice.
 * The columns are:
 * 1: reference number (0 to N-1, as per BAM ref_id)
 * 2: reference position of 1st read in slice (1..?)
 * 3: number of reads in slice
 * 4: offset of container start (relative to end of SAM header, so 1st
 *    container is offset 0).
 * 5: slice number within container (ie which landmark).
 *
 * In memory, we hold this in a nested containment list. Each list element is
 * a cram_index struct. Each element in turn can contain its own list of
 * cram_index structs.
 *
 * Any start..end range which is entirely contained within another (and
 * earlier as it is sorted) range will be held within it. This ensures that
 * the outer list will never have containments and we can safely do a
 * binary search to find the first range which overlaps any given coordinate.
 */

#define HTS_BUILDING_LIBRARY // Enables HTSLIB_EXPORT, see htslib/hts_defs.h
#include <config.h>

#include <stdio.h>
#include <errno.h>
#include <assert.h>
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>

#include "htslib/bgzf.h"
#include "htslib/hfile.h"
#include "hts_internal.h"
#include "cram/cram.h"
#include "cram/os.h"

static void cram_index_free_recurse(cram_index *e) {
    if (e->e) {
        int i;
        for (i = 0; i < e->nslice; i++) {
            cram_index_free_recurse(&e->e[i]);
        }
        free(e->e);
    }
}

void cram_index_free(cram_fd *fd) {
    int i;

    if (!fd->index)
        return;

    for (i = 0; i < fd->index_sz; i++) {
        cram_index_free_recurse(&fd->index[i]);
    }
    free(fd->index);

    fd->index = NULL;
}

/*
 * Searches the index for the first slice overlapping a reference ID
 * and position, or one immediately preceding it if none is found in
 * the index to overlap this position. (Our index may have missing
 * entries, but we require at least one per reference.)
 *
 * If the index finds multiple slices overlapping this position we
 * return the first one only. Subsequent calls should specify
 * "from" as the last slice we checked to find the next one. Otherwise
 * set "from" to be NULL to find the first one.
 *
 * Refid can also be any of the special HTS_IDX_ values.
 * For backwards compatibility, refid -1 is equivalent to HTS_IDX_NOCOOR.
 *
 * Returns the cram_index pointer on success
 *         NULL on failure
 */
cram_index *cram_index_query(cram_fd *fd, int refid, hts_pos_t pos,
                             cram_index *from) {
    int i, j, k;
    cram_index *e;

    if (from) {
        // Continue from a previous search.
        // We switch to just scanning the linked list, as the nested
        // lists are typically short.
        e = from->e_next;
        if (e && e->refid == refid && e->start <= pos)
            return e;
        else
            return NULL;
    }

    switch(refid) {
    case HTS_IDX_NONE:
    case HTS_IDX_REST:
        // fail, or already there, dealt with elsewhere.
        return NULL;

    case HTS_IDX_NOCOOR:
        refid = -1;
        pos = 0;
        break;

    case HTS_IDX_START: {
        int64_t min_idx = INT64_MAX;
        for (i = 0, j = -1; i < fd->index_sz; i++) {
            if (fd->index[i].e && fd->index[i].e[0].offset < min_idx) {
                min_idx = fd->index[i].e[0].offset;
                j = i;
            }
        }
        if (j < 0)
            return NULL;
        return fd->index[j].e;
    }

    default:
        if (refid < HTS_IDX_NONE || refid+1 >= fd->index_sz)
            return NULL;
    }

    from = &fd->index[refid+1];

    // Ref with nothing aligned against it.
    if (!from->e)
        return NULL;

    // This sequence is covered by the index, so binary search to find
    // the optimal starting block.
    i = 0, j = fd->index[refid+1].nslice-1;
    for (k = j/2; k != i; k = (j-i)/2 + i) {
        if (from->e[k].refid > refid) {
            j = k;
            continue;
        }

        if (from->e[k].refid < refid) {
            i = k;
            continue;
        }

        if (from->e[k].start >= pos) {
            j = k;
            continue;
        }

        if (from->e[k].start < pos) {
            i = k;
            continue;
        }
    }
    // i==j or i==j-1. Check if j is better.
    if (j >= 0 && from->e[j].start < pos && from->e[j].refid == refid)
        i = j;

    /* The above found *a* bin overlapping, but not necessarily the first */
    while (i > 0 && from->e[i-1].end >= pos)
        i--;

    /* We may be one bin before the optimum, so check */
    while (i+1 < from->nslice &&
           (from->e[i].refid < refid ||
            from->e[i].end < pos))
        i++;

    e = &from->e[i];

    return e;
}

// Return the index entry for last slice on a specific reference.
cram_index *cram_index_last(cram_fd *fd, int refid, cram_index *from) {
    int slice;

    if (refid+1 < 0 || refid+1 >= fd->index_sz)
        return NULL;

    if (!from)
        from = &fd->index[refid+1];

    // Ref with nothing aligned against it.
    if (!from->e)
        return NULL;

    slice = fd->index[refid+1].nslice - 1;

    // e is the last entry in the nested containment list, but it may
    // contain further slices within it.
    cram_index *e = &from->e[slice];
    while (e->e_next)
        e = e->e_next;

    return e;
}

/*
 * Find the last container overlapping pos 'end', and the file offset of
 * its end (equivalent to the start offset of the container following it).
 */
cram_index *cram_index_query_last(cram_fd *fd, int refid, hts_pos_t end) {
    cram_index *e = NULL, *prev_e;
    do {
        prev_e = e;
        e = cram_index_query(fd, refid, end, prev_e);
    } while (e);

    if (!prev_e)
        return NULL;
    e = prev_e;

    // Note: offset of e and e->e_next may be the same if we're using a
    // multi-ref container where a single container generates multiple
    // index entries.
    //
    // We need to keep iterating until offset differs in order to find
    // the genuine file offset for the end of container.
    do {
        prev_e = e;
        e = e->e_next;
    } while (e && e->offset == prev_e->offset);

    return prev_e;
}

/*
 * Skips to a container overlapping the start coordinate listed in
 * cram_range.
 *
 * In theory we call cram_index_query multiple times, once per slice
 * overlapping the range. However slices may be absent from the index
 * which makes this problematic. Instead we find the left-most slice
 * and then read from then on, skipping decoding of slices and/or
 * whole containers when they don't overlap the specified cram_range.
 *
 * This function also updates the cram_fd range field.
 *
 * Returns 0 on success
 *        -1 on general failure
 *        -2 on no-data (empty chromosome)
 */
int cram_seek_to_refpos(cram_fd *fd, cram_range *r) {
    int ret = 0;
    cram_index *e;

    if (r->refid == HTS_IDX_NONE) {
        ret = -2; goto err;
    }

    // Ideally use an index, so see if we have one.
    if ((e = cram_index_query(fd, r->refid, r->start, NULL))) {
        if (0 != cram_seek(fd, e->offset, SEEK_SET)) {
            if (0 != cram_seek(fd, e->offset - fd->first_container, SEEK_CUR)) {
                ret = -1; goto err;
            }
        }
    } else {
        // Absent from index, but this most likely means it simply has no data.
        ret = -2; goto err;
    }

    pthread_mutex_lock(&fd->range_lock);
    fd->range = *r;
    if (r->refid == HTS_IDX_NOCOOR) {
        fd->range.refid = -1;
        fd->range.start = 0;
    } else if (r->refid == HTS_IDX_START || r->refid == HTS_IDX_REST) {
        fd->range.refid = -2; // special case in cram_next_slice
    }
    pthread_mutex_unlock(&fd->range_lock);

    if (fd->ctr) {
        cram_free_container(fd->ctr);
        if (fd->ctr_mt && fd->ctr_mt != fd->ctr)
            cram_free_container(fd->ctr_mt);
        fd->ctr = NULL;
        fd->ctr_mt = NULL;
        fd->ooc = 0;
        fd->eof = 0;
    }

    return 0;

 err:
    // It's unlikely fd->range will be accessed after EOF or error,
    // but this maintains identical behaviour to the previous code.
    pthread_mutex_lock(&fd->range_lock);
    fd->range = *r;
    pthread_mutex_unlock(&fd->range_lock);
    return ret;
}