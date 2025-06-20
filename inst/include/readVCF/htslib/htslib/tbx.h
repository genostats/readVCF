/// @file htslib/tbx.h
/// Tabix API functions.
/*
    Copyright (C) 2009, 2012-2015, 2019 Genome Research Ltd.
    Copyright (C) 2010, 2012 Broad Institute.

    Author: Heng Li <lh3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#ifndef HTSLIB_TBX_H
#define HTSLIB_TBX_H

#include "hts.h"

#ifdef __cplusplus
extern "C" {
#endif

#define TBX_MAX_SHIFT 31

#define TBX_GENERIC 0
#define TBX_SAM     1
#define TBX_VCF     2
#define TBX_UCSC    0x10000

typedef struct tbx_conf_t {
    int32_t preset;
    int32_t sc, bc, ec; // seq col., beg col. and end col.
    int32_t meta_char, line_skip;
} tbx_conf_t;

typedef struct tbx_t {
    tbx_conf_t conf;
    hts_idx_t *idx;
    void *dict;
} tbx_t;

HTSLIB_EXPORT
extern const tbx_conf_t tbx_conf_gff, tbx_conf_bed, tbx_conf_psltbl, tbx_conf_sam, tbx_conf_vcf;

    HTSLIB_EXPORT
    int tbx_readrec(BGZF *fp, void *tbxv, void *sv, int *tid, hts_pos_t *beg, hts_pos_t *end);

/// Build an index of the lines in a BGZF-compressed file
/** The index struct returned by a successful call should be freed
    via tbx_destroy() when it is no longer needed.
*/
    HTSLIB_EXPORT
    tbx_t *tbx_index(BGZF *fp, int min_shift, const tbx_conf_t *conf);
/*
 * All tbx_index_build* methods return: 0 (success), -1 (general failure) or -2 (compression not BGZF)
 */
    HTSLIB_EXPORT
    int tbx_index_build(const char *fn, int min_shift, const tbx_conf_t *conf);

    HTSLIB_EXPORT
    int tbx_index_build2(const char *fn, const char *fnidx, int min_shift, const tbx_conf_t *conf);

    HTSLIB_EXPORT
    int tbx_index_build3(const char *fn, const char *fnidx, int min_shift, int n_threads, const tbx_conf_t *conf);


/// Load or stream a .tbi or .csi index
/** @param fn     Name of the data file corresponding to the index

    Equivalent to tbx_index_load3(fn, NULL, HTS_IDX_SAVE_REMOTE);
*/
    HTSLIB_EXPORT
    tbx_t *tbx_index_load(const char *fn);

/// Load or stream a .tbi or .csi index
/** @param fn     Name of the data file corresponding to the index
    @param fnidx  Name of the indexed file
    @return The index, or NULL if an error occurred

    If @p fnidx is NULL, the index name will be derived from @p fn.

    Equivalent to tbx_index_load3(fn, fnidx, HTS_IDX_SAVE_REMOTE);
*/
    HTSLIB_EXPORT
    tbx_t *tbx_index_load2(const char *fn, const char *fnidx);

/// Load or stream a .tbi or .csi index
/** @param fn     Name of the data file corresponding to the index
    @param fnidx  Name of the indexed file
    @param flags  Flags to alter behaviour (see description)
    @return The index, or NULL if an error occurred

    If @p fnidx is NULL, the index name will be derived from @p fn.

    The @p flags parameter can be set to a combination of the following
    values:

        HTS_IDX_SAVE_REMOTE   Save a local copy of any remote indexes
        HTS_IDX_SILENT_FAIL   Fail silently if the index is not present

    The index struct returned by a successful call should be freed
    via tbx_destroy() when it is no longer needed.
*/
    HTSLIB_EXPORT
    tbx_t *tbx_index_load3(const char *fn, const char *fnidx, int flags);

    HTSLIB_EXPORT
    const char **tbx_seqnames(tbx_t *tbx, int *n);  // free the array but not the values

    HTSLIB_EXPORT
    void tbx_destroy(tbx_t *tbx);

#ifdef __cplusplus
}
#endif

#endif
