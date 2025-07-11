/*  tbx.c -- tabix API functions.

    Copyright (C) 2009, 2010, 2012-2015, 2017-2020, 2022-2023 Genome Research Ltd.
    Copyright (C) 2010-2012 Broad Institute.

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

#define HTS_BUILDING_LIBRARY // Enables HTSLIB_EXPORT, see htslib/hts_defs.h
#include <config.h>

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <errno.h>
#include "htslib/tbx.h"
#include "htslib/bgzf.h"
#include "htslib/hts_endian.h"
#include "hts_internal.h"

#include "htslib/khash.h"
KHASH_DECLARE(s2i, kh_cstr_t, int64_t)

HTSLIB_EXPORT
const tbx_conf_t tbx_conf_gff = { 0, 1, 4, 5, '#', 0 };

HTSLIB_EXPORT
const tbx_conf_t tbx_conf_bed = { TBX_UCSC, 1, 2, 3, '#', 0 };

HTSLIB_EXPORT
const tbx_conf_t tbx_conf_psltbl = { TBX_UCSC, 15, 17, 18, '#', 0 };

HTSLIB_EXPORT
const tbx_conf_t tbx_conf_sam = { TBX_SAM, 3, 4, 0, '@', 0 };

HTSLIB_EXPORT
const tbx_conf_t tbx_conf_vcf = { TBX_VCF, 1, 2, 0, '#', 0 };

typedef struct {
    int64_t beg, end;
    char *ss, *se;
    int tid;
} tbx_intv_t;

static inline int get_tid(tbx_t *tbx, const char *ss, int is_add)
{
    khint_t k;
    khash_t(s2i) *d;
    if (tbx->dict == 0) tbx->dict = kh_init(s2i);
    if (!tbx->dict) return -1; // Out of memory
    d = (khash_t(s2i)*)tbx->dict;
    if (is_add) {
        int absent;
        k = kh_put(s2i, d, ss, &absent);
        if (absent < 0) {
            return -1; // Out of memory
        } else if (absent) {
            char *ss_dup = strdup(ss);
            if (ss_dup) {
                kh_key(d, k) = ss_dup;
                kh_val(d, k) = kh_size(d) - 1;
            } else {
                kh_del(s2i, d, k);
                return -1; // Out of memory
            }
        }
    } else k = kh_get(s2i, d, ss);
    return k == kh_end(d)? -1 : kh_val(d, k);
}

int tbx_name2id(tbx_t *tbx, const char *ss)
{
    return get_tid(tbx, ss, 0);
}

int tbx_parse1(const tbx_conf_t *conf, size_t len, char *line, tbx_intv_t *intv)
{
    size_t i, b = 0;
    int id = 1;
    char *s;
    intv->ss = intv->se = 0; intv->beg = intv->end = -1;
    for (i = 0; i <= len; ++i) {
        if (line[i] == '\t' || line[i] == 0) {
            if (id == conf->sc) {
                intv->ss = line + b; intv->se = line + i;
            } else if (id == conf->bc) {
                // here ->beg is 0-based.
                intv->beg = strtoll(line + b, &s, 0);

                if (conf->bc <= conf->ec) // don't overwrite an already set end point
                    intv->end = intv->beg;

                if ( s==line+b ) return -1; // expected int

                if (!(conf->preset&TBX_UCSC))
                    --intv->beg;
                else if (conf->bc <= conf->ec)
                    ++intv->end;

                if (intv->beg < 0) {
                    hts_log_warning("Coordinate <= 0 detected. "
                                    "Did you forget to use the -0 option?");
                    intv->beg = 0;
                }
                if (intv->end < 1) intv->end = 1;
            } else {
                if ((conf->preset&0xffff) == TBX_GENERIC) {
                    if (id == conf->ec)
                    {
                        intv->end = strtoll(line + b, &s, 0);
                        if ( s==line+b ) return -1; // expected int
                    }
                } else if ((conf->preset&0xffff) == TBX_SAM) {
                    if (id == 6) { // CIGAR
                        int l = 0;
                        char *t;
                        for (s = line + b; s < line + i;) {
                            long x = strtol(s, &t, 10);
                            char op = toupper_c(*t);
                            if (op == 'M' || op == 'D' || op == 'N') l += x;
                            s = t + 1;
                        }
                        if (l == 0) l = 1;
                        intv->end = intv->beg + l;
                    }
                } else if ((conf->preset&0xffff) == TBX_VCF) {
                    if (id == 4) {
                        if (b < i) intv->end = intv->beg + (i - b);
                    } else if (id == 8) { // look for "END="
                        int c = line[i];
                        line[i] = 0;
                        s = strstr(line + b, "END=");
                        if (s == line + b) s += 4;
                        else if (s) {
                            s = strstr(line + b, ";END=");
                            if (s) s += 5;
                        }
                        if (s && *s != '.') {
                            long long end = strtoll(s, &s, 0);
                            if (end <= intv->beg) {
                                static int reported = 0;
                                if (!reported) {
                                    int l = intv->ss ? (int) (intv->se - intv->ss) : 0;
                                    hts_log_warning("VCF INFO/END=%lld is smaller than POS at %.*s:%"PRIhts_pos"\n"
                                                    "This tag will be ignored. "
                                                    "Note: only one invalid END tag will be reported.",
                                                    end, l >= 0 ? l : 0,
                                                    intv->ss ? intv->ss : "",
                                                    intv->beg);
                                    reported = 1;
                                }
                            } else {
                                intv->end = end;
                            }
                        }
                        line[i] = c;
                    }
                }
            }
            b = i + 1;
            ++id;
        }
    }
    if (intv->ss == 0 || intv->se == 0 || intv->beg < 0 || intv->end < 0) return -1;
    return 0;
}

static inline int get_intv(tbx_t *tbx, kstring_t *str, tbx_intv_t *intv, int is_add)
{
    if (tbx_parse1(&tbx->conf, str->l, str->s, intv) == 0) {
        int c = *intv->se;
        *intv->se = '\0'; intv->tid = get_tid(tbx, intv->ss, is_add); *intv->se = c;
        if (intv->tid < 0) return -2;  // get_tid out of memory
        return (intv->beg >= 0 && intv->end >= 0)? 0 : -1;
    } else {
        char *type = NULL;
        switch (tbx->conf.preset&0xffff)
        {
            case TBX_SAM: type = "TBX_SAM"; break;
            case TBX_VCF: type = "TBX_VCF"; break;
            case TBX_UCSC: type = "TBX_UCSC"; break;
            default: type = "TBX_GENERIC"; break;
        }
        hts_log_error("Failed to parse %s, was wrong -p [type] used?\nThe offending line was: \"%s\"",
            type, str->s);
        return -1;
    }
}

/*
 * Called by tabix iterator to read the next record.
 * Returns    >=  0 on success
 *               -1 on EOF
 *            <= -2 on error
 */
int tbx_readrec(BGZF *fp, void *tbxv, void *sv, int *tid, hts_pos_t *beg, hts_pos_t *end)
{
    tbx_t *tbx = (tbx_t *) tbxv;
    kstring_t *s = (kstring_t *) sv;
    int ret;
    if ((ret = bgzf_getline(fp, '\n', s)) >= 0) {
        tbx_intv_t intv;
        if (get_intv(tbx, s, &intv, 0) < 0)
            return -2;
        *tid = intv.tid; *beg = intv.beg; *end = intv.end;
    }
    return ret;
}

static int tbx_set_meta(tbx_t *tbx)
{
    int i, l = 0, l_nm;
    uint32_t x[7];
    char **name;
    uint8_t *meta;
    khint_t k;
    khash_t(s2i) *d = (khash_t(s2i)*)tbx->dict;

    memcpy(x, &tbx->conf, 24);
    name = (char**)malloc(sizeof(char*) * kh_size(d));
    if (!name) return -1;
    for (k = kh_begin(d), l = 0; k != kh_end(d); ++k) {
        if (!kh_exist(d, k)) continue;
        name[kh_val(d, k)] = (char*)kh_key(d, k);
        l += strlen(kh_key(d, k)) + 1; // +1 to include '\0'
    }
    l_nm = x[6] = l;
    meta = (uint8_t*)malloc(l_nm + 28);
    if (!meta) { free(name); return -1; }
    if (ed_is_big())
        for (i = 0; i < 7; ++i)
            x[i] = ed_swap_4(x[i]);
    memcpy(meta, x, 28);
    for (l = 28, i = 0; i < (int)kh_size(d); ++i) {
        int x = strlen(name[i]) + 1;
        memcpy(meta + l, name[i], x);
        l += x;
    }
    free(name);
    hts_idx_set_meta(tbx->idx, l, meta, 0);
    return 0;
}

// Minimal effort parser to extract reference length out of VCF header line
// This is used only used to adjust the number of levels if necessary,
// so not a major problem if it doesn't always work.
static void adjust_max_ref_len_vcf(const char *str, int64_t *max_ref_len)
{
    const char *ptr;
    int64_t len;
    if (strncmp(str, "##contig", 8) != 0) return;
    ptr = strstr(str + 8, "length");
    if (!ptr) return;
    for (ptr += 6; *ptr == ' ' || *ptr == '='; ptr++) {}
    len = strtoll(ptr, NULL, 10);
    if (*max_ref_len < len) *max_ref_len = len;
}

// Same for sam files
static void adjust_max_ref_len_sam(const char *str, int64_t *max_ref_len)
{
    const char *ptr;
    int64_t len;
    if (strncmp(str, "@SQ", 3) != 0) return;
    ptr = strstr(str + 3, "\tLN:");
    if (!ptr) return;
    ptr += 4;
    len = strtoll(ptr, NULL, 10);
    if (*max_ref_len < len) *max_ref_len = len;
}

// Adjusts number of levels if not big enough.  This can happen for
// files with very large contigs.
static int adjust_n_lvls(int min_shift, int n_lvls, int64_t max_len)
{
    int64_t s = 1LL << (min_shift + n_lvls * 3);
    max_len += 256;
    for (; max_len > s; ++n_lvls, s <<= 3) {}
    return n_lvls;
}

tbx_t *tbx_index(BGZF *fp, int min_shift, const tbx_conf_t *conf)
{
    tbx_t *tbx;
    kstring_t str;
    int ret, first = 0, n_lvls, fmt;
    int64_t lineno = 0;
    uint64_t last_off = 0;
    tbx_intv_t intv;
    int64_t max_ref_len = 0;

    str.s = 0; str.l = str.m = 0;
    tbx = (tbx_t*)calloc(1, sizeof(tbx_t));
    if (!tbx) return NULL;
    tbx->conf = *conf;
    if (min_shift > 0) n_lvls = (TBX_MAX_SHIFT - min_shift + 2) / 3, fmt = HTS_FMT_CSI;
    else min_shift = 14, n_lvls = 5, fmt = HTS_FMT_TBI;
    while ((ret = bgzf_getline(fp, '\n', &str)) >= 0) {
        ++lineno;
        if (str.s[0] == tbx->conf.meta_char && fmt == HTS_FMT_CSI) {
            switch (tbx->conf.preset) {
                case TBX_SAM:
                    adjust_max_ref_len_sam(str.s, &max_ref_len); break;
                case TBX_VCF:
                    adjust_max_ref_len_vcf(str.s, &max_ref_len); break;
                default:
                    break;
            }
        }
        if (lineno <= tbx->conf.line_skip || str.s[0] == tbx->conf.meta_char) {
            last_off = bgzf_tell(fp);
            continue;
        }
        if (first == 0) {
            if (fmt == HTS_FMT_CSI) {
                if (!max_ref_len)
                    max_ref_len = (int64_t)100*1024*1024*1024; // 100G default
                n_lvls = adjust_n_lvls(min_shift, n_lvls, max_ref_len);
            }
            tbx->idx = hts_idx_init(0, fmt, last_off, min_shift, n_lvls);
            if (!tbx->idx) goto fail;
            first = 1;
        }
        ret = get_intv(tbx, &str, &intv, 1);
        if (ret < -1) goto fail;  // Out of memory
        if (ret < 0) continue; // Skip unparsable lines
        if (hts_idx_push(tbx->idx, intv.tid, intv.beg, intv.end,
                         bgzf_tell(fp), 1) < 0) {
            goto fail;
        }
    }
    if (ret < -1) goto fail;
    if ( !tbx->idx ) tbx->idx = hts_idx_init(0, fmt, last_off, min_shift, n_lvls);   // empty file
    if (!tbx->idx) goto fail;
    if ( !tbx->dict ) tbx->dict = kh_init(s2i);
    if (!tbx->dict) goto fail;
    if (hts_idx_finish(tbx->idx, bgzf_tell(fp)) != 0) goto fail;
    if (tbx_set_meta(tbx) != 0) goto fail;
    free(str.s);
    return tbx;

 fail:
    free(str.s);
    tbx_destroy(tbx);
    return NULL;
}

void tbx_destroy(tbx_t *tbx)
{
    khash_t(s2i) *d = (khash_t(s2i)*)tbx->dict;
    if (d != NULL)
    {
        khint_t k;
        for (k = kh_begin(d); k != kh_end(d); ++k)
            if (kh_exist(d, k)) free((char*)kh_key(d, k));
    }
    hts_idx_destroy(tbx->idx);
    kh_destroy(s2i, d);
    free(tbx);
}

static tbx_t *index_load(const char *fn, const char *fnidx, int flags)
{
    tbx_t *tbx;
    uint8_t *meta;
    char *nm, *p;
    uint32_t l_meta, l_nm;
    tbx = (tbx_t*)calloc(1, sizeof(tbx_t));
    if (!tbx)
        return NULL;
    tbx->idx = hts_idx_load3(fn, fnidx, HTS_FMT_TBI, flags);
    if ( !tbx->idx )
    {
        free(tbx);
        return NULL;
    }
    meta = hts_idx_get_meta(tbx->idx, &l_meta);
    if ( !meta || l_meta < 28) goto invalid;

    tbx->conf.preset = le_to_i32(&meta[0]);
    tbx->conf.sc = le_to_i32(&meta[4]);
    tbx->conf.bc = le_to_i32(&meta[8]);
    tbx->conf.ec = le_to_i32(&meta[12]);
    tbx->conf.meta_char = le_to_i32(&meta[16]);
    tbx->conf.line_skip = le_to_i32(&meta[20]);
    l_nm = le_to_u32(&meta[24]);
    if (l_nm > l_meta - 28) goto invalid;

    p = nm = (char*)meta + 28;
    // This assumes meta is NUL-terminated, so we can merrily strlen away.
    // hts_idx_load_local() assures this for us by adding a NUL on the end
    // of whatever it reads.
    for (; p - nm < l_nm; p += strlen(p) + 1) {
        if (get_tid(tbx, p, 1) < 0) {
            hts_log_error("%s", strerror(errno));
            goto fail;
        }
    }
    return tbx;

 invalid:
    hts_log_error("Invalid index header for %s", fnidx ? fnidx : fn);

 fail:
    tbx_destroy(tbx);
    return NULL;
}

tbx_t *tbx_index_load3(const char *fn, const char *fnidx, int flags)
{
    return index_load(fn, fnidx, flags);
}

tbx_t *tbx_index_load2(const char *fn, const char *fnidx)
{
    return index_load(fn, fnidx, 1);
}

tbx_t *tbx_index_load(const char *fn)
{
    return index_load(fn, NULL, 1);
}

const char **tbx_seqnames(tbx_t *tbx, int *n)
{
    khash_t(s2i) *d = (khash_t(s2i)*)tbx->dict;
    if (d == NULL)
    {
        *n = 0;
        return calloc(1, sizeof(char *));
    }
    int tid, m = kh_size(d);
    const char **names = (const char**) calloc(m,sizeof(const char*));
    khint_t k;
    if (!names) {
        *n = 0;
        return NULL;
    }
    for (k=kh_begin(d); k<kh_end(d); k++)
    {
        if ( !kh_exist(d,k) ) continue;
        tid = kh_val(d,k);
        assert( tid<m );
        names[tid] = kh_key(d,k);
    }
    // sanity check: there should be no gaps
    for (tid=0; tid<m; tid++)
        assert(names[tid]);
    *n = m;
    return names;
}

