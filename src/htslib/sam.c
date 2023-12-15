/*  sam.c -- SAM and BAM file I/O and manipulation.

    Copyright (C) 2008-2010, 2012-2023 Genome Research Ltd.
    Copyright (C) 2010, 2012, 2013 Broad Institute.

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

#include <strings.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <zlib.h>
#include <assert.h>
#include <signal.h>
#include <inttypes.h>
#include <unistd.h>

#ifdef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
#include "fuzz_settings.h"
#endif

// Suppress deprecation message for cigar_tab, which we initialise
#include "htslib/hts_defs.h"
#undef HTS_DEPRECATED
#define HTS_DEPRECATED(message)

#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include "cram/cram.h"
#include "hts_internal.h"
#include "sam_internal.h"
#include "htslib/hfile.h"
#include "htslib/hts_endian.h"
#include "htslib/hts_expr.h"
#include "header.h"

#include "htslib/khash.h"
KHASH_DECLARE(s2i, kh_cstr_t, int64_t)
KHASH_SET_INIT_INT(tag)

#ifndef EFTYPE
#define EFTYPE ENOEXEC
#endif
#ifndef EOVERFLOW
#define EOVERFLOW ERANGE
#endif

/**********************
 *** BAM header I/O ***
 **********************/

HTSLIB_EXPORT
const int8_t bam_cigar_table[256] = {
    // 0 .. 47
    -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,
    -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,
    -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,

    // 48 .. 63  (including =)
    -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, BAM_CEQUAL, -1, -1,

    // 64 .. 79  (including MIDNHB)
    -1, -1, BAM_CBACK, -1,  BAM_CDEL, -1, -1, -1,
        BAM_CHARD_CLIP, BAM_CINS, -1, -1,  -1, BAM_CMATCH, BAM_CREF_SKIP, -1,

    // 80 .. 95  (including SPX)
    BAM_CPAD, -1, -1, BAM_CSOFT_CLIP,  -1, -1, -1, -1,
        BAM_CDIFF, -1, -1, -1,  -1, -1, -1, -1,

    // 96 .. 127
    -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,
    -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,

    // 128 .. 255
    -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,
    -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,
    -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,
    -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,
    -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,
    -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,
    -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,
    -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1
};

sam_hdr_t *sam_hdr_init()
{
    sam_hdr_t *bh = (sam_hdr_t*)calloc(1, sizeof(sam_hdr_t));
    if (bh == NULL) return NULL;

    bh->cigar_tab = bam_cigar_table;
    return bh;
}

HTSLIB_EXPORT
void sam_hdr_destroy(sam_hdr_t *bh)
{
    int32_t i;

    if (bh == NULL) return;

    if (bh->ref_count > 0) {
        --bh->ref_count;
        return;
    }

    if (bh->target_name) {
        for (i = 0; i < bh->n_targets; ++i)
            free(bh->target_name[i]);
        free(bh->target_name);
        free(bh->target_len);
    }
    free(bh->text);
    if (bh->hrecs)
        sam_hrecs_free(bh->hrecs);
    if (bh->sdict)
        kh_destroy(s2i, (khash_t(s2i) *) bh->sdict);
    free(bh);
}

// Copy the sam_hdr_t::sdict hash, used to store the real lengths of long
// references before sam_hdr_t::hrecs is populated
int sam_hdr_dup_sdict(const sam_hdr_t *h0, sam_hdr_t *h)
{
    const khash_t(s2i) *src_long_refs = (khash_t(s2i) *) h0->sdict;
    khash_t(s2i) *dest_long_refs = kh_init(s2i);
    int i;
    if (!dest_long_refs) return -1;

    for (i = 0; i < h->n_targets; i++) {
        int ret;
        khiter_t ksrc, kdest;
        if (h->target_len[i] < UINT32_MAX) continue;
        ksrc = kh_get(s2i, src_long_refs, h->target_name[i]);
        if (ksrc == kh_end(src_long_refs)) continue;
        kdest = kh_put(s2i, dest_long_refs, h->target_name[i], &ret);
        if (ret < 0) {
            kh_destroy(s2i, dest_long_refs);
            return -1;
        }
        kh_val(dest_long_refs, kdest) = kh_val(src_long_refs, ksrc);
    }

    h->sdict = dest_long_refs;
    return 0;
}

sam_hdr_t *sam_hdr_dup(const sam_hdr_t *h0)
{
    if (h0 == NULL) return NULL;
    sam_hdr_t *h;
    if ((h = sam_hdr_init()) == NULL) return NULL;
    // copy the simple data
    h->n_targets = 0;
    h->ignore_sam_err = h0->ignore_sam_err;
    h->l_text = 0;

    // Then the pointery stuff

    if (!h0->hrecs) {
        h->target_len = (uint32_t*)calloc(h0->n_targets, sizeof(uint32_t));
        if (!h->target_len) goto fail;
        h->target_name = (char**)calloc(h0->n_targets, sizeof(char*));
        if (!h->target_name) goto fail;

        int i;
        for (i = 0; i < h0->n_targets; ++i) {
            h->target_len[i] = h0->target_len[i];
            h->target_name[i] = strdup(h0->target_name[i]);
            if (!h->target_name[i]) break;
        }
        h->n_targets = i;
        if (i < h0->n_targets) goto fail;

        if (h0->sdict) {
            if (sam_hdr_dup_sdict(h0, h) < 0) goto fail;
        }
    }

    if (h0->hrecs) {
        kstring_t tmp = { 0, 0, NULL };
        if (sam_hrecs_rebuild_text(h0->hrecs, &tmp) != 0) {
            free(ks_release(&tmp));
            goto fail;
        }

        h->l_text = tmp.l;
        h->text   = ks_release(&tmp);

        if (sam_hdr_update_target_arrays(h, h0->hrecs, 0) != 0)
            goto fail;
    } else {
        h->l_text = h0->l_text;
        h->text = malloc(h->l_text + 1);
        if (!h->text) goto fail;
        memcpy(h->text, h0->text, h->l_text);
        h->text[h->l_text] = '\0';
    }

    return h;

 fail:
    sam_hdr_destroy(h);
    return NULL;
}

sam_hdr_t *bam_hdr_read(BGZF *fp)
{
    sam_hdr_t *h;
    uint8_t buf[4];
    int magic_len, has_EOF;
    int32_t i, name_len, num_names = 0;
    size_t bufsize;
    ssize_t bytes;
    // check EOF
    has_EOF = bgzf_check_EOF(fp);
    if (has_EOF < 0) {
        perror("[W::bam_hdr_read] bgzf_check_EOF");
    } else if (has_EOF == 0) {
        hts_log_warning("EOF marker is absent. The input is probably truncated");
    }
    // read "BAM1"
    magic_len = bgzf_read(fp, buf, 4);
    if (magic_len != 4 || memcmp(buf, "BAM\1", 4)) {
        hts_log_error("Invalid BAM binary header");
        return 0;
    }
    h = sam_hdr_init();
    if (!h) goto nomem;

    // read plain text and the number of reference sequences
    bytes = bgzf_read(fp, buf, 4);
    if (bytes != 4) goto read_err;
    h->l_text = le_to_u32(buf);

    bufsize = h->l_text + 1;
    if (bufsize < h->l_text) goto nomem; // so large that adding 1 overflowed
#ifdef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
    if (bufsize > FUZZ_ALLOC_LIMIT) goto nomem;
#endif
    h->text = (char*)malloc(bufsize);
    if (!h->text) goto nomem;
    h->text[h->l_text] = 0; // make sure it is NULL terminated
    bytes = bgzf_read(fp, h->text, h->l_text);
    if (bytes != h->l_text) goto read_err;

    bytes = bgzf_read(fp, &h->n_targets, 4);
    if (bytes != 4) goto read_err;
    if (fp->is_be) ed_swap_4p(&h->n_targets);

    if (h->n_targets < 0) goto invalid;

    // read reference sequence names and lengths
#ifdef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
    if (h->n_targets > (FUZZ_ALLOC_LIMIT - bufsize)/(sizeof(char*)+sizeof(uint32_t)))
        goto nomem;
#endif
    if (h->n_targets > 0) {
        h->target_name = (char**)calloc(h->n_targets, sizeof(char*));
        if (!h->target_name) goto nomem;
        h->target_len = (uint32_t*)calloc(h->n_targets, sizeof(uint32_t));
        if (!h->target_len) goto nomem;
    }
    else {
        h->target_name = NULL;
        h->target_len = NULL;
    }

    for (i = 0; i != h->n_targets; ++i) {
        bytes = bgzf_read(fp, &name_len, 4);
        if (bytes != 4) goto read_err;
        if (fp->is_be) ed_swap_4p(&name_len);
        if (name_len <= 0) goto invalid;

        h->target_name[i] = (char*)malloc(name_len);
        if (!h->target_name[i]) goto nomem;
        num_names++;

        bytes = bgzf_read(fp, h->target_name[i], name_len);
        if (bytes != name_len) goto read_err;

        if (h->target_name[i][name_len - 1] != '\0') {
            /* Fix missing NUL-termination.  Is this being too nice?
               We could alternatively bail out with an error. */
            char *new_name;
            if (name_len == INT32_MAX) goto invalid;
            new_name = realloc(h->target_name[i], name_len + 1);
            if (new_name == NULL) goto nomem;
            h->target_name[i] = new_name;
            h->target_name[i][name_len] = '\0';
        }

        bytes = bgzf_read(fp, &h->target_len[i], 4);
        if (bytes != 4) goto read_err;
        if (fp->is_be) ed_swap_4p(&h->target_len[i]);
    }
    return h;

 nomem:
    hts_log_error("Out of memory");
    goto clean;

 read_err:
    if (bytes < 0) {
        hts_log_error("Error reading BGZF stream");
    } else {
        hts_log_error("Truncated BAM header");
    }
    goto clean;

 invalid:
    hts_log_error("Invalid BAM binary header");

 clean:
    if (h != NULL) {
        h->n_targets = num_names; // ensure we free only allocated target_names
        sam_hdr_destroy(h);
    }
    return NULL;
}

/*************************
 *** BAM alignment I/O ***
 *************************/

int sam_realloc_bam_data(bam1_t *b, size_t desired)
{
    uint32_t new_m_data;
    uint8_t *new_data;
    new_m_data = desired;
    kroundup32(new_m_data);
    if (new_m_data < desired) {
        errno = ENOMEM; // Not strictly true but we can't store the size
        return -1;
    }
#ifdef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
    if (new_m_data > FUZZ_ALLOC_LIMIT) {
        errno = ENOMEM;
        return -1;
    }
#endif
    if ((bam_get_mempolicy(b) & BAM_USER_OWNS_DATA) == 0) {
        new_data = realloc(b->data, new_m_data);
    } else {
        if ((new_data = malloc(new_m_data)) != NULL) {
            if (b->l_data > 0)
                memcpy(new_data, b->data,
                       b->l_data < b->m_data ? b->l_data : b->m_data);
            bam_set_mempolicy(b, bam_get_mempolicy(b) & (~BAM_USER_OWNS_DATA));
        }
    }
    if (!new_data) return -1;
    b->data = new_data;
    b->m_data = new_m_data;
    return 0;
}

void bam_destroy1(bam1_t *b)
{
    if (b == 0) return;
    if ((bam_get_mempolicy(b) & BAM_USER_OWNS_DATA) == 0) {
        free(b->data);
        if ((bam_get_mempolicy(b) & BAM_USER_OWNS_STRUCT) != 0) {
            // In case of reuse
            b->data = NULL;
            b->m_data = 0;
            b->l_data = 0;
        }
    }

    if ((bam_get_mempolicy(b) & BAM_USER_OWNS_STRUCT) == 0)
        free(b);
}

bam1_t *bam_copy1(bam1_t *bdst, const bam1_t *bsrc)
{
    if (realloc_bam_data(bdst, bsrc->l_data) < 0) return NULL;
    memcpy(bdst->data, bsrc->data, bsrc->l_data); // copy var-len data
    memcpy(&bdst->core, &bsrc->core, sizeof(bsrc->core)); // copy the rest
    bdst->l_data = bsrc->l_data;
    bdst->id = bsrc->id;
    return bdst;
}

static void bam_cigar2rqlens(int n_cigar, const uint32_t *cigar,
                             hts_pos_t *rlen, hts_pos_t *qlen)
{
    int k;
    *rlen = *qlen = 0;
    for (k = 0; k < n_cigar; ++k) {
        int type = bam_cigar_type(bam_cigar_op(cigar[k]));
        int len = bam_cigar_oplen(cigar[k]);
        if (type & 1) *qlen += len;
        if (type & 2) *rlen += len;
    }
}

static int subtract_check_underflow(size_t length, size_t *limit)
{
    if (length <= *limit) {
        *limit -= length;
        return 0;
    }

    return -1;
}

int bam_set1(bam1_t *bam,
             size_t l_qname, const char *qname,
             uint16_t flag, int32_t tid, hts_pos_t pos, uint8_t mapq,
             size_t n_cigar, const uint32_t *cigar,
             int32_t mtid, hts_pos_t mpos, hts_pos_t isize,
             size_t l_seq, const char *seq, const char *qual,
             size_t l_aux)
{
    // use a default qname "*" if none is provided
    if (l_qname == 0) {
        l_qname = 1;
        qname = "*";
    }

    // note: the qname is stored nul terminated and padded as described in the
    // documentation for the bam1_t struct.
    size_t qname_nuls = 4 - l_qname % 4;

    // the aligment length, needed for bam_reg2bin(), is calculated as in bam_endpos().
    // can't use bam_endpos() directly as some fields not yet set up.
    hts_pos_t rlen = 0, qlen = 0;
    if (!(flag & BAM_FUNMAP)) {
        bam_cigar2rqlens((int)n_cigar, cigar, &rlen, &qlen);
    }
    if (rlen == 0) {
        rlen = 1;
    }

    // validate parameters
    if (l_qname > 254) {
        hts_log_error("Query name too long");
        errno = EINVAL;
        return -1;
    }
    if (HTS_POS_MAX - rlen <= pos) {
        hts_log_error("Read ends beyond highest supported position");
        errno = EINVAL;
        return -1;
    }
    if (!(flag & BAM_FUNMAP) && l_seq > 0 && n_cigar == 0) {
        hts_log_error("Mapped query must have a CIGAR");
        errno = EINVAL;
        return -1;
    }
    if (!(flag & BAM_FUNMAP) && l_seq > 0 && l_seq != qlen) {
        hts_log_error("CIGAR and query sequence are of different length");
        errno = EINVAL;
        return -1;
    }

    size_t limit = INT32_MAX;
    int u = subtract_check_underflow(l_qname + qname_nuls, &limit);
    u    += subtract_check_underflow(n_cigar * 4, &limit);
    u    += subtract_check_underflow((l_seq + 1) / 2, &limit);
    u    += subtract_check_underflow(l_seq, &limit);
    u    += subtract_check_underflow(l_aux, &limit);
    if (u != 0) {
        hts_log_error("Size overflow");
        errno = EINVAL;
        return -1;
    }

    // re-allocate the data buffer as needed.
    size_t data_len = l_qname + qname_nuls + n_cigar * 4 + (l_seq + 1) / 2 + l_seq;
    if (realloc_bam_data(bam, data_len + l_aux) < 0) {
        return -1;
    }

    bam->l_data = (int)data_len;
    bam->core.pos = pos;
    bam->core.tid = tid;
    bam->core.bin = bam_reg2bin(pos, pos + rlen);
    bam->core.qual = mapq;
    bam->core.l_extranul = (uint8_t)(qname_nuls - 1);
    bam->core.flag = flag;
    bam->core.l_qname = (uint16_t)(l_qname + qname_nuls);
    bam->core.n_cigar = (uint32_t)n_cigar;
    bam->core.l_qseq = (int32_t)l_seq;
    bam->core.mtid = mtid;
    bam->core.mpos = mpos;
    bam->core.isize = isize;

    uint8_t *cp = bam->data;
    strncpy((char *)cp, qname, l_qname);
    int i;
    for (i = 0; i < qname_nuls; i++) {
        cp[l_qname + i] = '\0';
    }
    cp += l_qname + qname_nuls;

    if (n_cigar > 0) {
        memcpy(cp, cigar, n_cigar * 4);
    }
    cp += n_cigar * 4;

#define NN 16
    const uint8_t *useq = (uint8_t *)seq;
    for (i = 0; i + NN < l_seq; i += NN) {
        int j;
        const uint8_t *u2 = useq+i;
        for (j = 0; j < NN/2; j++)
            cp[j] = (seq_nt16_table[u2[j*2]]<<4) | seq_nt16_table[u2[j*2+1]];
        cp += NN/2;
    }
    for (; i + 1 < l_seq; i += 2) {
        *cp++ = (seq_nt16_table[useq[i]] << 4) | seq_nt16_table[useq[i + 1]];
    }

    for (; i < l_seq; i++) {
        *cp++ = seq_nt16_table[(unsigned char)seq[i]] << 4;
    }

    if (qual) {
        memcpy(cp, qual, l_seq);
    }
    else {
        memset(cp, '\xff', l_seq);
    }

    return (int)data_len;
}

hts_pos_t bam_cigar2qlen(int n_cigar, const uint32_t *cigar)
{
    int k;
    hts_pos_t l;
    for (k = l = 0; k < n_cigar; ++k)
        if (bam_cigar_type(bam_cigar_op(cigar[k]))&1)
            l += bam_cigar_oplen(cigar[k]);
    return l;
}

hts_pos_t bam_cigar2rlen(int n_cigar, const uint32_t *cigar)
{
    int k;
    hts_pos_t l;
    for (k = l = 0; k < n_cigar; ++k)
        if (bam_cigar_type(bam_cigar_op(cigar[k]))&2)
            l += bam_cigar_oplen(cigar[k]);
    return l;
}

hts_pos_t bam_endpos(const bam1_t *b)
{
    hts_pos_t rlen = (b->core.flag & BAM_FUNMAP)? 0 : bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b));
    if (rlen == 0) rlen = 1;
    return b->core.pos + rlen;
}

static int bam_tag2cigar(bam1_t *b, int recal_bin, int give_warning) // return 0 if CIGAR is untouched; 1 if CIGAR is updated with CG
{
    bam1_core_t *c = &b->core;
    uint32_t cigar_st, n_cigar4, CG_st, CG_en, ori_len = b->l_data, *cigar0, CG_len, fake_bytes;
    uint8_t *CG;

    // test where there is a real CIGAR in the CG tag to move
    if (c->n_cigar == 0 || c->tid < 0 || c->pos < 0) return 0;
    cigar0 = bam_get_cigar(b);
    if (bam_cigar_op(cigar0[0]) != BAM_CSOFT_CLIP || bam_cigar_oplen(cigar0[0]) != c->l_qseq) return 0;
    fake_bytes = c->n_cigar * 4;
    int saved_errno = errno;
    CG = bam_aux_get(b, "CG");
    if (!CG) {
        if (errno != ENOENT) return -1;  // Bad aux data
        errno = saved_errno; // restore errno on expected no-CG-tag case
        return 0;
    }
    if (CG[0] != 'B' || !(CG[1] == 'I' || CG[1] == 'i'))
        return 0; // not of type B,I
    CG_len = le_to_u32(CG + 2);
    if (CG_len < c->n_cigar || CG_len >= 1U<<29) return 0; // don't move if the real CIGAR length is shorter than the fake cigar length

    // move from the CG tag to the right position
    cigar_st = (uint8_t*)cigar0 - b->data;
    c->n_cigar = CG_len;
    n_cigar4 = c->n_cigar * 4;
    CG_st = CG - b->data - 2;
    CG_en = CG_st + 8 + n_cigar4;
    if (possibly_expand_bam_data(b, n_cigar4 - fake_bytes) < 0) return -1;
    b->l_data = b->l_data - fake_bytes + n_cigar4; // we need c->n_cigar-fake_bytes bytes to swap CIGAR to the right place
    memmove(b->data + cigar_st + n_cigar4, b->data + cigar_st + fake_bytes, ori_len - (cigar_st + fake_bytes)); // insert c->n_cigar-fake_bytes empty space to make room
    memcpy(b->data + cigar_st, b->data + (n_cigar4 - fake_bytes) + CG_st + 8, n_cigar4); // copy the real CIGAR to the right place; -fake_bytes for the fake CIGAR
    if (ori_len > CG_en) // move data after the CG tag
        memmove(b->data + CG_st + n_cigar4 - fake_bytes, b->data + CG_en + n_cigar4 - fake_bytes, ori_len - CG_en);
    b->l_data -= n_cigar4 + 8; // 8: CGBI (4 bytes) and CGBI length (4)
    if (recal_bin)
        b->core.bin = hts_reg2bin(b->core.pos, bam_endpos(b), 14, 5);
    if (give_warning)
        hts_log_error("%s encodes a CIGAR with %d operators at the CG tag", bam_get_qname(b), c->n_cigar);
    return 1;
}

static inline int aux_type2size(uint8_t type)
{
    switch (type) {
    case 'A': case 'c': case 'C':
        return 1;
    case 's': case 'S':
        return 2;
    case 'i': case 'I': case 'f':
        return 4;
    case 'd':
        return 8;
    case 'Z': case 'H': case 'B':
        return type;
    default:
        return 0;
    }
}

static void swap_data(const bam1_core_t *c, int l_data, uint8_t *data, int is_host)
{
    uint32_t *cigar = (uint32_t*)(data + c->l_qname);
    uint32_t i;
    for (i = 0; i < c->n_cigar; ++i) ed_swap_4p(&cigar[i]);
}

// Fix bad records where qname is not terminated correctly.
static int fixup_missing_qname_nul(bam1_t *b) {
    bam1_core_t *c = &b->core;

    // Note this is called before c->l_extranul is added to c->l_qname
    if (c->l_extranul > 0) {
        b->data[c->l_qname++] = '\0';
        c->l_extranul--;
    } else {
        if (b->l_data > INT_MAX - 4) return -1;
        if (realloc_bam_data(b, b->l_data + 4) < 0) return -1;
        b->l_data += 4;
        b->data[c->l_qname++] = '\0';
        c->l_extranul = 3;
    }
    return 0;
}

/*
 * Note a second interface that returns a bam pointer instead would avoid bam_copy1
 * in multi-threaded handling.  This may be worth considering for htslib2.
 */
int bam_read1(BGZF *fp, bam1_t *b)
{
    bam1_core_t *c = &b->core;
    int32_t block_len, ret, i;
    uint32_t x[8], new_l_data;

    b->l_data = 0;

    if ((ret = bgzf_read(fp, &block_len, 4)) != 4) {
        if (ret == 0) return -1; // normal end-of-file
        else return -2; // truncated
    }
    if (fp->is_be)
        ed_swap_4p(&block_len);
    if (block_len < 32) return -4;  // block_len includes core data
    if (bgzf_read(fp, x, 32) != 32) return -3;
    if (fp->is_be) {
        for (i = 0; i < 8; ++i) ed_swap_4p(x + i);
    }
    c->tid = x[0]; c->pos = (int32_t)x[1];
    c->bin = x[2]>>16; c->qual = x[2]>>8&0xff; c->l_qname = x[2]&0xff;
    c->l_extranul = (c->l_qname%4 != 0)? (4 - c->l_qname%4) : 0;
    c->flag = x[3]>>16; c->n_cigar = x[3]&0xffff;
    c->l_qseq = x[4];
    c->mtid = x[5]; c->mpos = (int32_t)x[6]; c->isize = (int32_t)x[7];

    new_l_data = block_len - 32 + c->l_extranul;
    if (new_l_data > INT_MAX || c->l_qseq < 0 || c->l_qname < 1) return -4;
    if (((uint64_t) c->n_cigar << 2) + c->l_qname + c->l_extranul
        + (((uint64_t) c->l_qseq + 1) >> 1) + c->l_qseq > (uint64_t) new_l_data)
        return -4;
    if (realloc_bam_data(b, new_l_data) < 0) return -4;
    b->l_data = new_l_data;

    if (bgzf_read(fp, b->data, c->l_qname) != c->l_qname) return -4;
    if (b->data[c->l_qname - 1] != '\0') { // Try to fix missing NUL termination
        if (fixup_missing_qname_nul(b) < 0) return -4;
    }
    for (i = 0; i < c->l_extranul; ++i) b->data[c->l_qname+i] = '\0';
    c->l_qname += c->l_extranul;
    if (b->l_data < c->l_qname ||
        bgzf_read(fp, b->data + c->l_qname, b->l_data - c->l_qname) != b->l_data - c->l_qname)
        return -4;
    if (fp->is_be) swap_data(c, b->l_data, b->data, 0);
    if (bam_tag2cigar(b, 0, 0) < 0)
        return -4;

    if (c->n_cigar > 0) { // recompute "bin" and check CIGAR-qlen consistency
        hts_pos_t rlen, qlen;
        bam_cigar2rqlens(c->n_cigar, bam_get_cigar(b), &rlen, &qlen);
        if ((b->core.flag & BAM_FUNMAP) || rlen == 0) rlen = 1;
        b->core.bin = hts_reg2bin(b->core.pos, b->core.pos + rlen, 14, 5);
        // Sanity check for broken CIGAR alignments
        if (c->l_qseq > 0 && !(c->flag & BAM_FUNMAP) && qlen != c->l_qseq) {
            hts_log_error("CIGAR and query sequence lengths differ for %s",
                    bam_get_qname(b));
            return -4;
        }
    }

    return 4 + block_len;
}

// Initialise fp->idx for the current format type.
// This must be called after the header has been written but no other data.
int sam_idx_init(htsFile *fp, sam_hdr_t *h, int min_shift, const char *fnidx) {
    fp->fnidx = fnidx;
    if (fp->format.format == bam || fp->format.format == bcf ||
        (fp->format.format == sam && fp->format.compression == bgzf)) {
        int n_lvls, fmt = HTS_FMT_CSI;
        if (min_shift > 0) {
            int64_t max_len = 0, s;
            int i;
            for (i = 0; i < h->n_targets; ++i)
                if (max_len < h->target_len[i]) max_len = h->target_len[i];
            max_len += 256;
            for (n_lvls = 0, s = 1<<min_shift; max_len > s; ++n_lvls, s <<= 3);

        } else min_shift = 14, n_lvls = 5, fmt = HTS_FMT_BAI;

        fp->idx = hts_idx_init(h->n_targets, fmt, bgzf_tell(fp->fp.bgzf), min_shift, n_lvls);
        return fp->idx ? 0 : -1;
    }

    if (fp->format.format == cram) {
        fp->fp.cram->idxfp = bgzf_open(fnidx, "wg");
        return fp->fp.cram->idxfp ? 0 : -1;
    }

    return -1;
}

// Bam record pointer and SAM header combined
typedef struct {
    const sam_hdr_t *h;
    const bam1_t *b;
} hb_pair;


/**********************
 *** SAM header I/O ***
 **********************/

#include "htslib/kseq.h"
#include "htslib/kstring.h"

sam_hdr_t *sam_hdr_parse(size_t l_text, const char *text)
{
    sam_hdr_t *bh = sam_hdr_init();
    if (!bh) return NULL;

    if (sam_hdr_add_lines(bh, text, l_text) != 0) {
        sam_hdr_destroy(bh);
        return NULL;
    }

    return bh;
}

static int valid_sam_header_type(const char *s) {
    if (s[0] != '@') return 0;
    switch (s[1]) {
    case 'H':
        return s[2] == 'D' && s[3] == '\t';
    case 'S':
        return s[2] == 'Q' && s[3] == '\t';
    case 'R':
    case 'P':
        return s[2] == 'G' && s[3] == '\t';
    case 'C':
        return s[2] == 'O';
    }
    return 0;
}

// Minimal sanitisation of a header to ensure.
// - null terminated string.
// - all lines start with @ (also implies no blank lines).
//
// Much more could be done, but currently is not, including:
// - checking header types are known (HD, SQ, etc).
// - syntax (eg checking tab separated fields).
// - validating n_targets matches @SQ records.
// - validating target lengths against @SQ records.
static sam_hdr_t *sam_hdr_sanitise(sam_hdr_t *h) {
    if (!h)
        return NULL;

    // Special case for empty headers.
    if (h->l_text == 0)
        return h;

    size_t i;
    unsigned int lnum = 0;
    char *cp = h->text, last = '\n';
    for (i = 0; i < h->l_text; i++) {
        // NB: l_text excludes terminating nul.  This finds early ones.
        if (cp[i] == 0)
            break;

        // Error on \n[^@], including duplicate newlines
        if (last == '\n') {
            lnum++;
            if (cp[i] != '@') {
                hts_log_error("Malformed SAM header at line %u", lnum);
                sam_hdr_destroy(h);
                return NULL;
            }
        }

        last = cp[i];
    }

    if (i < h->l_text) { // Early nul found.  Complain if not just padding.
        size_t j = i;
        while (j < h->l_text && cp[j] == '\0') j++;
        if (j < h->l_text)
            hts_log_warning("Unexpected NUL character in header. Possibly truncated");
    }

    // Add trailing newline and/or trailing nul if required.
    if (last != '\n') {
        hts_log_warning("Missing trailing newline on SAM header. Possibly truncated");

        if (h->l_text < 2 || i >= h->l_text - 2) {
            if (h->l_text >= SIZE_MAX - 2) {
                hts_log_error("No room for extra newline");
                sam_hdr_destroy(h);
                return NULL;
            }

            cp = realloc(h->text, (size_t) h->l_text+2);
            if (!cp) {
                sam_hdr_destroy(h);
                return NULL;
            }
            h->text = cp;
        }
        cp[i++] = '\n';

        // l_text may be larger already due to multiple nul padding
        if (h->l_text < i)
            h->l_text = i;
        cp[h->l_text] = '\0';
    }

    return h;
}

static void known_stderr(const char *tool, const char *advice) {
    hts_log_warning("SAM file corrupted by embedded %s error/log message", tool);
    hts_log_warning("%s", advice);
}

static void warn_if_known_stderr(const char *line) {
    if (strstr(line, "M::bwa_idx_load_from_disk") != NULL)
        known_stderr("bwa", "Use `bwa mem -o file.sam ...` or `bwa sampe -f file.sam ...` instead of `bwa ... > file.sam`");
    else if (strstr(line, "M::mem_pestat") != NULL)
        known_stderr("bwa", "Use `bwa mem -o file.sam ...` instead of `bwa mem ... > file.sam`");
    else if (strstr(line, "loaded/built the index") != NULL)
        known_stderr("minimap2", "Use `minimap2 -o file.sam ...` instead of `minimap2 ... > file.sam`");
}

static sam_hdr_t *sam_hdr_create(htsFile* fp) {
    kstring_t str = { 0, 0, NULL };
    khint_t k;
    sam_hdr_t* h = sam_hdr_init();
    const char *q, *r;
    char* sn = NULL;
    khash_t(s2i) *d = kh_init(s2i);
    khash_t(s2i) *long_refs = NULL;
    if (!h || !d)
        goto error;

    int ret, has_SQ = 0;
    int next_c = '@';
    while (next_c == '@' && (ret = hts_getline(fp, KS_SEP_LINE, &fp->line)) >= 0) {
        if (fp->line.s[0] != '@')
            break;

        if (fp->line.l > 3 && strncmp(fp->line.s, "@SQ", 3) == 0) {
            has_SQ = 1;
            hts_pos_t ln = -1;
            for (q = fp->line.s + 4;; ++q) {
                if (strncmp(q, "SN:", 3) == 0) {
                    q += 3;
                    for (r = q;*r != '\t' && *r != '\n' && *r != '\0';++r);

                    if (sn) {
                        hts_log_warning("SQ header line has more than one SN: tag");
                        free(sn);
                    }
                    sn = (char*)calloc(r - q + 1, 1);
                    if (!sn)
                        goto error;

                    strncpy(sn, q, r - q);
                    q = r;
                } else {
                    if (strncmp(q, "LN:", 3) == 0)
                        ln = strtoll(q + 3, (char**)&q, 10);
                }

                while (*q != '\t' && *q != '\n' && *q != '\0')
                    ++q;
                if (*q == '\0' || *q == '\n')
                    break;
            }
            if (sn) {
                if (ln >= 0) {
                    int absent;
                    k = kh_put(s2i, d, sn, &absent);
                    if (absent < 0)
                        goto error;

                    if (!absent) {
                        hts_log_warning("Duplicated sequence \"%s\" in file \"%s\"", sn, fp->fn);
                        free(sn);
                    } else {
                        sn = NULL;
                        if (ln >= UINT32_MAX) {
                            // Stash away ref length that
                            // doesn't fit in target_len array
                            int k2;
                            if (!long_refs) {
                                long_refs = kh_init(s2i);
                                if (!long_refs)
                                    goto error;
                            }
                            k2 = kh_put(s2i, long_refs, kh_key(d, k), &absent);
                            if (absent < 0)
                                goto error;
                            kh_val(long_refs, k2) = ln;
                            kh_val(d, k) = ((int64_t) (kh_size(d) - 1) << 32
                                            | UINT32_MAX);
                        } else {
                            kh_val(d, k) = (int64_t) (kh_size(d) - 1) << 32 | ln;
                        }
                    }
                } else {
                    hts_log_warning("Ignored @SQ SN:%s : bad or missing LN tag", sn);
                    warn_if_known_stderr(fp->line.s);
                    free(sn);
                }
            } else {
                hts_log_warning("Ignored @SQ line with missing SN: tag");
                warn_if_known_stderr(fp->line.s);
            }
            sn = NULL;
        }
        else if (!valid_sam_header_type(fp->line.s)) {
            hts_log_error("Invalid header line: must start with @HD/@SQ/@RG/@PG/@CO");
            warn_if_known_stderr(fp->line.s);
            goto error;
        }

        if (kputsn(fp->line.s, fp->line.l, &str) < 0)
            goto error;

        if (kputc('\n', &str) < 0)
            goto error;

        if (fp->is_bgzf) {
            next_c = bgzf_peek(fp->fp.bgzf);
        } else {
            unsigned char nc;
            ssize_t pret = hpeek(fp->fp.hfile, &nc, 1);
            next_c = pret > 0 ? nc : pret - 1;
        }
        if (next_c < -1)
            goto error;
    }
    if (next_c != '@')
        fp->line.l = 0;

    if (ret < -1)
        goto error;

    if (!has_SQ && fp->fn_aux) {
        kstring_t line = { 0, 0, NULL };

        /* The reference index (.fai) is actually needed here */
        char *fai_fn = fp->fn_aux;
        char *fn_delim = strstr(fp->fn_aux, HTS_IDX_DELIM);
        if (fn_delim)
            fai_fn = fn_delim + strlen(HTS_IDX_DELIM);

        hFILE* f = hopen(fai_fn, "r");
        int e = 0, absent;
        if (f == NULL)
            goto error;

        while (line.l = 0, kgetline(&line, (kgets_func*) hgets, f) >= 0) {
            char* tab = strchr(line.s, '\t');
            hts_pos_t ln;

            if (tab == NULL)
                continue;

            sn = (char*)calloc(tab-line.s+1, 1);
            if (!sn) {
                e = 1;
                break;
            }
            memcpy(sn, line.s, tab-line.s);
            k = kh_put(s2i, d, sn, &absent);
            if (absent < 0) {
                e = 1;
                break;
            }

            ln = strtoll(tab, NULL, 10);

            if (!absent) {
                hts_log_warning("Duplicated sequence \"%s\" in the file \"%s\"", sn, fai_fn);
                free(sn);
                sn = NULL;
            } else {
                sn = NULL;
                if (ln >= UINT32_MAX) {
                    // Stash away ref length that
                    // doesn't fit in target_len array
                    khint_t k2;
                    int absent = -1;
                    if (!long_refs) {
                        long_refs = kh_init(s2i);
                        if (!long_refs) {
                            e = 1;
                            break;
                        }
                    }
                    k2 = kh_put(s2i, long_refs, kh_key(d, k), &absent);
                    if (absent < 0) {
                         e = 1;
                         break;
                    }
                    kh_val(long_refs, k2) = ln;
                    kh_val(d, k) = ((int64_t) (kh_size(d) - 1) << 32
                                    | UINT32_MAX);
                } else {
                    kh_val(d, k) = (int64_t) (kh_size(d) - 1) << 32 | ln;
                }
                has_SQ = 1;
            }

            e |= kputs("@SQ\tSN:", &str) < 0;
            e |= kputsn(line.s, tab - line.s, &str) < 0;
            e |= kputs("\tLN:", &str) < 0;
            e |= kputll(ln, &str) < 0;
            e |= kputc('\n', &str) < 0;
            if (e)
                break;
        }

        ks_free(&line);
        if (hclose(f) != 0) {
            hts_log_error("Error on closing %s", fai_fn);
            e = 1;
        }
        if (e)
            goto error;
    }

    if (has_SQ) {
        // Populate the targets array
        h->n_targets = kh_size(d);

        h->target_name = (char**) malloc(sizeof(char*) * h->n_targets);
        if (!h->target_name) {
            h->n_targets = 0;
            goto error;
        }

        h->target_len = (uint32_t*) malloc(sizeof(uint32_t) * h->n_targets);
        if (!h->target_len) {
            h->n_targets = 0;
            goto error;
        }

        for (k = kh_begin(d); k != kh_end(d); ++k) {
            if (!kh_exist(d, k))
                continue;

            h->target_name[kh_val(d, k) >> 32] = (char*) kh_key(d, k);
            h->target_len[kh_val(d, k) >> 32] = kh_val(d, k) & 0xffffffffUL;
            kh_val(d, k) >>= 32;
        }
    }

    // Repurpose sdict to hold any references longer than UINT32_MAX
    h->sdict = long_refs;

    kh_destroy(s2i, d);

    if (str.l == 0)
        kputsn("", 0, &str);
    h->l_text = str.l;
    h->text = ks_release(&str);
    fp->bam_header = sam_hdr_sanitise(h);
    fp->bam_header->ref_count = 1;

    return fp->bam_header;

 error:
    if (h && d && (!h->target_name || !h->target_len)) {
        for (k = kh_begin(d); k != kh_end(d); ++k)
            if (kh_exist(d, k)) free((void *)kh_key(d, k));
    }
    sam_hdr_destroy(h);
    ks_free(&str);
    kh_destroy(s2i, d);
    kh_destroy(s2i, long_refs);
    if (sn) free(sn);
    return NULL;
}

sam_hdr_t *sam_hdr_read(htsFile *fp)
{
    if (!fp) {
        errno = EINVAL;
        return NULL;
    }

    switch (fp->format.format) {
    case bam:
        return sam_hdr_sanitise(bam_hdr_read(fp->fp.bgzf));

    case cram:
        return sam_hdr_sanitise(sam_hdr_dup(fp->fp.cram->header));

    case sam:
        return sam_hdr_create(fp);

    case fastq_format:
    case fasta_format:
        return sam_hdr_init();

    case empty_format:
        errno = EPIPE;
        return NULL;

    default:
        errno = EFTYPE;
        return NULL;
    }
}

static inline unsigned int parse_sam_flag(char *v, char **rv, int *overflow) {
    if (*v >= '1' && *v <= '9') {
        return hts_str2uint(v, rv, 16, overflow);
    }
    else if (*v == '0') {
        // handle single-digit "0" directly; otherwise it's hex or octal
        if (v[1] == '\t') { *rv = v+1; return 0; }
        else {
            unsigned long val = strtoul(v, rv, 0);
            if (val > 65535) { *overflow = 1; return 65535; }
            return val;
        }
    }
    else {
        // TODO implement symbolic flag letters
        *rv = v;
        return 0;
    }
}

/*
 * -----------------------------------------------------------------------------
 * SAM threading
 */
// Size of SAM text block (reading)
#define SAM_NBYTES 240000

// Number of BAM records (writing, up to NB_mem in size)
#define SAM_NBAM 1000

struct SAM_state;

// Output job - a block of BAM records
typedef struct sp_bams {
    struct sp_bams *next;
    int serial;

    bam1_t *bams;
    int nbams, abams; // used and alloc for bams[] array
    size_t bam_mem;   // very approximate total size

    struct SAM_state *fd;
} sp_bams;

// Input job - a block of SAM text
typedef struct sp_lines {
    struct sp_lines *next;
    int serial;

    char *data;
    int data_size;
    int alloc;

    struct SAM_state *fd;
    sp_bams *bams;
} sp_lines;

enum sam_cmd {
    SAM_NONE = 0,
    SAM_CLOSE,
    SAM_CLOSE_DONE,
};

typedef struct SAM_state {
    sam_hdr_t *h;

    hts_tpool *p;
    int own_pool;
    pthread_mutex_t lines_m;
    hts_tpool_process *q;
    pthread_t dispatcher;
    int dispatcher_set;

    sp_lines *lines;
    sp_bams *bams;

    sp_bams *curr_bam;
    int curr_idx;
    int serial;

    // Be warned: moving these mutexes around in this struct can reduce
    // threading performance by up to 70%!
    pthread_mutex_t command_m;
    pthread_cond_t command_c;
    enum sam_cmd command;

    // One of the E* errno codes
    int errcode;

    htsFile *fp;
} SAM_state;

static int sam_format1_append(const bam_hdr_t *h, const bam1_t *b, kstring_t *str);

typedef struct {
    kstring_t name;
    kstring_t comment; // NB: pointer into name, do not free
    kstring_t seq;
    kstring_t qual;
    int casava;
    int aux;
    int rnum;
    char BC[3];         // aux tag ID for barcode
    khash_t(tag) *tags; // which aux tags to use (if empty, use all).
    char nprefix;
    int sra_names;
} fastq_state;

// Initialise fastq state.
// Name char of '@' or '>' distinguishes fastq vs fasta variant
static fastq_state *fastq_state_init(int name_char) {
    fastq_state *x = (fastq_state *)calloc(1, sizeof(*x));
    if (!x)
        return NULL;
    strcpy(x->BC, "BC");
    x->nprefix = name_char;

    return x;
}

void fastq_state_destroy(htsFile *fp) {
    if (fp->state) {
        fastq_state *x = (fastq_state *)fp->state;
        if (x->tags)
            kh_destroy(tag, x->tags);
        ks_free(&x->name);
        ks_free(&x->seq);
        ks_free(&x->qual);
        free(fp->state);
    }
}

int fastq_state_set(samFile *fp, enum hts_fmt_option opt, ...) {
    va_list args;

    if (!fp)
        return -1;
    if (!fp->state)
        if (!(fp->state = fastq_state_init(fp->format.format == fastq_format
                                           ? '@' : '>')))
            return -1;

    fastq_state *x = (fastq_state *)fp->state;

    switch (opt) {
    case FASTQ_OPT_CASAVA:
        x->casava = 1;
        break;

    case FASTQ_OPT_NAME2:
        x->sra_names = 1;
        break;

    case FASTQ_OPT_AUX: {
        va_start(args, opt);
        x->aux = 1;
        char *tag = va_arg(args, char *);
        va_end(args);
        if (tag && strcmp(tag, "1") != 0) {
            if (!x->tags)
                if (!(x->tags = kh_init(tag)))
                    return -1;

            size_t i, tlen = strlen(tag);
            for (i = 0; i+3 <= tlen+1; i += 3) {
                if (tag[i+0] == ',' || tag[i+1] == ',' ||
                    !(tag[i+2] == ',' || tag[i+2] == '\0')) {
                    hts_log_warning("Bad tag format '%.3s'; skipping option", tag+i);
                    break;
                }
                int ret, tcode = tag[i+0]*256 + tag[i+1];
                kh_put(tag, x->tags, tcode, &ret);
                if (ret < 0)
                    return -1;
            }
        }
        break;
    }

    case FASTQ_OPT_BARCODE: {
        va_start(args, opt);
        char *bc = va_arg(args, char *);
        va_end(args);
        strncpy(x->BC, bc, 2);
        x->BC[2] = 0;
        break;
    }

    case FASTQ_OPT_RNUM:
        x->rnum = 1;
        break;

    default:
        break;
    }
    return 0;
}

static int sam_format1_append(const bam_hdr_t *h, const bam1_t *b, kstring_t *str)
{
    int i, r = 0;
    uint8_t *s, *end;
    const bam1_core_t *c = &b->core;

    if (c->l_qname == 0)
        return -1;
    r |= kputsn_(bam_get_qname(b), c->l_qname-1-c->l_extranul, str);
    r |= kputc_('\t', str); // query name
    r |= kputw(c->flag, str); r |= kputc_('\t', str); // flag
    if (c->tid >= 0) { // chr
        r |= kputs(h->target_name[c->tid] , str);
        r |= kputc_('\t', str);
    } else r |= kputsn_("*\t", 2, str);
    r |= kputll(c->pos + 1, str); r |= kputc_('\t', str); // pos
    r |= kputw(c->qual, str); r |= kputc_('\t', str); // qual
    if (c->n_cigar) { // cigar
        uint32_t *cigar = bam_get_cigar(b);
        for (i = 0; i < c->n_cigar; ++i) {
            r |= kputw(bam_cigar_oplen(cigar[i]), str);
            r |= kputc_(bam_cigar_opchr(cigar[i]), str);
        }
    } else r |= kputc_('*', str);
    r |= kputc_('\t', str);
    if (c->mtid < 0) r |= kputsn_("*\t", 2, str); // mate chr
    else if (c->mtid == c->tid) r |= kputsn_("=\t", 2, str);
    else {
        r |= kputs(h->target_name[c->mtid], str);
        r |= kputc_('\t', str);
    }
    r |= kputll(c->mpos + 1, str); r |= kputc_('\t', str); // mate pos
    r |= kputll(c->isize, str); r |= kputc_('\t', str); // template len
    if (c->l_qseq) { // seq and qual
        uint8_t *s = bam_get_seq(b);
        if (ks_resize(str, str->l+2+2*c->l_qseq) < 0) goto mem_err;
        char *cp = str->s + str->l;

        // Sequence, 2 bases at a time
        nibble2base(s, cp, c->l_qseq);
        cp[c->l_qseq] = '\t';
        cp += c->l_qseq+1;

        // Quality
        s = bam_get_qual(b);
        i = 0;
        if (s[0] == 0xff) {
            cp[i++] = '*';
        } else {
            // local copy of c->l_qseq to aid unrolling
            uint32_t lqseq = c->l_qseq;
            for (i = 0; i < lqseq; ++i)
                cp[i]=s[i]+33;
        }
        cp[i] = 0;
        cp += i;
        str->l = cp - str->s;
    } else r |= kputsn_("*\t*", 3, str);

    s = bam_get_aux(b); // aux
    end = b->data + b->l_data;

    while (end - s >= 4) {
        r |= kputc_('\t', str);
        if ((s = (uint8_t *)sam_format_aux1(s, s[2], s+3, end, str)) == NULL)
            goto bad_aux;
    }
    r |= kputsn("", 0, str); // nul terminate
    if (r < 0) goto mem_err;

    return str->l;

 bad_aux:
    hts_log_error("Corrupted aux data for read %.*s",
                  b->core.l_qname, bam_get_qname(b));
    errno = EINVAL;
    return -1;

 mem_err:
    hts_log_error("Out of memory");
    errno = ENOMEM;
    return -1;
}

int sam_format1(const bam_hdr_t *h, const bam1_t *b, kstring_t *str)
{
    str->l = 0;
    return sam_format1_append(h, b, str);
}

static inline uint8_t *skip_aux(uint8_t *s, uint8_t *end);
int fastq_format1(fastq_state *x, const bam1_t *b, kstring_t *str)
{
    unsigned flag = b->core.flag;
    int i, e = 0, len = b->core.l_qseq;
    uint8_t *seq, *qual;

    str->l = 0;

    // Name
    if (kputc(x->nprefix, str) == EOF || kputs(bam_get_qname(b), str) == EOF)
        return -1;

    // /1 or /2 suffix
    if (x && x->rnum && (flag & BAM_FPAIRED)) {
        int r12 = flag & (BAM_FREAD1 | BAM_FREAD2);
        if (r12 == BAM_FREAD1) {
            if (kputs("/1", str) == EOF)
                return -1;
        } else if (r12 == BAM_FREAD2) {
            if (kputs("/2", str) == EOF)
                return -1;
        }
    }

    // Illumina CASAVA tag.
    // This is <rnum>:<Y/N qcfail>:<control-bits>:<barcode-or-zero>
    if (x && x->casava) {
        int rnum = (flag & BAM_FREAD1)? 1 : (flag & BAM_FREAD2)? 2 : 0;
        char filtered = (flag & BAM_FQCFAIL)? 'Y' : 'N';
        uint8_t *bc = bam_aux_get(b, x->BC);
        if (ksprintf(str, " %d:%c:0:%s", rnum, filtered,
                     bc ? (char *)bc+1 : "0") < 0)
            return -1;

        if (bc && (*bc != 'Z' || (!isupper_c(bc[1]) && !islower_c(bc[1])))) {
            hts_log_warning("BC tag starts with non-sequence base; using '0'");
            str->l -= strlen((char *)bc)-2; // limit to 1 char
            str->s[str->l-1] = '0';
            str->s[str->l] = 0;
            bc = NULL;
        }

        // Replace any non-alpha with '+'.  Ie seq-seq to seq+seq
        if (bc) {
            int l = strlen((char *)bc+1);
            char *c = (char *)str->s + str->l - l;
            for (i = 0; i < l; i++) {
                if (!isalpha_c(c[i]))
                    c[i] = '+';
                else if (islower_c(c[i]))
                    c[i] = toupper_c(c[i]);
            }
        }
    }

    // Aux tags
    if (x && x->aux) {
        uint8_t *s = bam_get_aux(b), *end = b->data + b->l_data;
        while (s && end - s >= 4) {
            int tt = s[0]*256 + s[1];
            if (x->tags == NULL ||
                kh_get(tag, x->tags, tt) != kh_end(x->tags)) {
                e |= kputc_('\t', str) < 0;
                if (!(s = (uint8_t *)sam_format_aux1(s, s[2], s+3, end, str)))
                    return -1;
            } else {
                s = skip_aux(s+2, end);
            }
        }
        e |= kputsn("", 0, str) < 0; // nul terminate
    }

    if (ks_resize(str, str->l + 1 + len+1 + 2 + len+1 + 1) < 0) return -1;
    e |= kputc_('\n', str) < 0;

    // Seq line
    seq = bam_get_seq(b);
    if (flag & BAM_FREVERSE)
        for (i = len-1; i >= 0; i--)
            e |= kputc_("!TGKCYSBAWRDMHVN"[bam_seqi(seq, i)], str) < 0;
    else
        for (i = 0; i < len; i++)
            e |= kputc_(seq_nt16_str[bam_seqi(seq, i)], str) < 0;


    // Qual line
    if (x->nprefix == '@') {
        kputsn("\n+\n", 3, str);
        qual = bam_get_qual(b);
        if (qual[0] == 0xff)
            for (i = 0; i < len; i++)
                e |= kputc_('B', str) < 0;
        else if (flag & BAM_FREVERSE)
            for (i = len-1; i >= 0; i--)
                e |= kputc_(33 + qual[i], str) < 0;
        else
            for (i = 0; i < len; i++)
                e |= kputc_(33 + qual[i], str) < 0;

    }
    e |= kputc('\n', str) < 0;

    return e ? -1 : str->l;
}

/************************
 *** Auxiliary fields ***
 ************************/
#ifndef HTS_LITTLE_ENDIAN
static int aux_to_le(char type, uint8_t *out, const uint8_t *in, size_t len) {
    int tsz = aux_type2size(type);

    if (tsz >= 2 && tsz <= 8 && (len & (tsz - 1)) != 0) return -1;

    switch (tsz) {
        case 'H': case 'Z': case 1:  // Trivial
            memcpy(out, in, len);
            break;

#define aux_val_to_le(type_t, store_le) do {                            \
        type_t v;                                                       \
        size_t i;                                                       \
        for (i = 0; i < len; i += sizeof(type_t), out += sizeof(type_t)) { \
            memcpy(&v, in + i, sizeof(type_t));                         \
            store_le(v, out);                                           \
        }                                                               \
    } while (0)

        case 2: aux_val_to_le(uint16_t, u16_to_le); break;
        case 4: aux_val_to_le(uint32_t, u32_to_le); break;
        case 8: aux_val_to_le(uint64_t, u64_to_le); break;

#undef aux_val_to_le

        case 'B': { // Recurse!
            uint32_t n;
            if (len < 5) return -1;
            memcpy(&n, in + 1, 4);
            out[0] = in[0];
            u32_to_le(n, out + 1);
            return aux_to_le(in[0], out + 5, in + 5, len - 5);
        }

        default: // Unknown type code
            return -1;
    }



    return 0;
}
#endif

int bam_aux_append(bam1_t *b, const char tag[2], char type, int len, const uint8_t *data)
{
    uint32_t new_len;

    assert(b->l_data >= 0);
    new_len = b->l_data + 3 + len;
    if (new_len > INT32_MAX || new_len < b->l_data) goto nomem;

    if (realloc_bam_data(b, new_len) < 0) return -1;

    b->data[b->l_data] = tag[0];
    b->data[b->l_data + 1] = tag[1];
    b->data[b->l_data + 2] = type;

#ifdef HTS_LITTLE_ENDIAN
    memcpy(b->data + b->l_data + 3, data, len);
#else
    if (aux_to_le(type, b->data + b->l_data + 3, data, len) != 0) {
        errno = EINVAL;
        return -1;
    }
#endif

    b->l_data = new_len;

    return 0;

 nomem:
    errno = ENOMEM;
    return -1;
}

static inline uint8_t *skip_aux(uint8_t *s, uint8_t *end)
{
    int size;
    uint32_t n;
    if (s >= end) return end;
    size = aux_type2size(*s); ++s; // skip type
    switch (size) {
    case 'Z':
    case 'H':
        while (s < end && *s) ++s;
        return s < end ? s + 1 : end;
    case 'B':
        if (end - s < 5) return NULL;
        size = aux_type2size(*s); ++s;
        n = le_to_u32(s);
        s += 4;
        if (size == 0 || end - s < size * n) return NULL;
        return s + size * n;
    case 0:
        return NULL;
    default:
        if (end - s < size) return NULL;
        return s + size;
    }
}

uint8_t *bam_aux_first(const bam1_t *b)
{
    uint8_t *s = bam_get_aux(b);
    uint8_t *end = b->data + b->l_data;
    if (s >= end) { errno = ENOENT; return NULL; }
    return s+2;
}

uint8_t *bam_aux_next(const bam1_t *b, const uint8_t *s)
{
    uint8_t *end = b->data + b->l_data;
    uint8_t *next = s? skip_aux((uint8_t *) s, end) : end;
    if (next == NULL) goto bad_aux;
    if (next >= end) { errno = ENOENT; return NULL; }
    return next+2;

 bad_aux:
    hts_log_error("Corrupted aux data for read %s", bam_get_qname(b));
    errno = EINVAL;
    return NULL;
}

uint8_t *bam_aux_get(const bam1_t *b, const char tag[2])
{
    uint8_t *s;
    for (s = bam_aux_first(b); s; s = bam_aux_next(b, s))
        if (s[-2] == tag[0] && s[-1] == tag[1]) {
            // Check the tag value is valid and complete
            uint8_t *e = skip_aux(s, b->data + b->l_data);
            if (e == NULL) goto bad_aux;
            if ((*s == 'Z' || *s == 'H') && *(e - 1) != '\0') goto bad_aux;

            return s;
        }

    // errno now as set by bam_aux_first()/bam_aux_next()
    return NULL;

 bad_aux:
    hts_log_error("Corrupted aux data for read %s", bam_get_qname(b));
    errno = EINVAL;
    return NULL;
}

static inline int64_t get_int_aux_val(uint8_t type, const uint8_t *s,
                                      uint32_t idx)
{
    switch (type) {
        case 'c': return le_to_i8(s + idx);
        case 'C': return s[idx];
        case 's': return le_to_i16(s + 2 * idx);
        case 'S': return le_to_u16(s + 2 * idx);
        case 'i': return le_to_i32(s + 4 * idx);
        case 'I': return le_to_u32(s + 4 * idx);
        default:
            errno = EINVAL;
            return 0;
    }
}

int64_t bam_aux2i(const uint8_t *s)
{
    int type;
    type = *s++;
    return get_int_aux_val(type, s, 0);
}

double bam_aux2f(const uint8_t *s)
{
    int type;
    type = *s++;
    if (type == 'd') return le_to_double(s);
    else if (type == 'f') return le_to_float(s);
    else return get_int_aux_val(type, s, 0);
}

#define STRNCMP(a,b,n) (strncasecmp((a),(b),(n)) || strlen(a)!=(n))
int bam_str2flag(const char *str)
{
    char *end, *beg = (char*) str;
    long int flag = strtol(str, &end, 0);
    if ( end!=str ) return flag;    // the conversion was successful
    flag = 0;
    while ( *str )
    {
        end = beg;
        while ( *end && *end!=',' ) end++;
        if ( !STRNCMP("PAIRED",beg,end-beg) ) flag |= BAM_FPAIRED;
        else if ( !STRNCMP("PROPER_PAIR",beg,end-beg) ) flag |= BAM_FPROPER_PAIR;
        else if ( !STRNCMP("UNMAP",beg,end-beg) ) flag |= BAM_FUNMAP;
        else if ( !STRNCMP("MUNMAP",beg,end-beg) ) flag |= BAM_FMUNMAP;
        else if ( !STRNCMP("REVERSE",beg,end-beg) ) flag |= BAM_FREVERSE;
        else if ( !STRNCMP("MREVERSE",beg,end-beg) ) flag |= BAM_FMREVERSE;
        else if ( !STRNCMP("READ1",beg,end-beg) ) flag |= BAM_FREAD1;
        else if ( !STRNCMP("READ2",beg,end-beg) ) flag |= BAM_FREAD2;
        else if ( !STRNCMP("SECONDARY",beg,end-beg) ) flag |= BAM_FSECONDARY;
        else if ( !STRNCMP("QCFAIL",beg,end-beg) ) flag |= BAM_FQCFAIL;
        else if ( !STRNCMP("DUP",beg,end-beg) ) flag |= BAM_FDUP;
        else if ( !STRNCMP("SUPPLEMENTARY",beg,end-beg) ) flag |= BAM_FSUPPLEMENTARY;
        else return -1;
        if ( !*end ) break;
        beg = end + 1;
    }
    return flag;
}

/**************************
 *** Pileup and Mpileup ***
 **************************/

#if !defined(BAM_NO_PILEUP)

#include <assert.h>

/*******************
 *** Memory pool ***
 *******************/

typedef struct {
    int k, y;
    hts_pos_t x, end;
} cstate_t;

typedef struct __linkbuf_t {
    bam1_t b;
    hts_pos_t beg, end;
    cstate_t s;
    struct __linkbuf_t *next;
    bam_pileup_cd cd;
} lbnode_t;

typedef struct {
    int cnt, n, max;
    lbnode_t **buf;
} mempool_t;

static inline lbnode_t *mp_alloc(mempool_t *mp)
{
    ++mp->cnt;
    if (mp->n == 0) return (lbnode_t*)calloc(1, sizeof(lbnode_t));
    else return mp->buf[--mp->n];
}
static inline void mp_free(mempool_t *mp, lbnode_t *p)
{
    --mp->cnt; p->next = 0; // clear lbnode_t::next here
    if (mp->n == mp->max) {
        mp->max = mp->max? mp->max<<1 : 256;
        mp->buf = (lbnode_t**)realloc(mp->buf, sizeof(lbnode_t*) * mp->max);
    }
    mp->buf[mp->n++] = p;
}

/**********************
 *** CIGAR resolver ***
 **********************/

/* s->k: the index of the CIGAR operator that has just been processed.
   s->x: the reference coordinate of the start of s->k
   s->y: the query coordinate of the start of s->k
 */
static inline int resolve_cigar2(bam_pileup1_t *p, hts_pos_t pos, cstate_t *s)
{
#define _cop(c) ((c)&BAM_CIGAR_MASK)
#define _cln(c) ((c)>>BAM_CIGAR_SHIFT)

    bam1_t *b = p->b;
    bam1_core_t *c = &b->core;
    uint32_t *cigar = bam_get_cigar(b);
    int k;
    // determine the current CIGAR operation
    //fprintf(stderr, "%s\tpos=%ld\tend=%ld\t(%d,%ld,%d)\n", bam_get_qname(b), pos, s->end, s->k, s->x, s->y);
    if (s->k == -1) { // never processed
        p->qpos = 0;
        if (c->n_cigar == 1) { // just one operation, save a loop
          if (_cop(cigar[0]) == BAM_CMATCH || _cop(cigar[0]) == BAM_CEQUAL || _cop(cigar[0]) == BAM_CDIFF) s->k = 0, s->x = c->pos, s->y = 0;
        } else { // find the first match or deletion
            for (k = 0, s->x = c->pos, s->y = 0; k < c->n_cigar; ++k) {
                int op = _cop(cigar[k]);
                int l = _cln(cigar[k]);
                if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP ||
                    op == BAM_CEQUAL || op == BAM_CDIFF) break;
                else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) s->y += l;
            }
            assert(k < c->n_cigar);
            s->k = k;
        }
    } else { // the read has been processed before
        int op, l = _cln(cigar[s->k]);
        if (pos - s->x >= l) { // jump to the next operation
            assert(s->k < c->n_cigar); // otherwise a bug: this function should not be called in this case
            op = _cop(cigar[s->k+1]);
            if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP || op == BAM_CEQUAL || op == BAM_CDIFF) { // jump to the next without a loop
              if (_cop(cigar[s->k]) == BAM_CMATCH|| _cop(cigar[s->k]) == BAM_CEQUAL || _cop(cigar[s->k]) == BAM_CDIFF) s->y += l;
                s->x += l;
                ++s->k;
            } else { // find the next M/D/N/=/X
              if (_cop(cigar[s->k]) == BAM_CMATCH|| _cop(cigar[s->k]) == BAM_CEQUAL || _cop(cigar[s->k]) == BAM_CDIFF) s->y += l;
                s->x += l;
                for (k = s->k + 1; k < c->n_cigar; ++k) {
                    op = _cop(cigar[k]), l = _cln(cigar[k]);
                    if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP || op == BAM_CEQUAL || op == BAM_CDIFF) break;
                    else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) s->y += l;
                }
                s->k = k;
            }
            assert(s->k < c->n_cigar); // otherwise a bug
        } // else, do nothing
    }
    { // collect pileup information
        int op, l;
        op = _cop(cigar[s->k]); l = _cln(cigar[s->k]);
        p->is_del = p->indel = p->is_refskip = 0;
        if (s->x + l - 1 == pos && s->k + 1 < c->n_cigar) { // peek the next operation
            int op2 = _cop(cigar[s->k+1]);
            int l2 = _cln(cigar[s->k+1]);
            if (op2 == BAM_CDEL && op != BAM_CDEL) {
                // At start of a new deletion, merge e.g. 1D2D to 3D.
                // Within a deletion (the 2D in 1D2D) we keep p->indel=0
                // and rely on is_del=1 as we would for 3D.
                p->indel = -(int)l2;
                for (k = s->k+2; k < c->n_cigar; ++k) {
                    op2 = _cop(cigar[k]); l2 = _cln(cigar[k]);
                    if (op2 == BAM_CDEL) p->indel -= l2;
                    else break;
                }
            } else if (op2 == BAM_CINS) {
                p->indel = l2;
                for (k = s->k+2; k < c->n_cigar; ++k) {
                    op2 = _cop(cigar[k]); l2 = _cln(cigar[k]);
                    if (op2 == BAM_CINS) p->indel += l2;
                    else if (op2 != BAM_CPAD) break;
                }
            } else if (op2 == BAM_CPAD && s->k + 2 < c->n_cigar) {
                int l3 = 0;
                for (k = s->k + 2; k < c->n_cigar; ++k) {
                    op2 = _cop(cigar[k]); l2 = _cln(cigar[k]);
                    if (op2 == BAM_CINS) l3 += l2;
                    else if (op2 == BAM_CDEL || op2 == BAM_CMATCH || op2 == BAM_CREF_SKIP || op2 == BAM_CEQUAL || op2 == BAM_CDIFF) break;
                }
                if (l3 > 0) p->indel = l3;
            }
        }
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            p->qpos = s->y + (pos - s->x);
        } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
            p->is_del = 1; p->qpos = s->y; // FIXME: distinguish D and N!!!!!
            p->is_refskip = (op == BAM_CREF_SKIP);
        } // cannot be other operations; otherwise a bug
        p->is_head = (pos == c->pos); p->is_tail = (pos == s->end);
    }
    p->cigar_ind = s->k;
    return 1;
}

/***********************
 *** Pileup iterator ***
 ***********************/

// Dictionary of overlapping reads
KHASH_MAP_INIT_STR(olap_hash, lbnode_t *)
typedef khash_t(olap_hash) olap_hash_t;

struct bam_plp_s {
    mempool_t *mp;
    lbnode_t *head, *tail;
    int32_t tid, max_tid;
    hts_pos_t pos, max_pos;
    int is_eof, max_plp, error, maxcnt;
    uint64_t id;
    bam_pileup1_t *plp;
    // for the "auto" interface only
    bam1_t *b;
    bam_plp_auto_f func;
    void *data;
    olap_hash_t *overlaps;

    // For notification of creation and destruction events
    // and associated client-owned pointer.
    int (*plp_construct)(void *data, const bam1_t *b, bam_pileup_cd *cd);
    int (*plp_destruct )(void *data, const bam1_t *b, bam_pileup_cd *cd);
};


/************************
 *** Mpileup iterator ***
 ************************/

struct bam_mplp_s {
    int n;
    int32_t min_tid, *tid;
    hts_pos_t min_pos, *pos;
    bam_plp_t *iter;
    int *n_plp;
    const bam_pileup1_t **plp;
};

#endif // ~!defined(BAM_NO_PILEUP)
