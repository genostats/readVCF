/// @file htslib/sam.h
/// High-level SAM/BAM/CRAM sequence file operations.
/*
    Copyright (C) 2008, 2009, 2013-2023 Genome Research Ltd.
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

#ifndef HTSLIB_SAM_H
#define HTSLIB_SAM_H

#include <errno.h>
#include <stdint.h>
#include <sys/types.h>
#include "hts.h"
#include "hts_endian.h"

// Ensure ssize_t exists within this header. All #includes must precede this,
// and ssize_t must be undefined again at the end of this header.
#if defined _MSC_VER && defined _INTPTR_T_DEFINED && !defined _SSIZE_T_DEFINED && !defined ssize_t
#define HTSLIB_SSIZE_T
#define ssize_t intptr_t
#endif

#ifdef __cplusplus
extern "C" {
#endif

/// Highest SAM format version supported by this library
#define SAM_FORMAT_VERSION "1.6"

/***************************
 *** SAM/BAM/CRAM header ***
 ***************************/

/*! @typedef
 * @abstract Header extension structure, grouping a collection
 *  of hash tables that contain the parsed header data.
 */

typedef struct sam_hrecs_t sam_hrecs_t;

/*! @typedef
 @abstract Structure for the alignment header.
 @field n_targets   number of reference sequences
 @field l_text      length of the plain text in the header (may be zero if
                    the header has been edited)
 @field target_len  lengths of the reference sequences
 @field target_name names of the reference sequences
 @field text        plain text (may be NULL if the header has been edited)
 @field sdict       header dictionary
 @field hrecs       pointer to the extended header struct (internal use only)
 @field ref_count   reference count

 @note The text and l_text fields are included for backwards compatibility.
 These fields may be set to NULL and zero respectively as a side-effect
 of calling some header API functions.  New code that needs to access the
 header text should use the sam_hdr_str() and sam_hdr_length() functions
 instead of these fields.
 */

typedef struct sam_hdr_t {
    int32_t n_targets, ignore_sam_err;
    size_t l_text;
    uint32_t *target_len;
    const int8_t *cigar_tab HTS_DEPRECATED("Use bam_cigar_table[] instead");
    char **target_name;
    char *text;
    void *sdict;
    sam_hrecs_t *hrecs;
    uint32_t ref_count;
} sam_hdr_t;

/*! @typedef
 * @abstract Old name for compatibility with existing code.
 */
typedef sam_hdr_t bam_hdr_t;

/****************************
 *** CIGAR related macros ***
 ****************************/

#define BAM_CMATCH      0
#define BAM_CINS        1
#define BAM_CDEL        2
#define BAM_CREF_SKIP   3
#define BAM_CSOFT_CLIP  4
#define BAM_CHARD_CLIP  5
#define BAM_CPAD        6
#define BAM_CEQUAL      7
#define BAM_CDIFF       8
#define BAM_CBACK       9

#define BAM_CIGAR_STR   "MIDNSHP=XB"
#define BAM_CIGAR_SHIFT 4
#define BAM_CIGAR_MASK  0xf
#define BAM_CIGAR_TYPE  0x3C1A7

/*! @abstract Table for converting a CIGAR operator character to BAM_CMATCH etc.
Result is operator code or -1. Be sure to cast the index if it is a plain char:
    int op = bam_cigar_table[(unsigned char) ch];
*/
HTSLIB_EXPORT
extern const int8_t bam_cigar_table[256];

#define bam_cigar_op(c) ((c)&BAM_CIGAR_MASK)
#define bam_cigar_oplen(c) ((c)>>BAM_CIGAR_SHIFT)
// Note that BAM_CIGAR_STR is padded to length 16 bytes below so that
// the array look-up will not fall off the end.  '?' is chosen as the
// padding character so it's easy to spot if one is emitted, and will
// result in a parsing failure (in sam_parse1(), at least) if read.
#define bam_cigar_opchr(c) (BAM_CIGAR_STR "??????" [bam_cigar_op(c)])
#define bam_cigar_gen(l, o) ((l)<<BAM_CIGAR_SHIFT|(o))

/* bam_cigar_type returns a bit flag with:
 *   bit 1 set if the cigar operation consumes the query
 *   bit 2 set if the cigar operation consumes the reference
 *
 * For reference, the unobfuscated truth table for this function is:
 * BAM_CIGAR_TYPE  QUERY  REFERENCE
 * --------------------------------
 * BAM_CMATCH      1      1
 * BAM_CINS        1      0
 * BAM_CDEL        0      1
 * BAM_CREF_SKIP   0      1
 * BAM_CSOFT_CLIP  1      0
 * BAM_CHARD_CLIP  0      0
 * BAM_CPAD        0      0
 * BAM_CEQUAL      1      1
 * BAM_CDIFF       1      1
 * BAM_CBACK       0      0
 * --------------------------------
 */
#define bam_cigar_type(o) (BAM_CIGAR_TYPE>>((o)<<1)&3) // bit 1: consume query; bit 2: consume reference

/*! @abstract the read is paired in sequencing, no matter whether it is mapped in a pair */
#define BAM_FPAIRED        1
/*! @abstract the read is mapped in a proper pair */
#define BAM_FPROPER_PAIR   2
/*! @abstract the read itself is unmapped; conflictive with BAM_FPROPER_PAIR */
#define BAM_FUNMAP         4
/*! @abstract the mate is unmapped */
#define BAM_FMUNMAP        8
/*! @abstract the read is mapped to the reverse strand */
#define BAM_FREVERSE      16
/*! @abstract the mate is mapped to the reverse strand */
#define BAM_FMREVERSE     32
/*! @abstract this is read1 */
#define BAM_FREAD1        64
/*! @abstract this is read2 */
#define BAM_FREAD2       128
/*! @abstract not primary alignment */
#define BAM_FSECONDARY   256
/*! @abstract QC failure */
#define BAM_FQCFAIL      512
/*! @abstract optical or PCR duplicate */
#define BAM_FDUP        1024
/*! @abstract supplementary alignment */
#define BAM_FSUPPLEMENTARY 2048

/*************************
 *** Alignment records ***
 *************************/

/*
 * Assumptions made here.  While pos can be 64-bit, no sequence
 * itself is that long, but due to ref skip CIGAR fields it
 * may span more than that.  (CIGAR itself is 28-bit len + 4 bit
 * type, but in theory we can combine multiples together.)
 *
 * Mate position and insert size also need to be 64-bit, but
 * we won't accept more than 32-bit for tid.
 *
 * The bam1_core_t structure is the *in memory* layout and not
 * the same as the on-disk format.  64-bit changes here permit
 * SAM to work with very long chromosomes and permit BAM and CRAM
 * to seamlessly update in the future without further API/ABI
 * revisions.
 */

/*! @typedef
 @abstract Structure for core alignment information.
 @field  pos     0-based leftmost coordinate
 @field  tid     chromosome ID, defined by sam_hdr_t
 @field  bin     bin calculated by bam_reg2bin()
 @field  qual    mapping quality
 @field  l_extranul length of extra NULs between qname & cigar (for alignment)
 @field  flag    bitwise flag
 @field  l_qname length of the query name
 @field  n_cigar number of CIGAR operations
 @field  l_qseq  length of the query sequence (read)
 @field  mtid    chromosome ID of next read in template, defined by sam_hdr_t
 @field  mpos    0-based leftmost coordinate of next read in template
 @field  isize   observed template length ("insert size")
 */
typedef struct bam1_core_t {
    hts_pos_t pos;
    int32_t tid;
    uint16_t bin; // NB: invalid on 64-bit pos
    uint8_t qual;
    uint8_t l_extranul;
    uint16_t flag;
    uint16_t l_qname;
    uint32_t n_cigar;
    int32_t l_qseq;
    int32_t mtid;
    hts_pos_t mpos;
    hts_pos_t isize;
} bam1_core_t;

/*! @typedef
 @abstract Structure for one alignment.
 @field  core       core information about the alignment
 @field  id
 @field  data       all variable-length data, concatenated; structure: qname-cigar-seq-qual-aux
 @field  l_data     current length of bam1_t::data
 @field  m_data     maximum length of bam1_t::data
 @field  mempolicy  memory handling policy, see bam_set_mempolicy()

 @discussion Notes:

 1. The data blob should be accessed using bam_get_qname, bam_get_cigar,
    bam_get_seq, bam_get_qual and bam_get_aux macros.  These returns pointers
    to the start of each type of data.
 2. qname is terminated by one to four NULs, so that the following
    cigar data is 32-bit aligned; core.l_qname includes these trailing NULs,
    while core.l_extranul counts the excess NULs (so 0 <= l_extranul <= 3).
 3. Cigar data is encoded 4 bytes per CIGAR operation.
    See the bam_cigar_* macros for manipulation.
 4. seq is nibble-encoded according to bam_nt16_table.
    See the bam_seqi macro for retrieving individual bases.
 5. Per base qualities are stored in the Phred scale with no +33 offset.
    Ie as per the BAM specification and not the SAM ASCII printable method.
 */
typedef struct bam1_t {
    bam1_core_t core;
    uint64_t id;
    uint8_t *data;
    int l_data;
    uint32_t m_data;
    uint32_t mempolicy:2, :30 /* Reserved */;
} bam1_t;

/*! @function
 @abstract  Get whether the query is on the reverse strand
 @param  b  pointer to an alignment
 @return    boolean true if query is on the reverse strand
 */
#define bam_is_rev(b) (((b)->core.flag&BAM_FREVERSE) != 0)
/*! @function
 @abstract  Get whether the query's mate is on the reverse strand
 @param  b  pointer to an alignment
 @return    boolean true if query's mate on the reverse strand
 */
#define bam_is_mrev(b) (((b)->core.flag&BAM_FMREVERSE) != 0)
/*! @function
 @abstract  Get the name of the query
 @param  b  pointer to an alignment
 @return    pointer to the name string, null terminated
 */
#define bam_get_qname(b) ((char*)(b)->data)
/*! @function
 @abstract  Get the CIGAR array
 @param  b  pointer to an alignment
 @return    pointer to the CIGAR array

 @discussion In the CIGAR array, each element is a 32-bit integer. The
 lower 4 bits gives a CIGAR operation and the higher 28 bits keep the
 length of a CIGAR.
 */
#define bam_get_cigar(b) ((uint32_t*)((b)->data + (b)->core.l_qname))
/*! @function
 @abstract  Get query sequence
 @param  b  pointer to an alignment
 @return    pointer to sequence

 @discussion Each base is encoded in 4 bits: 1 for A, 2 for C, 4 for G,
 8 for T and 15 for N. Two bases are packed in one byte with the base
 at the higher 4 bits having smaller coordinate on the read. It is
 recommended to use bam_seqi() macro to get the base.
 */
#define bam_get_seq(b)   ((b)->data + ((b)->core.n_cigar<<2) + (b)->core.l_qname)
/*! @function
 @abstract  Get query quality
 @param  b  pointer to an alignment
 @return    pointer to quality string
 */
#define bam_get_qual(b)  ((b)->data + ((b)->core.n_cigar<<2) + (b)->core.l_qname + (((b)->core.l_qseq + 1)>>1))
/*! @function
 @abstract  Get auxiliary data
 @param  b  pointer to an alignment
 @return    pointer to the concatenated auxiliary data
 */
#define bam_get_aux(b)   ((b)->data + ((b)->core.n_cigar<<2) + (b)->core.l_qname + (((b)->core.l_qseq + 1)>>1) + (b)->core.l_qseq)
/*! @function
 @abstract  Get length of auxiliary data
 @param  b  pointer to an alignment
 @return    length of the concatenated auxiliary data
 */
#define bam_get_l_aux(b) ((b)->l_data - ((b)->core.n_cigar<<2) - (b)->core.l_qname - (b)->core.l_qseq - (((b)->core.l_qseq + 1)>>1))
/*! @function
 @abstract  Get a base on read
 @param  s  Query sequence returned by bam_get_seq()
 @param  i  The i-th position, 0-based
 @return    4-bit integer representing the base.
 */
#define bam_seqi(s, i) ((s)[(i)>>1] >> ((~(i)&1)<<2) & 0xf)
/*!
 @abstract  Modifies a single base in the bam structure.
 @param s   Query sequence returned by bam_get_seq()
 @param i   The i-th position, 0-based
 @param b   Base in nt16 nomenclature (see seq_nt16_table)
*/
#define bam_set_seqi(s,i,b) ((s)[(i)>>1] = ((s)[(i)>>1] & (0xf0 >> ((~(i)&1)<<2))) | ((b)<<((~(i)&1)<<2)))

/**************************
 *** Exported functions ***
 **************************/

/***************
 *** BAM I/O ***
 ***************/

/* Header */

/// Generates a new unpopulated header structure.
/*!
 *
 * @return  A valid pointer to new header on success, NULL on failure
 *
 * The sam_hdr_t struct returned by a successful call should be freed
 * via sam_hdr_destroy() when it is no longer needed.
 */
HTSLIB_EXPORT
sam_hdr_t *sam_hdr_init(void);

/// Read the header from a BAM compressed file.
/*!
 * @param fp  File pointer
 * @return    A valid pointer to new header on success, NULL on failure
 *
 * This function only works with BAM files.  It is usually better to use
 * sam_hdr_read(), which works on SAM, BAM and CRAM files.
 *
 * The sam_hdr_t struct returned by a successful call should be freed
 * via sam_hdr_destroy() when it is no longer needed.
 */
HTSLIB_EXPORT
sam_hdr_t *bam_hdr_read(BGZF *fp);

/// Writes the header to a BAM file.
/*!
 * @param fp  File pointer
 * @param h   Header pointer
 * @return    0 on success, -1 on failure
 *
 * This function only works with BAM files.  Use sam_hdr_write() to
 * write in any of the SAM, BAM or CRAM formats.
 */
HTSLIB_EXPORT
int bam_hdr_write(BGZF *fp, const sam_hdr_t *h) HTS_RESULT_USED;

/*!
 * Frees the resources associated with a header.
 */

/// Duplicate a header structure.
/*!
 * @return  A valid pointer to new header on success, NULL on failure
 *
 * The sam_hdr_t struct returned by a successful call should be freed
 * via sam_hdr_destroy() when it is no longer needed.
 */
HTSLIB_EXPORT
sam_hdr_t *sam_hdr_dup(const sam_hdr_t *h0);

/*!
 * @abstract Old names for compatibility with existing code.
 */
static inline sam_hdr_t *bam_hdr_init(void) { return sam_hdr_init(); }
static inline sam_hdr_t *bam_hdr_dup(const sam_hdr_t *h0) { return sam_hdr_dup(h0); }

typedef htsFile samFile;

/// Create a header from existing text.
/*!
 * @param l_text    Length of text
 * @param text      Header text
 * @return A populated sam_hdr_t structure on success; NULL on failure.
 * @note The text field of the returned header will be NULL, and the l_text
 * field will be zero.
 *
 * The sam_hdr_t struct returned by a successful call should be freed
 * via sam_hdr_destroy() when it is no longer needed.
 */
HTSLIB_EXPORT
sam_hdr_t *sam_hdr_parse(size_t l_text, const char *text);

/// Read a header from a SAM, BAM or CRAM file.
/*!
 * @param fp    Pointer to a SAM, BAM or CRAM file handle
 * @return  A populated sam_hdr_t struct on success; NULL on failure.
 *
 * The sam_hdr_t struct returned by a successful call should be freed
 * via sam_hdr_destroy() when it is no longer needed.
 */
HTSLIB_EXPORT
sam_hdr_t *sam_hdr_read(samFile *fp);

/// Returns the current length of the header text.
/*!
 * @return  >= 0 on success, SIZE_MAX on failure
 */
HTSLIB_EXPORT
size_t sam_hdr_length(sam_hdr_t *h);

/// Returns the text representation of the header.
/*!
 * @return  valid char pointer on success, NULL on failure
 *
 * The returned string is part of the header structure.  It will remain
 * valid until a call to a header API function causes the string to be
 * invalidated, or the header is destroyed.
 *
 * The caller should not attempt to free or realloc this pointer.
 */
HTSLIB_EXPORT
const char *sam_hdr_str(sam_hdr_t *h);

/* ==== Line level methods ==== */

/// Add formatted lines to an existing header.
/*!
 * @param lines  Full SAM header record, eg "@SQ\tSN:foo\tLN:100", with
 *               optional new-line. If it contains more than 1 line then
 *               multiple lines will be added in order
 * @param len    The maximum length of lines (if an early NUL is not
 *               encountered). len may be 0 if unknown, in which case
 *               lines must be NUL-terminated
 * @return       0 on success, -1 on failure
 *
 * The lines will be appended to the end of the existing header
 * (apart from HD, which always comes first).
 */
HTSLIB_EXPORT
int sam_hdr_add_lines(sam_hdr_t *h, const char *lines, size_t len);

/// Adds a single line to an existing header.
/*!
 * Specify type and one or more key,value pairs, ending with the NULL key.
 * Eg. sam_hdr_add_line(h, "SQ", "SN", "foo", "LN", "100", NULL).
 *
 * @param type  Type of the added line. Eg. "SQ"
 * @return      0 on success, -1 on failure
 *
 * The new line will be added immediately after any others of the same
 * type, or at the end of the existing header if no lines of the
 * given type currently exist.  The exception is HD lines, which always
 * come first.  If an HD line already exists, it will be replaced.
 */
HTSLIB_EXPORT
int sam_hdr_add_line(sam_hdr_t *h, const char *type, ...);

/// Returns a complete line of formatted text for a given type and ID.
/*!
 * @param type      Type of the searched line. Eg. "SQ"
 * @param ID_key    Tag key defining the line. Eg. "SN"
 * @param ID_value  Tag value associated with the key above. Eg. "ref1"
 * @param ks        kstring to hold the result
 * @return          0 on success;
 *                 -1 if no matching line is found
 *                 -2 on other failures
 *
 * Puts a complete line of formatted text for a specific header type/ID
 * combination into @p ks. If ID_key is NULL then it returns the first line of
 * the specified type.
 *
 * Any existing content in @p ks will be overwritten.
 */
HTSLIB_EXPORT
int sam_hdr_find_line_id(sam_hdr_t *h, const char *type,
                      const char *ID_key, const char *ID_val, kstring_t *ks);

/// Returns a complete line of formatted text for a given type and index.
/*!
 * @param type      Type of the searched line. Eg. "SQ"
 * @param position  Index in lines of this type (zero-based)
 * @param ks        kstring to hold the result
 * @return          0 on success;
 *                 -1 if no matching line is found
 *                 -2 on other failures
 *
 * Puts a complete line of formatted text for a specific line into @p ks.
 * The header line is selected using the @p type and @p position parameters.
 *
 * Any existing content in @p ks will be overwritten.
 */
HTSLIB_EXPORT
int sam_hdr_find_line_pos(sam_hdr_t *h, const char *type,
                          int pos, kstring_t *ks);

/// Remove a line with given type / id from a header
/*!
 * @param type      Type of the searched line. Eg. "SQ"
 * @param ID_key    Tag key defining the line. Eg. "SN"
 * @param ID_value  Tag value associated with the key above. Eg. "ref1"
 * @return          0 on success, -1 on error
 *
 * Remove a line from the header by specifying a tag:value that uniquely
 * identifies the line, i.e. the @SQ line containing "SN:ref1".
 *
 * \@SQ line is uniquely identified by the SN tag.
 * \@RG line is uniquely identified by the ID tag.
 * \@PG line is uniquely identified by the ID tag.
 * Eg. sam_hdr_remove_line_id(h, "SQ", "SN", "ref1")
 *
 * If no key:value pair is specified, the type MUST be followed by a NULL argument and
 * the first line of the type will be removed, if any.
 * Eg. sam_hdr_remove_line_id(h, "SQ", NULL, NULL)
 *
 * @note Removing \@PG lines is currently unsupported.
 */
HTSLIB_EXPORT
int sam_hdr_remove_line_id(sam_hdr_t *h, const char *type, const char *ID_key, const char *ID_value);

/// Remove nth line of a given type from a header
/*!
 * @param type     Type of the searched line. Eg. "SQ"
 * @param position Index in lines of this type (zero-based). E.g. 3
 * @return         0 on success, -1 on error
 *
 * Remove a line from the header by specifying the position in the type
 * group, i.e. 3rd @SQ line.
 */
HTSLIB_EXPORT
int sam_hdr_remove_line_pos(sam_hdr_t *h, const char *type, int position);

/// Add or update tag key,value pairs in a header line.
/*!
 * @param type      Type of the searched line. Eg. "SQ"
 * @param ID_key    Tag key defining the line. Eg. "SN"
 * @param ID_value  Tag value associated with the key above. Eg. "ref1"
 * @return          0 on success, -1 on error
 *
 * Adds or updates tag key,value pairs in a header line.
 * Eg. for adding M5 tags to @SQ lines or updating sort order for the
 * @HD line.
 *
 * Specify multiple key,value pairs ending in NULL. Eg.
 * sam_hdr_update_line(h, "RG", "ID", "rg1", "DS", "description", "PG", "samtools", NULL)
 *
 * Attempting to update the record name (i.e. @SQ SN or @RG ID) will
 * work as long as the new name is not already in use, however doing this
 * on a file opened for reading may produce unexpected results.
 *
 * Renaming an @RG record in this way will only change the header.  Alignment
 * records written later will not be updated automatically even if they
 * reference the old read group name.
 *
 * Attempting to change an @PG ID tag is not permitted.
 */
HTSLIB_EXPORT
int sam_hdr_update_line(sam_hdr_t *h, const char *type,
        const char *ID_key, const char *ID_value, ...);

/// Remove all lines of a given type from a header, except the one matching an ID
/*!
 * @param type      Type of the searched line. Eg. "SQ"
 * @param ID_key    Tag key defining the line. Eg. "SN"
 * @param ID_value  Tag value associated with the key above. Eg. "ref1"
 * @return          0 on success, -1 on failure
 *
 * Remove all lines of type <type> from the header, except the one
 * specified by tag:value, i.e. the @SQ line containing "SN:ref1".
 *
 * If no line matches the key:value ID, all lines of the given type are removed.
 * To remove all lines of a given type, use NULL for both ID_key and ID_value.
 */
HTSLIB_EXPORT
int sam_hdr_remove_except(sam_hdr_t *h, const char *type, const char *ID_key, const char *ID_value);

HTSLIB_EXPORT
int sam_hdr_remove_lines(sam_hdr_t *h, const char *type, const char *id, void *rh);

/// Count the number of lines for a given header type
/*!
 * @param h     BAM header
 * @param type  Header type to count. Eg. "RG"
 * @return  Number of lines of this type on success; -1 on failure
 */
HTSLIB_EXPORT
int sam_hdr_count_lines(sam_hdr_t *h, const char *type);

/// Index of the line for the types that have dedicated look-up tables (SQ, RG, PG)
/*!
 * @param h     BAM header
 * @param type  Type of the searched line. Eg. "RG"
 * @param key   The value of the identifying key. Eg. "rg1"
 * @return  0-based index on success; -1 if line does not exist; -2 on failure
 */
HTSLIB_EXPORT
int sam_hdr_line_index(sam_hdr_t *bh, const char *type, const char *key);

/// Id key of the line for the types that have dedicated look-up tables (SQ, RG, PG)
/*!
 * @param h     BAM header
 * @param type  Type of the searched line. Eg. "RG"
 * @param pos   Zero-based index inside the type group. Eg. 2 (for the third RG line)
 * @return  Valid key string on success; NULL on failure
 */
HTSLIB_EXPORT
const char *sam_hdr_line_name(sam_hdr_t *bh, const char *type, int pos);

/* ==== Key:val level methods ==== */

/// Return the value associated with a key for a header line identified by ID_key:ID_val
/*!
 * @param type      Type of the line to which the tag belongs. Eg. "SQ"
 * @param ID_key    Tag key defining the line. Eg. "SN". Can be NULL, if looking for the first line.
 * @param ID_value  Tag value associated with the key above. Eg. "ref1". Can be NULL, if ID_key is NULL.
 * @param key       Key of the searched tag. Eg. "LN"
 * @param ks        kstring where the value will be written
 * @return          0 on success
 *                 -1 if the requested tag does not exist
 *                 -2 on other errors
 *
 * Looks for a specific key in a single SAM header line and writes the
 * associated value into @p ks.  The header line is selected using the ID_key
 * and ID_value parameters.  Any pre-existing content in @p ks will be
 * overwritten.
 */
HTSLIB_EXPORT
int sam_hdr_find_tag_id(sam_hdr_t *h, const char *type, const char *ID_key, const char *ID_value, const char *key, kstring_t *ks);

/// Return the value associated with a key for a header line identified by position
/*!
 * @param type      Type of the line to which the tag belongs. Eg. "SQ"
 * @param position  Index in lines of this type (zero-based). E.g. 3
 * @param key       Key of the searched tag. Eg. "LN"
 * @param ks        kstring where the value will be written
 * @return          0 on success
 *                 -1 if the requested tag does not exist
 *                 -2 on other errors
 *
 * Looks for a specific key in a single SAM header line and writes the
 * associated value into @p ks.  The header line is selected using the @p type
 * and @p position parameters.  Any pre-existing content in @p ks will be
 * overwritten.
 */
HTSLIB_EXPORT
int sam_hdr_find_tag_pos(sam_hdr_t *h, const char *type, int pos, const char *key, kstring_t *ks);

/// Remove the key from the line identified by type, ID_key and ID_value.
/*!
 * @param type      Type of the line to which the tag belongs. Eg. "SQ"
 * @param ID_key    Tag key defining the line. Eg. "SN"
 * @param ID_value  Tag value associated with the key above. Eg. "ref1"
 * @param key       Key of the targeted tag. Eg. "M5"
 * @return          1 if the key was removed; 0 if it was not present; -1 on error
 */
HTSLIB_EXPORT
int sam_hdr_remove_tag_id(sam_hdr_t *h, const char *type, const char *ID_key, const char *ID_value, const char *key);

/// Get the target id for a given reference sequence name
/*!
 * @param ref  Reference name
 * @return     Positive value on success,
 *             -1 if unknown reference,
 *             -2 if the header could not be parsed
 *
 * Looks up a reference sequence by name in the reference hash table
 * and returns the numerical target id.
 */
HTSLIB_EXPORT
int sam_hdr_name2tid(sam_hdr_t *h, const char *ref);

/// Alias of sam_hdr_name2tid(), for backwards compatibility.
/*!
 * @param ref  Reference name
 * @return     Positive value on success,
 *             -1 if unknown reference,
 *             -2 if the header could not be parsed
 */
static inline int bam_name2id(sam_hdr_t *h, const char *ref) { return sam_hdr_name2tid(h, ref); }

/// Generate a unique \@PG ID: value
/*!
 * @param name  Name of the program. Eg. samtools
 * @return      Valid ID on success, NULL on failure
 *
 * Returns a unique ID from a base name.  The string returned will remain
 * valid until the next call to this function, or the header is destroyed.
 * The caller should not attempt to free() or realloc() it.
 */
HTSLIB_EXPORT
const char *sam_hdr_pg_id(sam_hdr_t *h, const char *name);

/// Add an \@PG line.
/*!
 * @param name  Name of the program. Eg. samtools
 * @return      0 on success, -1 on failure
 *
 * If we wish complete control over this use sam_hdr_add_line() directly. This
 * function uses that, but attempts to do a lot of tedious house work for
 * you too.
 *
 * - It will generate a suitable ID if the supplied one clashes.
 * - It will generate multiple \@PG records if we have multiple PG chains.
 *
 * Call it as per sam_hdr_add_line() with a series of key,value pairs ending
 * in NULL.
 */
HTSLIB_EXPORT
int sam_hdr_add_pg(sam_hdr_t *h, const char *name, ...);


/*
 * Macros for changing the \@HD line. They eliminate the need to use NULL method arguments.
 */

/// Returns the SAM formatted text of the \@HD header line
#define sam_hdr_find_hd(h, ks) sam_hdr_find_line_id((h), "HD", NULL, NULL, (ks))
/// Returns the value associated with a given \@HD line tag
#define sam_hdr_find_tag_hd(h, key, ks) sam_hdr_find_tag_id((h), "HD", NULL, NULL, (key), (ks))
/// Adds or updates tags on the header \@HD line
#define sam_hdr_update_hd(h, ...) sam_hdr_update_line((h), "HD", NULL, NULL, __VA_ARGS__, NULL)
/// Removes the \@HD line tag with the given key
#define sam_hdr_remove_tag_hd(h, key) sam_hdr_remove_tag_id((h), "HD", NULL, NULL, (key))

/* Alignment */

/// Create a new bam1_t alignment structure
/**
   @return An empty bam1_t structure on success, NULL on failure

   The bam1_t struct returned by a successful call should be freed
   via bam_destroy1() when it is no longer needed.
 */
HTSLIB_EXPORT
bam1_t *bam_init1(void);

#define BAM_USER_OWNS_STRUCT 1
#define BAM_USER_OWNS_DATA   2

static inline void bam_set_mempolicy(bam1_t *b, uint32_t policy) {
    b->mempolicy = policy;
}

/// Get alignment record memory policy
/** @param b    Alignment record

    See bam_set_mempolicy()
 */
static inline uint32_t bam_get_mempolicy(bam1_t *b) {
    return b->mempolicy;
}

/// Sets all components of an alignment structure
/**
   @param bam      Target alignment structure. Must be initialized by a call to bam_init1().
                   The data field will be reallocated automatically as needed.
   @param l_qname  Length of the query name. If set to 0, the placeholder query name "*" will be used.
   @param qname    Query name, may be NULL if l_qname = 0
   @param flag     Bitwise flag, a combination of the BAM_F* constants.
   @param tid      Chromosome ID, defined by sam_hdr_t (a.k.a. RNAME).
   @param pos      0-based leftmost coordinate.
   @param mapq     Mapping quality.
   @param n_cigar  Number of CIGAR operations.
   @param cigar    CIGAR data, may be NULL if n_cigar = 0.
   @param mtid     Chromosome ID of next read in template, defined by sam_hdr_t (a.k.a. RNEXT).
   @param mpos     0-based leftmost coordinate of next read in template (a.k.a. PNEXT).
   @param isize    Observed template length ("insert size") (a.k.a. TLEN).
   @param l_seq    Length of the query sequence (read) and sequence quality string.
   @param seq      Sequence, may be NULL if l_seq = 0.
   @param qual     Sequence quality, may be NULL.
   @param l_aux    Length to be reserved for auxiliary field data, may be 0.

   @return >= 0 on success (number of bytes written to bam->data), negative (with errno set) on failure.
*/

/*! @function
 @abstract  Set the name of the query
 @param  b  pointer to an alignment
 @return    0 on success, -1 on failure
 */
HTSLIB_EXPORT
int bam_set_qname(bam1_t *b, const char *qname);

/*! @function
 @abstract  Parse a CIGAR string into a bam1_t struct
 @param  in      [in]  pointer to the source string
 @param  end     [out] address of the pointer to the new end of the input string
                       can be NULL
 @param  b       [in/out]  address of the destination bam1_t struct
 @return         number of processed CIGAR operators; -1 on error

 @discussion The BAM record may be partial and empty of existing cigar, seq
 and quality, as is the case during SAM parsing, or it may be an existing
 BAM record in which case this function replaces the existing CIGAR field
 and shuffles data accordingly.  A CIGAR of "*" will remove the CIGAR,
 returning zero.
 */
HTSLIB_EXPORT
ssize_t bam_parse_cigar(const char *in, char **end, bam1_t *b);

/*************************
 *** BAM/CRAM indexing ***
 *************************/

// These BAM iterator functions work only on BAM files.  To work with either
// BAM or CRAM files use the sam_index_load() & sam_itr_*() functions.
#define bam_itr_destroy(iter) hts_itr_destroy(iter)

// Load/build .csi or .bai BAM index file.  Does not work with CRAM.
// It is recommended to use the sam_index_* functions below instead.
#define bam_index_load(fn) hts_idx_load((fn), HTS_FMT_BAI)

/// Writes the index initialised with sam_idx_init to disk.
/** @param fp        File handle for the data file being written.
    @return          0 on success, <0 on failure.
*/
HTSLIB_EXPORT
int sam_idx_save(htsFile *fp) HTS_RESULT_USED;

    /***************
     *** SAM I/O ***
     ***************/

    #define sam_open(fn, mode) (hts_open((fn), (mode)))
    #define sam_open_format(fn, mode, fmt) (hts_open_format((fn), (mode), (fmt)))
    #define sam_flush(fp) hts_flush((fp))
    #define sam_close(fp) hts_close(fp)

    HTSLIB_EXPORT
    int sam_hdr_change_HD(sam_hdr_t *h, const char *key, const char *val);


// Forward declaration, see hts_expr.h for full.
struct hts_filter_t;

/// sam_passes_filter - Checks whether a record passes an hts_filter.
/** @param h      Pointer to the header structure previously read
 *  @param b      Pointer to the BAM record to be checked
 *  @param filt   Pointer to the filter, created from hts_filter_init.
 *  @return       1 if passes, 0 if not, and <0 on error.
 */
HTSLIB_EXPORT
int sam_passes_filter(const sam_hdr_t *h, const bam1_t *b,
                      struct hts_filter_t *filt);

    /*************************************
     *** Manipulating auxiliary fields ***
     *************************************/

/// Converts a BAM aux tag to SAM format
/*
 * @param key  Two letter tag key
 * @param type Single letter type code: ACcSsIifHZB.
 * @param tag  Tag data pointer, in BAM format
 * @param end  Pointer to end of bam record (largest extent of tag)
 * @param ks   kstring to write the formatted tag to
 *
 * @return pointer to end of tag on success,
 *         NULL on failure.
 *
 * @discussion The three separate parameters key, type, tag may be
 * derived from a s=bam_aux_get() query as s-2, *s and s+1.  However
 * it is recommended to use bam_aux_get_str in this situation.
 * The desire to split these parameters up is for potential processing
 * of non-BAM formats that encode using a BAM type mechanism
 * (such as the internal CRAM representation).
 */
static inline const uint8_t *sam_format_aux1(const uint8_t *key,
                                             const uint8_t type,
                                             const uint8_t *tag,
                                             const uint8_t *end,
                                             kstring_t *ks) {
    int r = 0;
    const uint8_t *s = tag; // brevity and consistency with other code.
    r |= kputsn_((char*)key, 2, ks) < 0;
    r |= kputc_(':', ks) < 0;
    if (type == 'C') {
        r |= kputsn_("i:", 2, ks) < 0;
        r |= kputw(*s, ks) < 0;
        ++s;
    } else if (type == 'c') {
        r |= kputsn_("i:", 2, ks) < 0;
        r |= kputw(*(int8_t*)s, ks) < 0;
        ++s;
    } else if (type == 'S') {
        if (end - s >= 2) {
            r |= kputsn_("i:", 2, ks) < 0;
            r |= kputuw(le_to_u16(s), ks) < 0;
            s += 2;
        } else goto bad_aux;
    } else if (type == 's') {
        if (end - s >= 2) {
            r |= kputsn_("i:", 2, ks) < 0;
            r |= kputw(le_to_i16(s), ks) < 0;
            s += 2;
        } else goto bad_aux;
    } else if (type == 'I') {
        if (end - s >= 4) {
            r |= kputsn_("i:", 2, ks) < 0;
            r |= kputuw(le_to_u32(s), ks) < 0;
            s += 4;
        } else goto bad_aux;
    } else if (type == 'i') {
        if (end - s >= 4) {
            r |= kputsn_("i:", 2, ks) < 0;
            r |= kputw(le_to_i32(s), ks) < 0;
            s += 4;
        } else goto bad_aux;
    } else if (type == 'A') {
        r |= kputsn_("A:", 2, ks) < 0;
        r |= kputc_(*s, ks) < 0;
        ++s;
    } else if (type == 'f') {
        if (end - s >= 4) {
            // cast to avoid triggering -Wdouble-promotion
            ksprintf(ks, "f:%g", (double)le_to_float(s));
            s += 4;
        } else goto bad_aux;

    } else if (type == 'd') {
        // NB: "d" is not an official type in the SAM spec.
        // However for unknown reasons samtools has always supported this.
        // We believe, HOPE, it is not in general usage and we do not
        // encourage it.
        if (end - s >= 8) {
            ksprintf(ks, "d:%g", le_to_double(s));
            s += 8;
        } else goto bad_aux;
    } else if (type == 'Z' || type == 'H') {
        r |= kputc_(type, ks) < 0;
        r |= kputc_(':', ks) < 0;
        while (s < end && *s) r |= kputc_(*s++, ks) < 0;
        r |= kputsn("", 0, ks) < 0;     //ensures NUL termination
        if (s >= end)
            goto bad_aux;
        ++s;
    } else if (type == 'B') {
        uint8_t sub_type = *(s++);
        unsigned sub_type_size;

        // or externalise sam.c's aux_type2size function?
        switch (sub_type) {
        case 'A': case 'c': case 'C':
            sub_type_size = 1;
            break;
        case 's': case 'S':
            sub_type_size = 2;
            break;
        case 'i': case 'I': case 'f':
            sub_type_size = 4;
            break;
        default:
            sub_type_size = 0;
            break;
        }

        uint32_t i, n;
        if (sub_type_size == 0 || end - s < 4)
            goto bad_aux;
        n = le_to_u32(s);
        s += 4; // now points to the start of the array
        if ((size_t)(end - s) / sub_type_size < n)
            goto bad_aux;
        r |= kputsn_("B:", 2, ks) < 0;
        r |= kputc(sub_type, ks) < 0; // write the type
        switch (sub_type) {
        case 'c':
            if (ks_expand(ks, n*2) < 0) goto mem_err;
            for (i = 0; i < n; ++i) {
                ks->s[ks->l++] = ',';
                r |= kputw(*(int8_t*)s, ks) < 0;
                ++s;
            }
            break;
        case 'C':
            if (ks_expand(ks, n*2) < 0) goto mem_err;
            for (i = 0; i < n; ++i) {
                ks->s[ks->l++] = ',';
                r |= kputuw(*(uint8_t*)s, ks) < 0;
                ++s;
            }
            break;
        case 's':
            if (ks_expand(ks, n*4) < 0) goto mem_err;
            for (i = 0; i < n; ++i) {
                ks->s[ks->l++] = ',';
                r |= kputw(le_to_i16(s), ks) < 0;
                s += 2;
            }
            break;
        case 'S':
            if (ks_expand(ks, n*4) < 0) goto mem_err;
            for (i = 0; i < n; ++i) {
                ks->s[ks->l++] = ',';
                r |= kputuw(le_to_u16(s), ks) < 0;
                s += 2;
            }
            break;
        case 'i':
            if (ks_expand(ks, n*6) < 0) goto mem_err;
            for (i = 0; i < n; ++i) {
                ks->s[ks->l++] = ',';
                r |= kputw(le_to_i32(s), ks) < 0;
                s += 4;
            }
            break;
        case 'I':
            if (ks_expand(ks, n*6) < 0) goto mem_err;
            for (i = 0; i < n; ++i) {
                ks->s[ks->l++] = ',';
                r |= kputuw(le_to_u32(s), ks) < 0;
                s += 4;
            }
            break;
        case 'f':
            if (ks_expand(ks, n*8) < 0) goto mem_err;
            for (i = 0; i < n; ++i) {
                ks->s[ks->l++] = ',';
                // cast to avoid triggering -Wdouble-promotion
                r |= kputd((double)le_to_float(s), ks) < 0;
                s += 4;
            }
            break;
        default:
            goto bad_aux;
        }
    } else { // Unknown type
        goto bad_aux;
    }
    return r ? NULL : s;

 bad_aux:
    errno = EINVAL;
    return NULL;

 mem_err:
    hts_log_error("Out of memory");
    errno = ENOMEM;
    return NULL;
}

/// Delete tag data from a bam record
/** @param b   The BAM record to update
    @param s   Pointer to the aux field to delete, as returned by bam_aux_get()
               Must not be NULL
    @return    0 on success; -1 on failure

If the BAM record's aux data is corrupt, errno is set to EINVAL and this
function returns -1.
*/
HTSLIB_EXPORT
int bam_aux_del(bam1_t *b, uint8_t *s);

/// Delete an aux field from a BAM record
/** @param b   The BAM record to update
    @param s   Pointer to the aux field to delete, as returned by
               bam_aux_first()/_next()/_get(); must not be NULL
    @return    Pointer to the following aux field, or NULL if none or on error

Identical to @c bam_aux_del() apart from the return value, which is an
aux iterator suitable for use with @c bam_aux_next()/etc.

Whenever NULL is returned, errno will also be set: ENOENT if the aux field
deleted was the record's last one; otherwise EINVAL, indicating that the
BAM record's aux data is corrupt.
 */
HTSLIB_EXPORT
uint8_t *bam_aux_remove(bam1_t *b, uint8_t *s);

/**************************
 *** Pileup and Mpileup ***
 **************************/

#if !defined(BAM_NO_PILEUP)

/*! @typedef
 @abstract Generic pileup 'client data'.

 @discussion The pileup iterator allows setting a constructor and
 destructor function, which will be called every time a sequence is
 fetched and discarded.  This permits caching of per-sequence data in
 a tidy manner during the pileup process.  This union is the cached
 data to be manipulated by the "client" (the caller of pileup).
*/
typedef union {
    void *p;
    int64_t i;
    double f;
} bam_pileup_cd;

/*! @typedef
 @abstract Structure for one alignment covering the pileup position.
 @field  b          pointer to the alignment
 @field  qpos       position of the read base at the pileup site, 0-based
 @field  indel      indel length; 0 for no indel, positive for ins and negative for del
 @field  level      the level of the read in the "viewer" mode
 @field  is_del     1 iff the base on the padded read is a deletion
 @field  is_head    1 iff this is the first base in the query sequence
 @field  is_tail    1 iff this is the last base in the query sequence
 @field  is_refskip 1 iff the base on the padded read is part of CIGAR N op
 @field  aux        (used by bcf_call_gap_prep())
 @field  cigar_ind  index of the CIGAR operator that has just been processed

 @discussion See also bam_plbuf_push() and bam_lplbuf_push(). The
 difference between the two functions is that the former does not
 set bam_pileup1_t::level, while the later does. Level helps the
 implementation of alignment viewers, but calculating this has some
 overhead.
 */
typedef struct bam_pileup1_t {
    bam1_t *b;
    int32_t qpos;
    int indel, level;
    uint32_t is_del:1, is_head:1, is_tail:1, is_refskip:1, /* reserved */ :1, aux:27;
    bam_pileup_cd cd; // generic per-struct data, owned by caller.
    int cigar_ind;
} bam_pileup1_t;

typedef int (*bam_plp_auto_f)(void *data, bam1_t *b);

struct bam_plp_s;
typedef struct bam_plp_s *bam_plp_t;

struct bam_mplp_s;
typedef struct bam_mplp_s *bam_mplp_t;

    /// Get pileup padded insertion sequence
    /**
     * @param p       pileup data
     * @param ins     the kstring where the insertion sequence will be written
     * @param del_len location for deletion length
     * @return the length of insertion string on success; -1 on failure.
     *
     * Fills out the kstring with the padded insertion sequence for the current
     * location in 'p'.  If this is not an insertion site, the string is blank.
     *
     * If del_len is not NULL, the location pointed to is set to the length of
     * any deletion immediately following the insertion, or zero if none.
     */
    HTSLIB_EXPORT
    int bam_plp_insertion(const bam_pileup1_t *p, kstring_t *ins, int *del_len) HTS_RESULT_USED;


    /*! @typedef
     @abstract An opaque type used for caching base modification state between
     successive calls to bam_mods_* functions.
    */
    typedef struct hts_base_mod_state hts_base_mod_state;

    /// Create a new bam_mplp_t structure
    /** The struct returned by a successful call should be freed
     *  via bam_mplp_destroy() when it is no longer needed.
     */
    HTSLIB_EXPORT
    bam_mplp_t bam_mplp_init(int n, bam_plp_auto_f func, void **data);

    /// Set up mpileup overlap detection
    /**
     * @param iter    mpileup iterator
     * @return 0 on success; a negative value on error
     *
     *  If called, mpileup will detect overlapping
     *  read pairs and for each base pair set the base quality of the
     *  lower-quality base to zero, thus effectively discarding it from
     *  calling. If the two bases are identical, the quality of the other base
     *  is increased to the sum of their qualities (capped at 200), otherwise
     *  it is multiplied by 0.8.
     */
    HTSLIB_EXPORT
    int bam_mplp_init_overlaps(bam_mplp_t iter);

    HTSLIB_EXPORT
    void bam_mplp_destroy(bam_mplp_t iter);

    HTSLIB_EXPORT
    void bam_mplp_set_maxcnt(bam_mplp_t iter, int maxcnt);

    HTSLIB_EXPORT
    int bam_mplp_auto(bam_mplp_t iter, int *_tid, int *_pos, int *n_plp, const bam_pileup1_t **plp);

    HTSLIB_EXPORT
    int bam_mplp64_auto(bam_mplp_t iter, int *_tid, hts_pos_t *_pos, int *n_plp, const bam_pileup1_t **plp);

    HTSLIB_EXPORT
    void bam_mplp_reset(bam_mplp_t iter);

    HTSLIB_EXPORT
    void bam_mplp_constructor(bam_mplp_t iter,
                              int (*func)(void *data, const bam1_t *b, bam_pileup_cd *cd));

    HTSLIB_EXPORT
    void bam_mplp_destructor(bam_mplp_t iter,
                             int (*func)(void *data, const bam1_t *b, bam_pileup_cd *cd));

#endif // ~!defined(BAM_NO_PILEUP)


/***********************************
 * BAQ calculation and realignment *
 ***********************************/

HTSLIB_EXPORT
int sam_cap_mapq(bam1_t *b, const char *ref, hts_pos_t ref_len, int thres);

// Used as flag parameter in sam_prob_realn.
enum htsRealnFlags {
    BAQ_APPLY = 1,
    BAQ_EXTEND = 2,
    BAQ_REDO = 4,

    // Platform subfield, in bit position 3 onwards
    BAQ_AUTO = 0<<3,
    BAQ_ILLUMINA = 1<<3,
    BAQ_PACBIOCCS = 2<<3,
    BAQ_PACBIO = 3<<3,
    BAQ_ONT = 4<<3,
    BAQ_GENAPSYS = 5<<3
};

// ---------------------------
// Base modification retrieval

/*! @typedef
 @abstract Holds a single base modification.
 @field modified_base     The short base code (m, h, etc) or -ChEBI (negative)
 @field canonical_base    The canonical base referred to in the MM tag.
                          One of A, C, G, T or N.  Note this may not be the
                          explicit base recorded in the SEQ column (esp. if N).
 @field stran             0 or 1, indicating + or - strand from MM tag.
 @field qual              Quality code (256*probability), or -1 if unknown

 @discussion
 Note this doesn't hold any location data or information on which other
 modifications may be possible at this site.
*/
typedef struct hts_base_mod {
    int modified_base;
    int canonical_base;
    int strand;
    int qual;
} hts_base_mod;

#define HTS_MOD_UNKNOWN   -1  // In MM but no ML
#define HTS_MOD_UNCHECKED -2  // Not in MM and in explicit mode

// Flags for hts_parse_basemod2
#define HTS_MOD_REPORT_UNCHECKED 1

/// Allocates an hts_base_mode_state.
/**
 * @return An hts_base_mode_state pointer on success,
 *         NULL on failure.
 *
 * This just allocates the memory.  The initialisation of the contents is
 * done using bam_parse_basemod.  Successive calls may be made to that
 * without the need to free and allocate a new state.
 *
 * The state be destroyed using the hts_base_mode_state_free function.
 */
HTSLIB_EXPORT
hts_base_mod_state *hts_base_mod_state_alloc(void);

/// Destroys an  hts_base_mode_state.
/**
 * @param state    The base modification state pointer.
 *
 * The should have previously been created by hts_base_mode_state_alloc.
 */
HTSLIB_EXPORT
void hts_base_mod_state_free(hts_base_mod_state *state);

/// Parses the Mm and Ml tags out of a bam record.
/**
 * @param b        BAM alignment record
 * @param state    The base modification state pointer.
 * @return 0 on success,
 *         -1 on failure.
 *
 * This fills out the contents of the modification state, resetting the
 * iterator location to the first sequence base.
 */
HTSLIB_EXPORT
int bam_parse_basemod(const bam1_t *b, hts_base_mod_state *state);

/// Parses the Mm and Ml tags out of a bam record.
/**
 * @param b        BAM alignment record
 * @param state    The base modification state pointer.
 * @param flags    A bit-field controlling base modification processing
 *
 * @return 0 on success,
 *         -1 on failure.
 *
 * This fills out the contents of the modification state, resetting the
 * iterator location to the first sequence base.
 */
HTSLIB_EXPORT
int bam_parse_basemod2(const bam1_t *b, hts_base_mod_state *state,
                       uint32_t flags);

/// Finds the next location containing base modifications and returns them
/**
 * @param b        BAM alignment record
 * @param state    The base modification state pointer.
 * @param mods     A supplied array for returning base modifications
 * @param n_mods   The size of the mods array
 * @param pos      Pointer holding position of modification in sequence
 * @return The number of modifications found on success,
 *         0 if no more modifications are present,
 *         -1 on failure.
 *
 * Unlike bam_mods_at_next_pos this skips ahead to the next site
 * with modifications.
 *
 * If more than n_mods modifications are found, the total found is returned.
 * Note this means the caller needs to check whether this is higher than
 * n_mods.
 */
HTSLIB_EXPORT
int bam_next_basemod(const bam1_t *b, hts_base_mod_state *state,
                     hts_base_mod *mods, int n_mods, int *pos);


#ifdef __cplusplus
}
#endif

#ifdef HTSLIB_SSIZE_T
#undef HTSLIB_SSIZE_T
#undef ssize_t
#endif

#endif
