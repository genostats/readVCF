/* The MIT License

   Copyright (c) 2008 Broad Institute / Massachusetts Institute of Technology
                 2011, 2012 Attractive Chaos <attractor@live.co.uk>
   Copyright (C) 2009, 2013-2022 Genome Research Ltd

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
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
*/

#define HTS_BUILDING_LIBRARY // Enables HTSLIB_EXPORT, see htslib/hts_defs.h
#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <assert.h>
#include <pthread.h>
#include <sys/types.h>
#include <inttypes.h>
#include <zlib.h>

#include "htslib/hts.h"
#include "htslib/bgzf.h"
#include "htslib/hfile.h"
#include "htslib/thread_pool.h"
#include "htslib/hts_endian.h"
#include "cram/pooled_alloc.h"
#include "hts_internal.h"

#ifndef EFTYPE
#define EFTYPE ENOEXEC
#endif

#define BGZF_CACHE
#define BGZF_MT

#define BLOCK_HEADER_LENGTH 18
#define BLOCK_FOOTER_LENGTH 8


/* BGZF/GZIP header (specialized from RFC 1952; little endian):
 +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
 | 31|139|  8|  4|              0|  0|255|      6| 66| 67|      2|BLK_LEN|
 +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
  BGZF extension:
                ^                              ^   ^   ^
                |                              |   |   |
               FLG.EXTRA                     XLEN  B   C

  BGZF format is compatible with GZIP. It limits the size of each compressed
  block to 2^16 bytes and adds and an extra "BC" field in the gzip header which
  records the size.

*/
static const uint8_t g_magic[19] = "\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\0\0";

#ifdef BGZF_CACHE
typedef struct {
    int size;
    uint8_t *block;
    int64_t end_offset;
} cache_t;

#include "htslib/khash.h"
KHASH_MAP_INIT_INT64(cache, cache_t)
#endif

struct bgzf_cache_t {
    khash_t(cache) *h;
    khint_t last_pos;
};

#ifdef BGZF_MT

typedef struct bgzf_job {
    BGZF *fp;
    unsigned char comp_data[BGZF_MAX_BLOCK_SIZE];
    size_t comp_len;
    unsigned char uncomp_data[BGZF_MAX_BLOCK_SIZE];
    size_t uncomp_len;
    int errcode;
    int64_t block_address;
    int hit_eof;
} bgzf_job;

enum mtaux_cmd {
    NONE = 0,
    SEEK,
    SEEK_DONE,
    HAS_EOF,
    HAS_EOF_DONE,
    CLOSE,
};

// When multi-threaded bgzf_tell won't work, so we delay the hts_idx_push
// until we've written the last block.
typedef struct {
    hts_pos_t beg, end;
    int tid, is_mapped;  // args for hts_idx_push
    uint64_t offset, block_number;
} hts_idx_cache_entry;

typedef struct {
    int nentries, mentries; // used and allocated
    hts_idx_cache_entry *e; // hts_idx elements
} hts_idx_cache_t;

typedef struct bgzf_mtaux_t {
    // Memory pool for bgzf_job structs, to avoid many malloc/free
    pool_alloc_t *job_pool;
    bgzf_job *curr_job;

    // Thread pool
    int n_threads;
    int own_pool;
    hts_tpool *pool;

    // Output queue holding completed bgzf_jobs
    hts_tpool_process *out_queue;

    // I/O thread.
    pthread_t io_task;
    pthread_mutex_t job_pool_m;
    int jobs_pending; // number of jobs waiting
    int flush_pending;
    void *free_block;
    int hit_eof;  // r/w entirely within main thread

    // Message passing to the reader thread; eg seek requests
    int errcode;
    uint64_t block_address;
    int eof;
    pthread_mutex_t command_m; // Set whenever fp is being updated
    pthread_cond_t command_c;
    enum mtaux_cmd command;

    // For multi-threaded on-the-fly indexing. See bgzf_idx_push below.
    pthread_mutex_t idx_m;
    hts_idx_t *hts_idx;
    uint64_t block_number, block_written;
    hts_idx_cache_t idx_cache;
} mtaux_t;
#endif

typedef struct
{
    uint64_t uaddr;  // offset w.r.t. uncompressed data
    uint64_t caddr;  // offset w.r.t. compressed data
}
bgzidx1_t;

struct bgzidx_t
{
    int noffs, moffs;       // the size of the index, n:used, m:allocated
    bgzidx1_t *offs;        // offsets
    uint64_t ublock_addr;   // offset of the current block (uncompressed data)
};

static int bgzf_idx_flush(BGZF *fp) {
    mtaux_t *mt = fp->mt;

    if (!mt->idx_cache.e) {
        mt->block_written++;
        return 0;
    }

    pthread_mutex_lock(&mt->idx_m);

    hts_idx_cache_entry *e = mt->idx_cache.e;
    int i;

    assert(mt->idx_cache.nentries == 0 || mt->block_written <= e[0].block_number);

    for (i = 0; i < mt->idx_cache.nentries && e[i].block_number == mt->block_written; i++) {
        if (hts_idx_push(mt->hts_idx, e[i].tid, e[i].beg, e[i].end,
                         (mt->block_address << 16) + e[i].offset,
                         e[i].is_mapped) < 0) {
            pthread_mutex_unlock(&mt->idx_m);
            return -1;
        }
    }

    memmove(&e[0], &e[i], (mt->idx_cache.nentries - i) * sizeof(*e));
    mt->idx_cache.nentries -= i;
    mt->block_written++;

    pthread_mutex_unlock(&mt->idx_m);
    return 0;
}

void bgzf_index_destroy(BGZF *fp);
int bgzf_index_add_block(BGZF *fp);
static int mt_destroy(mtaux_t *mt);

static inline void packInt16(uint8_t *buffer, uint16_t value)
{
    buffer[0] = value;
    buffer[1] = value >> 8;
}

static inline int unpackInt16(const uint8_t *buffer)
{
    return buffer[0] | buffer[1] << 8;
}

static inline void packInt32(uint8_t *buffer, uint32_t value)
{
    buffer[0] = value;
    buffer[1] = value >> 8;
    buffer[2] = value >> 16;
    buffer[3] = value >> 24;
}

static void razf_info(hFILE *hfp, const char *filename)
{
    uint64_t usize, csize;
    off_t sizes_pos;

    if (filename == NULL || strcmp(filename, "-") == 0) filename = "FILE";

    // RAZF files end with USIZE,CSIZE stored as big-endian uint64_t
    if ((sizes_pos = hseek(hfp, -16, SEEK_END)) < 0) goto no_sizes;
    if (hread(hfp, &usize, 8) != 8 || hread(hfp, &csize, 8) != 8) goto no_sizes;
    if (!ed_is_big()) ed_swap_8p(&usize), ed_swap_8p(&csize);
    if (csize >= sizes_pos) goto no_sizes; // Very basic validity check

    hts_log_error(
"To decompress this file, use the following commands:\n"
"    truncate -s %" PRIu64 " %s\n"
"    gunzip %s\n"
"The resulting uncompressed file should be %" PRIu64 " bytes in length.\n"
"If you do not have a truncate command, skip that step (though gunzip will\n"
"likely produce a \"trailing garbage ignored\" message, which can be ignored).",
                  csize, filename, filename, usize);
    return;

no_sizes:
    hts_log_error(
"To decompress this file, use the following command:\n"
"    gunzip %s\n"
"This will likely produce a \"trailing garbage ignored\" message, which can\n"
"usually be safely ignored.", filename);
}

static const char *bgzf_zerr(int errnum, z_stream *zs)
{
    static char buffer[32];

    /* Return zs->msg if available.
       zlib doesn't set this very reliably.  Looking at the source suggests
       that it may get set to a useful message for deflateInit2, inflateInit2
       and inflate when it returns Z_DATA_ERROR. For inflate with other
       return codes, deflate, deflateEnd and inflateEnd it doesn't appear
       to be useful.  For the likely non-useful cases, the caller should
       pass NULL into zs. */

    if (zs && zs->msg) return zs->msg;

    // gzerror OF((gzFile file, int *errnum)
    switch (errnum) {
    case Z_ERRNO:
        return strerror(errno);
    case Z_STREAM_ERROR:
        return "invalid parameter/compression level, or inconsistent stream state";
    case Z_DATA_ERROR:
        return "invalid or incomplete IO";
    case Z_MEM_ERROR:
        return "out of memory";
    case Z_BUF_ERROR:
        return "progress temporarily not possible, or in() / out() returned an error";
    case Z_VERSION_ERROR:
        return "zlib version mismatch";
    case Z_NEED_DICT:
        return "data was compressed using a dictionary";
    case Z_OK: // 0: maybe gzgets error Z_NULL
    default:
        snprintf(buffer, sizeof(buffer), "[%d] unknown", errnum);
        return buffer;  // FIXME: Not thread-safe.
    }
}

static BGZF *bgzf_read_init(hFILE *hfpr, const char *filename)
{
    BGZF *fp;
    uint8_t magic[18];
    ssize_t n = hpeek(hfpr, magic, 18);
    if (n < 0) return NULL;

    fp = (BGZF*)calloc(1, sizeof(BGZF));
    if (fp == NULL) return NULL;

    fp->is_write = 0;
    fp->uncompressed_block = malloc(2 * BGZF_MAX_BLOCK_SIZE);
    if (fp->uncompressed_block == NULL) { free(fp); return NULL; }
    fp->compressed_block = (char *)fp->uncompressed_block + BGZF_MAX_BLOCK_SIZE;
    fp->is_compressed = (n==18 && magic[0]==0x1f && magic[1]==0x8b);
    fp->is_gzip = ( !fp->is_compressed || ((magic[3]&4) && memcmp(&magic[12], "BC\2\0",4)==0) ) ? 0 : 1;
    if (fp->is_compressed && (magic[3]&4) && memcmp(&magic[12], "RAZF", 4)==0) {
        hts_log_error("Cannot decompress legacy RAZF format");
        razf_info(hfpr, filename);
        free(fp->uncompressed_block);
        free(fp);
        errno = EFTYPE;
        return NULL;
    }
#ifdef BGZF_CACHE
    if (!(fp->cache = malloc(sizeof(*fp->cache)))) {
        free(fp->uncompressed_block);
        free(fp);
        return NULL;
    }
    if (!(fp->cache->h = kh_init(cache))) {
        free(fp->uncompressed_block);
        free(fp->cache);
        free(fp);
        return NULL;
    }
    fp->cache->last_pos = 0;
#endif
    return fp;
}

// get the compress level from the mode string: compress_level==-1 for the default level, -2 plain uncompressed
static int mode2level(const char *mode)
{
    int i, compress_level = -1;
    for (i = 0; mode[i]; ++i)
        if (mode[i] >= '0' && mode[i] <= '9') break;
    if (mode[i]) compress_level = (int)mode[i] - '0';
    if (strchr(mode, 'u')) compress_level = -2;
    return compress_level;
}
static BGZF *bgzf_write_init(const char *mode)
{
    BGZF *fp;
    fp = (BGZF*)calloc(1, sizeof(BGZF));
    if (fp == NULL) goto mem_fail;
    fp->is_write = 1;
    int compress_level = mode2level(mode);
    if ( compress_level==-2 )
    {
        fp->is_compressed = 0;
        return fp;
    }
    fp->is_compressed = 1;

    fp->uncompressed_block = malloc(2 * BGZF_MAX_BLOCK_SIZE);
    if (fp->uncompressed_block == NULL) goto mem_fail;
    fp->compressed_block = (char *)fp->uncompressed_block + BGZF_MAX_BLOCK_SIZE;

    fp->compress_level = compress_level < 0? Z_DEFAULT_COMPRESSION : compress_level; // Z_DEFAULT_COMPRESSION==-1
    if (fp->compress_level > 9) fp->compress_level = Z_DEFAULT_COMPRESSION;
    if ( strchr(mode,'g') )
    {
        // gzip output
        fp->is_gzip = 1;
        fp->gz_stream = (z_stream*)calloc(1,sizeof(z_stream));
        if (fp->gz_stream == NULL) goto mem_fail;
        fp->gz_stream->zalloc = NULL;
        fp->gz_stream->zfree  = NULL;
        fp->gz_stream->msg    = NULL;

        int ret = deflateInit2(fp->gz_stream, fp->compress_level, Z_DEFLATED, 15|16, 8, Z_DEFAULT_STRATEGY);
        if (ret!=Z_OK) {
            hts_log_error("Call to deflateInit2 failed: %s", bgzf_zerr(ret, fp->gz_stream));
            goto fail;
        }
    }
    return fp;

mem_fail:
    hts_log_error("%s", strerror(errno));

fail:
    if (fp != NULL) {
        free(fp->uncompressed_block);
        free(fp->gz_stream);
        free(fp);
    }
    return NULL;
}

BGZF *bgzf_open(const char *path, const char *mode)
{
    BGZF *fp = 0;
    if (strchr(mode, 'r')) {
        hFILE *fpr;
        if ((fpr = hopen(path, mode)) == 0) return 0;
        fp = bgzf_read_init(fpr, path);
        if (fp == 0) { hclose_abruptly(fpr); return NULL; }
        fp->fp = fpr;
    } else if (strchr(mode, 'w') || strchr(mode, 'a')) {
        hFILE *fpw;
        if ((fpw = hopen(path, mode)) == 0) return 0;
        fp = bgzf_write_init(mode);
        if (fp == NULL) return NULL;
        fp->fp = fpw;
    }
    else { errno = EINVAL; return 0; }

    fp->is_be = ed_is_big();
    return fp;
}

BGZF *bgzf_dopen(int fd, const char *mode)
{
    BGZF *fp = 0;
    if (strchr(mode, 'r')) {
        hFILE *fpr;
        if ((fpr = hdopen(fd, mode)) == 0) return 0;
        fp = bgzf_read_init(fpr, NULL);
        if (fp == 0) { hclose_abruptly(fpr); return NULL; } // FIXME this closes fd
        fp->fp = fpr;
    } else if (strchr(mode, 'w') || strchr(mode, 'a')) {
        hFILE *fpw;
        if ((fpw = hdopen(fd, mode)) == 0) return 0;
        fp = bgzf_write_init(mode);
        if (fp == NULL) return NULL;
        fp->fp = fpw;
    }
    else { errno = EINVAL; return 0; }

    fp->is_be = ed_is_big();
    return fp;
}

BGZF *bgzf_hopen(hFILE *hfp, const char *mode)
{
    BGZF *fp = NULL;
    if (strchr(mode, 'r')) {
        fp = bgzf_read_init(hfp, NULL);
        if (fp == NULL) return NULL;
    } else if (strchr(mode, 'w') || strchr(mode, 'a')) {
        fp = bgzf_write_init(mode);
        if (fp == NULL) return NULL;
    }
    else { errno = EINVAL; return 0; }

    fp->fp = hfp;
    fp->is_be = ed_is_big();
    return fp;
}

int bgzf_compress(void *_dst, size_t *dlen, const void *src, size_t slen, int level)
{
    uint32_t crc;
    z_stream zs;
    uint8_t *dst = (uint8_t*)_dst;

    if (level == 0) {
    uncomp:
        // Uncompressed data
        if (*dlen < slen+5 + BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH) return -1;
        dst[BLOCK_HEADER_LENGTH] = 1; // BFINAL=1, BTYPE=00; see RFC1951
        u16_to_le(slen,  &dst[BLOCK_HEADER_LENGTH+1]); // length
        u16_to_le(~slen, &dst[BLOCK_HEADER_LENGTH+3]); // ones-complement length
        memcpy(dst + BLOCK_HEADER_LENGTH+5, src, slen);
        *dlen = slen+5 + BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH;
    } else {
        // compress the body
        zs.zalloc = NULL; zs.zfree = NULL;
        zs.msg = NULL;
        zs.next_in  = (Bytef*)src;
        zs.avail_in = slen;
        zs.next_out = dst + BLOCK_HEADER_LENGTH;
        zs.avail_out = *dlen - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH;
        int ret = deflateInit2(&zs, level, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY); // -15 to disable zlib header/footer
        if (ret!=Z_OK) {
            hts_log_error("Call to deflateInit2 failed: %s", bgzf_zerr(ret, &zs));
            return -1;
        }
        if ((ret = deflate(&zs, Z_FINISH)) != Z_STREAM_END) {
            if (ret == Z_OK && zs.avail_out == 0) {
                deflateEnd(&zs);
                goto uncomp;
            } else {
                hts_log_error("Deflate operation failed: %s", bgzf_zerr(ret, ret == Z_DATA_ERROR ? &zs : NULL));
            }
            return -1;
        }
        // If we used up the entire output buffer, then we either ran out of
        // room or we *just* fitted, but either way we may as well store
        // uncompressed for faster decode.
        if (zs.avail_out == 0) {
            deflateEnd(&zs);
            goto uncomp;
        }
        if ((ret = deflateEnd(&zs)) != Z_OK) {
            hts_log_error("Call to deflateEnd failed: %s", bgzf_zerr(ret, NULL));
            return -1;
        }
        *dlen = zs.total_out + BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH;
    }

    // write the header
    memcpy(dst, g_magic, BLOCK_HEADER_LENGTH); // the last two bytes are a place holder for the length of the block
    packInt16(&dst[16], *dlen - 1); // write the compressed length; -1 to fit 2 bytes
    // write the footer
    crc = crc32(crc32(0L, NULL, 0L), (Bytef*)src, slen);
    packInt32((uint8_t*)&dst[*dlen - 8], crc);
    packInt32((uint8_t*)&dst[*dlen - 4], slen);
    return 0;
}

static int bgzf_gzip_compress(BGZF *fp, void *_dst, size_t *dlen, const void *src, size_t slen, int level)
{
    uint8_t *dst = (uint8_t*)_dst;
    z_stream *zs = fp->gz_stream;
    int flush = slen ? Z_PARTIAL_FLUSH : Z_FINISH;
    zs->next_in   = (Bytef*)src;
    zs->avail_in  = slen;
    zs->next_out  = dst;
    zs->avail_out = *dlen;
    int ret = deflate(zs, flush);
    if (ret == Z_STREAM_ERROR) {
        hts_log_error("Deflate operation failed: %s", bgzf_zerr(ret, NULL));
        return -1;
    }
    if (zs->avail_in != 0) {
        hts_log_error("Deflate block too large for output buffer");
        return -1;
    }
    *dlen = *dlen - zs->avail_out;
    return 0;
}

// Deflate the block in fp->uncompressed_block into fp->compressed_block. Also adds an extra field that stores the compressed block length.
static int deflate_block(BGZF *fp, int block_length)
{
    size_t comp_size = BGZF_MAX_BLOCK_SIZE;
    int ret;
    if ( !fp->is_gzip )
        ret = bgzf_compress(fp->compressed_block, &comp_size, fp->uncompressed_block, block_length, fp->compress_level);
    else
        ret = bgzf_gzip_compress(fp, fp->compressed_block, &comp_size, fp->uncompressed_block, block_length, fp->compress_level);

    if ( ret != 0 )
    {
        hts_log_debug("Compression error %d", ret);
        fp->errcode |= BGZF_ERR_ZLIB;
        return -1;
    }
    fp->block_offset = 0;
    return comp_size;
}

static int bgzf_uncompress(uint8_t *dst, size_t *dlen,
                           const uint8_t *src, size_t slen,
                           uint32_t expected_crc) {
    z_stream zs = {
        .zalloc = NULL,
        .zfree = NULL,
        .msg = NULL,
        .next_in = (Bytef*)src,
        .avail_in = slen,
        .next_out = (Bytef*)dst,
        .avail_out = *dlen
    };

    int ret = inflateInit2(&zs, -15);
    if (ret != Z_OK) {
        hts_log_error("Call to inflateInit2 failed: %s", bgzf_zerr(ret, &zs));
        return -1;
    }
    if ((ret = inflate(&zs, Z_FINISH)) != Z_STREAM_END) {
        hts_log_error("Inflate operation failed: %s", bgzf_zerr(ret, ret == Z_DATA_ERROR ? &zs : NULL));
        if ((ret = inflateEnd(&zs)) != Z_OK) {
            hts_log_warning("Call to inflateEnd failed: %s", bgzf_zerr(ret, NULL));
        }
        return -1;
    }
    if ((ret = inflateEnd(&zs)) != Z_OK) {
        hts_log_error("Call to inflateEnd failed: %s", bgzf_zerr(ret, NULL));
        return -1;
    }
    *dlen = *dlen - zs.avail_out;

    uint32_t crc = crc32(crc32(0L, NULL, 0L), (unsigned char *)dst, *dlen);
#ifdef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
    // Pretend the CRC was OK so the fuzzer doesn't have to get it right
    crc = expected_crc;
#endif
    if (crc != expected_crc) {
        hts_log_error("CRC32 checksum mismatch");
        return -2;
    }

    return 0;
}

// Inflate the block in fp->compressed_block into fp->uncompressed_block
static int inflate_block(BGZF* fp, int block_length)
{
    size_t dlen = BGZF_MAX_BLOCK_SIZE;
    uint32_t crc = le_to_u32((uint8_t *)fp->compressed_block + block_length-8);
    int ret = bgzf_uncompress(fp->uncompressed_block, &dlen,
                              (Bytef*)fp->compressed_block + 18,
                              block_length - 18, crc);
    if (ret < 0) {
        if (ret == -2)
            fp->errcode |= BGZF_ERR_CRC;
        else
            fp->errcode |= BGZF_ERR_ZLIB;
        return -1;
    }

    return dlen;
}

// Decompress the next part of a non-blocked GZIP file.
// Return the number of uncompressed bytes read, 0 on EOF, or a negative number on error.
// Will fill the output buffer unless the end of the GZIP file is reached.
static int inflate_gzip_block(BGZF *fp)
{
    // we will set this to true when we detect EOF, so we don't bang against the EOF more than once per call
    int input_eof = 0;

    // write to the part of the output buffer after block_offset
    fp->gz_stream->next_out = (Bytef*)fp->uncompressed_block + fp->block_offset;
    fp->gz_stream->avail_out = BGZF_MAX_BLOCK_SIZE - fp->block_offset;

    while ( fp->gz_stream->avail_out != 0 ) {
        // until we fill the output buffer (or hit EOF)

        if ( !input_eof && fp->gz_stream->avail_in == 0 ) {
            // we are out of input data in the buffer. Get more.
            fp->gz_stream->next_in = fp->compressed_block;
            int ret = hread(fp->fp, fp->compressed_block, BGZF_BLOCK_SIZE);
            if ( ret < 0 ) {
                // hread had an error. Pass it on.
                return ret;
            }
            fp->gz_stream->avail_in = ret;
            if ( fp->gz_stream->avail_in < BGZF_BLOCK_SIZE ) {
                // we have reached EOF but the decompressor hasn't necessarily
                input_eof = 1;
            }
        }

        fp->gz_stream->msg = NULL;
        // decompress as much data as we can
        int ret = inflate(fp->gz_stream, Z_SYNC_FLUSH);

        if ( (ret < 0 && ret != Z_BUF_ERROR) || ret == Z_NEED_DICT ) {
            // an error occurred, other than running out of space
            hts_log_error("Inflate operation failed: %s", bgzf_zerr(ret, ret == Z_DATA_ERROR ? fp->gz_stream : NULL));
            fp->errcode |= BGZF_ERR_ZLIB;
            return -1;
        } else if ( ret == Z_STREAM_END ) {
            // we finished a GZIP member

            // scratch for peeking to see if the file is over
            char c;
            if (fp->gz_stream->avail_in > 0 || hpeek(fp->fp, &c, 1) == 1) {
                // there is more data; try and read another GZIP member in the remaining data
                int reset_ret = inflateReset(fp->gz_stream);
                if (reset_ret != Z_OK) {
                    hts_log_error("Call to inflateReset failed: %s", bgzf_zerr(reset_ret, NULL));
                    fp->errcode |= BGZF_ERR_ZLIB;
                    return -1;
                }
            } else {
                // we consumed all the input data and hit Z_STREAM_END
                // so stop looping, even if we never fill the output buffer
                break;
            }
        } else if ( ret == Z_BUF_ERROR && input_eof && fp->gz_stream->avail_out > 0 ) {
            // the gzip file has ended prematurely
            hts_log_error("Gzip file truncated");
            fp->errcode |= BGZF_ERR_IO;
            return -1;
        }
    }

    // when we get here, the buffer is full or there is an EOF after a complete gzip member
    return BGZF_MAX_BLOCK_SIZE - fp->gz_stream->avail_out;
}

// Returns: 0 on success (BGZF header); -1 on non-BGZF GZIP header; -2 on error
static int check_header(const uint8_t *header)
{
    if ( header[0] != 31 || header[1] != 139 || header[2] != 8 ) return -2;
    return ((header[3] & 4) != 0
            && unpackInt16((uint8_t*)&header[10]) == 6
            && header[12] == 'B' && header[13] == 'C'
            && unpackInt16((uint8_t*)&header[14]) == 2) ? 0 : -1;
}

#ifdef BGZF_CACHE
static void free_cache(BGZF *fp)
{
    khint_t k;
    if (fp->is_write) return;
    khash_t(cache) *h = fp->cache->h;
    for (k = kh_begin(h); k < kh_end(h); ++k)
        if (kh_exist(h, k)) free(kh_val(h, k).block);
    kh_destroy(cache, h);
    free(fp->cache);
}

static int load_block_from_cache(BGZF *fp, int64_t block_address)
{
    khint_t k;
    cache_t *p;

    khash_t(cache) *h = fp->cache->h;
    k = kh_get(cache, h, block_address);
    if (k == kh_end(h)) return 0;
    p = &kh_val(h, k);
    if (fp->block_length != 0) fp->block_offset = 0;
    fp->block_address = block_address;
    fp->block_length = p->size;
    memcpy(fp->uncompressed_block, p->block, p->size);
    if ( hseek(fp->fp, p->end_offset, SEEK_SET) < 0 )
    {
        // todo: move the error up
        hts_log_error("Could not hseek to %" PRId64, p->end_offset);
        // commented to comply with cran's warning
        // exit(1);
    }
    return p->size;
}

static void cache_block(BGZF *fp, int size)
{
    int ret;
    khint_t k, k_orig;
    uint8_t *block = NULL;
    cache_t *p;
    //fprintf(stderr, "Cache block at %llx\n", (int)fp->block_address);
    khash_t(cache) *h = fp->cache->h;
    if (BGZF_MAX_BLOCK_SIZE >= fp->cache_size) return;
    if (fp->block_length < 0 || fp->block_length > BGZF_MAX_BLOCK_SIZE) return;
    if ((kh_size(h) + 1) * BGZF_MAX_BLOCK_SIZE > (uint32_t)fp->cache_size) {
        /* Remove uniformly from any position in the hash by a simple
         * round-robin approach.  An alternative strategy would be to
         * remove the least recently accessed block, but the round-robin
         * removal is simpler and is not expected to have a big impact
         * on performance */
        if (fp->cache->last_pos >= kh_end(h)) fp->cache->last_pos = kh_begin(h);
        k_orig = k = fp->cache->last_pos;
        if (++k >= kh_end(h)) k = kh_begin(h);
        while (k != k_orig) {
            if (kh_exist(h, k))
                break;
            if (++k == kh_end(h))
                k = kh_begin(h);
        }
        fp->cache->last_pos = k;

        if (k != k_orig) {
            block = kh_val(h, k).block;
            kh_del(cache, h, k);
        }
    } else {
        block = (uint8_t*)malloc(BGZF_MAX_BLOCK_SIZE);
    }
    if (!block) return;
    k = kh_put(cache, h, fp->block_address, &ret);
    if (ret <= 0) { // kh_put failed, or in there already (shouldn't happen)
        free(block);
        return;
    }
    p = &kh_val(h, k);
    p->size = fp->block_length;
    p->end_offset = fp->block_address + size;
    p->block = block;
    memcpy(p->block, fp->uncompressed_block, p->size);
}
#else
static void free_cache(BGZF *fp) {}
static int load_block_from_cache(BGZF *fp, int64_t block_address) {return 0;}
static void cache_block(BGZF *fp, int size) {}
#endif

/*
 * Absolute htell in this compressed file.
 *
 * Do not confuse with the external bgzf_tell macro which returns the virtual
 * offset.
 */
static off_t bgzf_htell(BGZF *fp) {
    if (fp->mt) {
        pthread_mutex_lock(&fp->mt->job_pool_m);
        off_t pos = fp->block_address + fp->block_clength;
        pthread_mutex_unlock(&fp->mt->job_pool_m);
        return pos;
    } else {
        return htell(fp->fp);
    }
}

int bgzf_read_block(BGZF *fp)
{
    hts_tpool_result *r;

    if (fp->errcode) return -1;

    if (fp->mt) {
    again:
        if (fp->mt->hit_eof) {
            // Further reading at EOF will always return 0
            fp->block_length = 0;
            return 0;
        }
        r = hts_tpool_next_result_wait(fp->mt->out_queue);
        bgzf_job *j = r ? (bgzf_job *)hts_tpool_result_data(r) : NULL;

        if (!j || j->errcode == BGZF_ERR_MT) {
            if (!fp->mt->free_block) {
                fp->uncompressed_block = malloc(2 * BGZF_MAX_BLOCK_SIZE);
                if (fp->uncompressed_block == NULL) return -1;
                fp->compressed_block = (char *)fp->uncompressed_block + BGZF_MAX_BLOCK_SIZE;
            } // else it's already allocated with malloc, maybe even in-use.
            if (mt_destroy(fp->mt) < 0) {
                fp->errcode = BGZF_ERR_IO;
            }
            fp->mt = NULL;
            hts_tpool_delete_result(r, 0);
            if (fp->errcode) {
                return -1;
            }
            goto single_threaded;
        }

        if (j->errcode) {
            fp->errcode = j->errcode;
            hts_log_error("BGZF decode jobs returned error %d "
                          "for block offset %"PRId64,
                          j->errcode, j->block_address);
            hts_tpool_delete_result(r, 0);
            return -1;
        }

        if (j->hit_eof) {
            if (!fp->last_block_eof && !fp->no_eof_block) {
                fp->no_eof_block = 1;
                hts_log_warning("EOF marker is absent. The input may be truncated");
            }
            fp->mt->hit_eof = 1;
        }

        // Zero length blocks in the middle of a file are (wrongly)
        // considered as EOF by many callers.  We work around this by
        // trying again to see if we hit a genuine EOF.
        if (!j->hit_eof && j->uncomp_len == 0) {
            fp->last_block_eof = 1;
            hts_tpool_delete_result(r, 0);
            goto again;
        }

        // block_length=0 and block_offset set by bgzf_seek.
        if (fp->block_length != 0) fp->block_offset = 0;
        if (!j->hit_eof) fp->block_address = j->block_address;
        fp->block_clength = j->comp_len;
        fp->block_length = j->uncomp_len;
        // bgzf_read() can change fp->block_length
        fp->last_block_eof = (fp->block_length == 0);

        if ( j->uncomp_len && j->fp->idx_build_otf )
        {
            bgzf_index_add_block(j->fp);
            j->fp->idx->ublock_addr += j->uncomp_len;
        }

        // Steal the data block as it's quicker than a memcpy.
        // We just need to make sure we delay the pool free.
        if (fp->mt->curr_job) {
            pthread_mutex_lock(&fp->mt->job_pool_m);
            pool_free(fp->mt->job_pool, fp->mt->curr_job);
            pthread_mutex_unlock(&fp->mt->job_pool_m);
        }
        fp->uncompressed_block = j->uncomp_data;
        fp->mt->curr_job = j;
        if (fp->mt->free_block) {
            free(fp->mt->free_block); // clear up last non-mt block
            fp->mt->free_block = NULL;
        }

        hts_tpool_delete_result(r, 0);
        return 0;
    }

    uint8_t header[BLOCK_HEADER_LENGTH], *compressed_block;
    int count, size, block_length, remaining;

 single_threaded:
    size = 0;

    int64_t block_address;
    block_address = bgzf_htell(fp);

    // Reading an uncompressed file
    if ( !fp->is_compressed )
    {
        count = hread(fp->fp, fp->uncompressed_block, BGZF_MAX_BLOCK_SIZE);
        if (count < 0)  // Error
        {
            hts_log_error("Failed to read uncompressed data "
                          "at offset %"PRId64"%s%s",
                          block_address, errno ? ": " : "", strerror(errno));
            fp->errcode |= BGZF_ERR_IO;
            return -1;
        }
        else if (count == 0)  // EOF
        {
            fp->block_length = 0;
            return 0;
        }
        if (fp->block_length != 0) fp->block_offset = 0;
        fp->block_address = block_address;
        fp->block_length = count;
        return 0;
    }

    // Reading compressed file
    if ( fp->is_gzip && fp->gz_stream ) // is this is an initialized gzip stream?
    {
        count = inflate_gzip_block(fp);
        if ( count<0 )
        {
            hts_log_error("Reading GZIP stream failed at offset %"PRId64,
                          block_address);
            fp->errcode |= BGZF_ERR_ZLIB;
            return -1;
        }
        fp->block_length = count;
        fp->block_address = block_address;
        return 0;
    }
    if (fp->cache_size && load_block_from_cache(fp, block_address)) return 0;

    // loop to skip empty bgzf blocks
    while (1)
    {
        count = hread(fp->fp, header, sizeof(header));
        if (count == 0) { // no data read
            if (!fp->last_block_eof && !fp->no_eof_block && !fp->is_gzip) {
                fp->no_eof_block = 1;
                hts_log_warning("EOF marker is absent. The input may be truncated");
            }
            fp->block_length = 0;
            return 0;
        }
        int ret = 0;
        if ( count != sizeof(header) || (ret=check_header(header))==-2 )
        {
            fp->errcode |= BGZF_ERR_HEADER;
            hts_log_error("%s BGZF header at offset %"PRId64,
                          ret ? "Invalid" : "Failed to read",
                          block_address);
            return -1;
        }
        if ( ret==-1 )
        {
            // GZIP, not BGZF
            uint8_t *cblock = (uint8_t*)fp->compressed_block;
            memcpy(cblock, header, sizeof(header));
            count = hread(fp->fp, cblock+sizeof(header), BGZF_BLOCK_SIZE - sizeof(header)) + sizeof(header);

            fp->is_gzip = 1;
            fp->gz_stream = (z_stream*) calloc(1,sizeof(z_stream));
            // Set up zlib, using a window size of 15, and its built-in GZIP header processing (+16).
            int ret = inflateInit2(fp->gz_stream, 15 + 16);
            if (ret != Z_OK)
            {
                hts_log_error("Call to inflateInit2 failed: %s", bgzf_zerr(ret, fp->gz_stream));
                fp->errcode |= BGZF_ERR_ZLIB;
                return -1;
            }
            fp->gz_stream->avail_in = count;
            fp->gz_stream->next_in  = cblock;
            count = inflate_gzip_block(fp);
            if ( count<0 )
            {
                hts_log_error("Reading GZIP stream failed at offset %"PRId64,
                              block_address);
                fp->errcode |= BGZF_ERR_ZLIB;
                return -1;
            }
            fp->block_length = count;
            fp->block_address = block_address;
            if ( fp->idx_build_otf ) return -1; // cannot build index for gzip
            return 0;
        }
        size = count;
        block_length = unpackInt16((uint8_t*)&header[16]) + 1; // +1 because when writing this number, we used "-1"
        if (block_length < BLOCK_HEADER_LENGTH)
        {
            hts_log_error("Invalid BGZF block length at offset %"PRId64,
                          block_address);
            fp->errcode |= BGZF_ERR_HEADER;
            return -1;
        }
        compressed_block = (uint8_t*)fp->compressed_block;
        memcpy(compressed_block, header, BLOCK_HEADER_LENGTH);
        remaining = block_length - BLOCK_HEADER_LENGTH;
        count = hread(fp->fp, &compressed_block[BLOCK_HEADER_LENGTH], remaining);
        if (count != remaining) {
            hts_log_error("Failed to read BGZF block data at offset %"PRId64
                          " expected %d bytes; hread returned %d",
                          block_address, remaining, count);
            fp->errcode |= BGZF_ERR_IO;
            return -1;
        }
        size += count;
        if ((count = inflate_block(fp, block_length)) < 0) {
            hts_log_debug("Inflate block operation failed for "
                          "block at offset %"PRId64": %s",
                          block_address, bgzf_zerr(count, NULL));
            fp->errcode |= BGZF_ERR_ZLIB;
            return -1;
        }
        fp->last_block_eof = (count == 0);
        if ( count ) break;     // otherwise an empty bgzf block
        block_address = bgzf_htell(fp); // update for new block start
    }
    if (fp->block_length != 0) fp->block_offset = 0; // Do not reset offset if this read follows a seek.
    fp->block_address = block_address;
    fp->block_length = count;
    if ( fp->idx_build_otf )
    {
        bgzf_index_add_block(fp);
        fp->idx->ublock_addr += count;
    }
    cache_block(fp, size);
    return 0;
}

ssize_t bgzf_read(BGZF *fp, void *data, size_t length)
{
    ssize_t bytes_read = 0;
    uint8_t *output = (uint8_t*)data;
    if (length <= 0) return 0;
    assert(fp->is_write == 0);
    while (bytes_read < length) {
        int copy_length, available = fp->block_length - fp->block_offset;
        uint8_t *buffer;
        if (available <= 0) {
            int ret = bgzf_read_block(fp);
            if (ret != 0) {
                hts_log_error("Read block operation failed with error %d after %zd of %zu bytes", fp->errcode, bytes_read, length);
                fp->errcode |= BGZF_ERR_ZLIB;
                return -1;
            }
            available = fp->block_length - fp->block_offset;
            if (available == 0) {
                if (fp->block_length == 0)
                    break; // EOF

                // Offset was at end of block (see commit e9863a0)
                fp->block_address = bgzf_htell(fp);
                fp->block_offset = fp->block_length = 0;
                continue;
            } else if (available < 0) {
                // Block offset was set to an invalid coordinate
                hts_log_error("BGZF block offset %d set beyond block size %d",
                              fp->block_offset, fp->block_length);
                fp->errcode |= BGZF_ERR_MISUSE;
                return -1;
            }
        }
        copy_length = length - bytes_read < available? length - bytes_read : available;
        buffer = (uint8_t*)fp->uncompressed_block;
        memcpy(output, buffer + fp->block_offset, copy_length);
        fp->block_offset += copy_length;
        output += copy_length;
        bytes_read += copy_length;

        // For raw gzip streams this avoids short reads.
        if (fp->block_offset == fp->block_length) {
            fp->block_address = bgzf_htell(fp);
            fp->block_offset = fp->block_length = 0;
        }
    }

    fp->uncompressed_address += bytes_read;

    return bytes_read;
}

// -1 for EOF, -2 for error, 0-255 for byte.
int bgzf_peek(BGZF *fp) {
    int available = fp->block_length - fp->block_offset;
    if (available <= 0) {
        if (bgzf_read_block(fp) < 0) {
            hts_log_error("Read block operation failed with error %d", fp->errcode);
            fp->errcode = BGZF_ERR_ZLIB;
            return -2;
        }
    }
    available = fp->block_length - fp->block_offset;
    if (available)
        return ((unsigned char *)fp->uncompressed_block)[fp->block_offset];

    return -1;
}

ssize_t bgzf_raw_read(BGZF *fp, void *data, size_t length)
{
    ssize_t ret = hread(fp->fp, data, length);
    if (ret < 0) fp->errcode |= BGZF_ERR_IO;
    return ret;
}

#ifdef BGZF_MT

/* Function to clean up when jobs are discarded (e.g. during seek)
 * This works for results too, as results are the same struct with
 * decompressed data stored in it. */
static void job_cleanup(void *arg) {
    bgzf_job *j = (bgzf_job *)arg;
    mtaux_t *mt = j->fp->mt;
    pthread_mutex_lock(&mt->job_pool_m);
    pool_free(mt->job_pool, j);
    pthread_mutex_unlock(&mt->job_pool_m);
}

static void *bgzf_encode_func(void *arg) {
    bgzf_job *j = (bgzf_job *)arg;

    j->comp_len = BGZF_MAX_BLOCK_SIZE;
    int ret = bgzf_compress(j->comp_data, &j->comp_len,
                            j->uncomp_data, j->uncomp_len,
                            j->fp->compress_level);
    if (ret != 0)
        j->errcode |= BGZF_ERR_ZLIB;

    return arg;
}

// Optimisation for compression level 0 (uncompressed deflate blocks)
// Avoids memcpy of the data from uncompressed to compressed buffer.
static void *bgzf_encode_level0_func(void *arg) {
    bgzf_job *j = (bgzf_job *)arg;
    uint32_t crc;
    j->comp_len = j->uncomp_len + BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH + 5;

    // Data will have already been copied in to
    // j->comp_data + BLOCK_HEADER_LENGTH + 5

    // Add preamble
    memcpy(j->comp_data, g_magic, BLOCK_HEADER_LENGTH);
    u16_to_le(j->comp_len-1, j->comp_data + 16);

    // Deflate uncompressed data header
    j->comp_data[BLOCK_HEADER_LENGTH] = 1; // BFINAL=1, BTYPE=00; see RFC1951
    u16_to_le(j->uncomp_len, j->comp_data + BLOCK_HEADER_LENGTH + 1);
    u16_to_le(~j->uncomp_len, j->comp_data + BLOCK_HEADER_LENGTH + 3);

    // Trailer (CRC, uncompressed length)
    crc = crc32(crc32(0L, NULL, 0L),
                (Bytef*)j->comp_data + BLOCK_HEADER_LENGTH + 5, j->uncomp_len);
    u32_to_le(crc, j->comp_data +  j->comp_len - 8);
    u32_to_le(j->uncomp_len, j->comp_data + j->comp_len - 4);

    return arg;
}

// Our input block has already been decoded by bgzf_mt_read_block().
// We need to split that into a fetch block (compressed) and make this
// do the actual decompression step.
static void *bgzf_decode_func(void *arg) {
    bgzf_job *j = (bgzf_job *)arg;

    j->uncomp_len = BGZF_MAX_BLOCK_SIZE;
    uint32_t crc = le_to_u32((uint8_t *)j->comp_data + j->comp_len-8);
    int ret = bgzf_uncompress(j->uncomp_data, &j->uncomp_len,
                              j->comp_data+18, j->comp_len-18, crc);
    if (ret != 0)
        j->errcode |= BGZF_ERR_ZLIB;

    return arg;
}

/*
 * Nul function so we can dispatch a job with the correct serial
 * to mark failure or to indicate an empty read (EOF).
 */
static void *bgzf_nul_func(void *arg) { return arg; }

/*
 * Takes compressed blocks off the results queue and calls hwrite to
 * punt them to the output stream.
 *
 * Returns NULL when no more are left, or -1 on error
 */
static void *bgzf_mt_writer(void *vp) {
    BGZF *fp = (BGZF *)vp;
    mtaux_t *mt = fp->mt;
    hts_tpool_result *r;

    if (fp->idx_build_otf) {
        fp->idx->moffs = fp->idx->noffs = 1;
        fp->idx->offs = (bgzidx1_t*) calloc(fp->idx->moffs, sizeof(bgzidx1_t));
        if (!fp->idx->offs) goto err;
    }

    // Iterates until result queue is shutdown, where it returns NULL.
    while ((r = hts_tpool_next_result_wait(mt->out_queue))) {
        bgzf_job *j = (bgzf_job *)hts_tpool_result_data(r);
        assert(j);

        if (fp->idx_build_otf) {
            fp->idx->noffs++;
            if ( fp->idx->noffs > fp->idx->moffs )
            {
                fp->idx->moffs = fp->idx->noffs;
                kroundup32(fp->idx->moffs);
                fp->idx->offs = (bgzidx1_t*) realloc(fp->idx->offs, fp->idx->moffs*sizeof(bgzidx1_t));
                if ( !fp->idx->offs ) goto err;
            }
            fp->idx->offs[ fp->idx->noffs-1 ].uaddr = fp->idx->offs[ fp->idx->noffs-2 ].uaddr + j->uncomp_len;
            fp->idx->offs[ fp->idx->noffs-1 ].caddr = fp->idx->offs[ fp->idx->noffs-2 ].caddr + j->comp_len;
        }

        // Flush any cached hts_idx_push calls
        if (bgzf_idx_flush(fp) < 0)
            goto err;

        if (hwrite(fp->fp, j->comp_data, j->comp_len) != j->comp_len)
            goto err;

        // Update our local block_address.  Cannot be fp->block_address due to no
        // locking in bgzf_tell.
        pthread_mutex_lock(&mt->idx_m);
        mt->block_address += j->comp_len;
        pthread_mutex_unlock(&mt->idx_m);

        /*
         * Periodically call hflush (which calls fsync when on a file).
         * This avoids the fsync being done at the bgzf_close stage,
         * which can sometimes cause significant delays.  As this is in
         * a separate thread, spreading the sync delays throughout the
         * program execution seems better.
         * Frequency of 1/512 has been chosen by experimentation
         * across local XFS, NFS and Lustre tests.
         */
        if (++mt->flush_pending % 512 == 0)
            if (hflush(fp->fp) != 0)
                goto err;


        hts_tpool_delete_result(r, 0);

        // Also updated by main thread
        pthread_mutex_lock(&mt->job_pool_m);
        pool_free(mt->job_pool, j);
        mt->jobs_pending--;
        pthread_mutex_unlock(&mt->job_pool_m);
    }

    if (hflush(fp->fp) != 0)
        goto err;

    hts_tpool_process_destroy(mt->out_queue);

    return NULL;

 err:
    hts_tpool_process_destroy(mt->out_queue);
    return (void *)-1;
}


/*
 * Reads a compressed block of data using hread and dispatches it to
 * the thread pool for decompression.  This is the analogue of the old
 * non-threaded bgzf_read_block() function, but without modifying fp
 * in any way (except for the read offset).  All output goes via the
 * supplied bgzf_job struct.
 *
 * Returns NULL when no more are left, or -1 on error
 */
int bgzf_mt_read_block(BGZF *fp, bgzf_job *j)
{
    uint8_t header[BLOCK_HEADER_LENGTH], *compressed_block;
    int count, block_length, remaining;

    // NOTE: Guaranteed to be compressed as we block multi-threading in
    // uncompressed mode.  However it may be gzip compression instead
    // of bgzf.

    // Reading compressed file
    int64_t block_address;
    block_address = htell(fp->fp);

    j->block_address = block_address;  // in case we exit with j->errcode

    if (fp->cache_size && load_block_from_cache(fp, block_address)) return 0;
    count = hpeek(fp->fp, header, sizeof(header));
    if (count == 0) // no data read
        return -1;
    int ret;
    if ( count != sizeof(header) || (ret=check_header(header))==-2 )
    {
        j->errcode |= BGZF_ERR_HEADER;
        return -1;
    }
    if (ret == -1) {
        j->errcode |= BGZF_ERR_MT;
        return -1;
    }

    count = hread(fp->fp, header, sizeof(header));
    if (count != sizeof(header)) // no data read
        return -1;

    block_length = unpackInt16((uint8_t*)&header[16]) + 1; // +1 because when writing this number, we used "-1"
    if (block_length < BLOCK_HEADER_LENGTH) {
        j->errcode |= BGZF_ERR_HEADER;
        return -1;
    }
    compressed_block = (uint8_t*)j->comp_data;
    memcpy(compressed_block, header, BLOCK_HEADER_LENGTH);
    remaining = block_length - BLOCK_HEADER_LENGTH;
    count = hread(fp->fp, &compressed_block[BLOCK_HEADER_LENGTH], remaining);
    if (count != remaining) {
        j->errcode |= BGZF_ERR_IO;
        return -1;
    }
    j->comp_len = block_length;
    j->uncomp_len = BGZF_MAX_BLOCK_SIZE;
    j->block_address = block_address;
    j->fp = fp;
    j->errcode = 0;

    return 0;
}


static int bgzf_check_EOF_common(BGZF *fp)
{
    uint8_t buf[28];
    off_t offset = htell(fp->fp);
    if (hseek(fp->fp, -28, SEEK_END) < 0) {
        if (errno == ESPIPE) { hclearerr(fp->fp); return 2; }
#ifdef _WIN32
        if (errno == EINVAL) { hclearerr(fp->fp); return 2; }
#else
        // Assume that EINVAL was due to the file being less than 28 bytes
        // long, rather than being a random error return from an hfile backend.
        // This should be reported as "no EOF block" rather than an error.
        if (errno == EINVAL) { hclearerr(fp->fp); return 0; }
#endif
        return -1;
    }
    if ( hread(fp->fp, buf, 28) != 28 ) return -1;
    if ( hseek(fp->fp, offset, SEEK_SET) < 0 ) return -1;
    return (memcmp("\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\033\0\3\0\0\0\0\0\0\0\0\0", buf, 28) == 0)? 1 : 0;
}

/*
 * Checks EOF from the reader thread.
 */
static void bgzf_mt_eof(BGZF *fp) {
    mtaux_t *mt = fp->mt;

    pthread_mutex_lock(&mt->job_pool_m);
    mt->eof = bgzf_check_EOF_common(fp);
    pthread_mutex_unlock(&mt->job_pool_m);
    mt->command = HAS_EOF_DONE;
    pthread_cond_signal(&mt->command_c);
}


/*
 * Performs the seek (called by reader thread).
 *
 * This simply drains the entire queue, throwing away blocks, seeks,
 * and starts it up again.  Brute force, but maybe sufficient.
 */
static void bgzf_mt_seek(BGZF *fp) {
    mtaux_t *mt = fp->mt;

    hts_tpool_process_reset(mt->out_queue, 0);
    pthread_mutex_lock(&mt->job_pool_m);
    mt->errcode = 0;

    if (hseek(fp->fp, mt->block_address, SEEK_SET) < 0)
        mt->errcode = BGZF_ERR_IO;

    pthread_mutex_unlock(&mt->job_pool_m);
    mt->command = SEEK_DONE;
    pthread_cond_signal(&mt->command_c);
}

static void *bgzf_mt_reader(void *vp) {
    BGZF *fp = (BGZF *)vp;
    mtaux_t *mt = fp->mt;

restart:
    pthread_mutex_lock(&mt->job_pool_m);
    bgzf_job *j = pool_alloc(mt->job_pool);
    pthread_mutex_unlock(&mt->job_pool_m);
    if (!j) goto err;
    j->errcode = 0;
    j->comp_len = 0;
    j->uncomp_len = 0;
    j->hit_eof = 0;
    j->fp = fp;

    while (bgzf_mt_read_block(fp, j) == 0) {
        // Dispatch
        if (hts_tpool_dispatch3(mt->pool, mt->out_queue, bgzf_decode_func, j,
                                job_cleanup, job_cleanup, 0) < 0) {
            job_cleanup(j);
            goto err;
        }

        // Check for command
        pthread_mutex_lock(&mt->command_m);
        switch (mt->command) {
        case SEEK:
            bgzf_mt_seek(fp);  // Sets mt->command to SEEK_DONE
            pthread_mutex_unlock(&mt->command_m);
            goto restart;

        case HAS_EOF:
            bgzf_mt_eof(fp);   // Sets mt->command to HAS_EOF_DONE
            break;

        case SEEK_DONE:
        case HAS_EOF_DONE:
            pthread_cond_signal(&mt->command_c);
            break;

        case CLOSE:
            pthread_cond_signal(&mt->command_c);
            pthread_mutex_unlock(&mt->command_m);
            hts_tpool_process_destroy(mt->out_queue);
            return NULL;

        default:
            break;
        }
        pthread_mutex_unlock(&mt->command_m);

        // Allocate buffer for next block
        pthread_mutex_lock(&mt->job_pool_m);
        j = pool_alloc(mt->job_pool);
        pthread_mutex_unlock(&mt->job_pool_m);
        if (!j) {
            hts_tpool_process_destroy(mt->out_queue);
            return NULL;
        }
        j->errcode = 0;
        j->comp_len = 0;
        j->uncomp_len = 0;
        j->hit_eof = 0;
        j->fp = fp;
    }

    if (j->errcode == BGZF_ERR_MT) {
        // Attempt to multi-thread decode a raw gzip stream cannot be done.
        // We tear down the multi-threaded decoder and revert to the old code.
        if (hts_tpool_dispatch3(mt->pool, mt->out_queue, bgzf_nul_func, j,
                                job_cleanup, job_cleanup, 0) < 0) {
            job_cleanup(j);
            hts_tpool_process_destroy(mt->out_queue);
            return NULL;
        }
        hts_tpool_process_ref_decr(mt->out_queue);
        return &j->errcode;
    }

    // Dispatch an empty block so EOF is spotted.
    // We also use this mechanism for returning errors, in which case
    // j->errcode is set already.

    j->hit_eof = 1;
    if (hts_tpool_dispatch3(mt->pool, mt->out_queue, bgzf_nul_func, j,
                            job_cleanup, job_cleanup, 0) < 0) {
        job_cleanup(j);
        hts_tpool_process_destroy(mt->out_queue);
        return NULL;
    }
    if (j->errcode != 0) {
        hts_tpool_process_destroy(mt->out_queue);
        return &j->errcode;
    }

    // We hit EOF so can stop reading, but we may get a subsequent
    // seek request.  In this case we need to restart the reader.
    //
    // To handle this we wait on a condition variable and then
    // monitor the command. (This could be either seek or close.)
    for (;;) {
        pthread_mutex_lock(&mt->command_m);
        if (mt->command == NONE)
            pthread_cond_wait(&mt->command_c, &mt->command_m);
        switch(mt->command) {
        default:
            pthread_mutex_unlock(&mt->command_m);
            break;

        case SEEK:
            bgzf_mt_seek(fp);
            pthread_mutex_unlock(&mt->command_m);
            goto restart;

        case HAS_EOF:
            bgzf_mt_eof(fp);   // Sets mt->command to HAS_EOF_DONE
            pthread_mutex_unlock(&mt->command_m);
            break;

        case SEEK_DONE:
        case HAS_EOF_DONE:
            pthread_cond_signal(&mt->command_c);
            pthread_mutex_unlock(&mt->command_m);
            break;

        case CLOSE:
            pthread_cond_signal(&mt->command_c);
            pthread_mutex_unlock(&mt->command_m);
            hts_tpool_process_destroy(mt->out_queue);
            return NULL;
        }
    }

 err:
    pthread_mutex_lock(&mt->command_m);
    mt->command = CLOSE;
    pthread_cond_signal(&mt->command_c);
    pthread_mutex_unlock(&mt->command_m);
    hts_tpool_process_destroy(mt->out_queue);
    return NULL;
}

int bgzf_thread_pool(BGZF *fp, hts_tpool *pool, int qsize) {
    // No gain from multi-threading when not compressed
    if (!fp->is_compressed)
        return 0;

    mtaux_t *mt;
    mt = (mtaux_t*)calloc(1, sizeof(mtaux_t));
    if (!mt) return -1;
    fp->mt = mt;

    mt->pool = pool;
    mt->n_threads = hts_tpool_size(pool);
    if (!qsize)
        qsize = mt->n_threads*2;
    if (!(mt->out_queue = hts_tpool_process_init(mt->pool, qsize, 0)))
        goto err;
    hts_tpool_process_ref_incr(mt->out_queue);

    mt->job_pool = pool_create(sizeof(bgzf_job));
    if (!mt->job_pool)
        goto err;

    pthread_mutex_init(&mt->job_pool_m, NULL);
    pthread_mutex_init(&mt->command_m, NULL);
    pthread_mutex_init(&mt->idx_m, NULL);
    pthread_cond_init(&mt->command_c, NULL);
    mt->flush_pending = 0;
    mt->jobs_pending = 0;
    mt->free_block = fp->uncompressed_block; // currently in-use block
    mt->block_address = fp->block_address;
    pthread_create(&mt->io_task, NULL,
                   fp->is_write ? bgzf_mt_writer : bgzf_mt_reader, fp);

    return 0;

 err:
    free(mt);
    fp->mt = NULL;
    return -1;
}

int bgzf_mt(BGZF *fp, int n_threads, int n_sub_blks)
{
    // No gain from multi-threading when not compressed
    if (!fp->is_compressed || fp->is_gzip)
        return 0;

    if (n_threads < 1) return -1;
    hts_tpool *p = hts_tpool_init(n_threads);
    if (!p)
        return -1;

    if (bgzf_thread_pool(fp, p, 0) != 0) {
        hts_tpool_destroy(p);
        return -1;
    }

    fp->mt->own_pool = 1;

    return 0;
}

static int mt_destroy(mtaux_t *mt)
{
    int ret = 0;

    // Tell the reader to shut down
    pthread_mutex_lock(&mt->command_m);
    mt->command = CLOSE;
    pthread_cond_signal(&mt->command_c);
    hts_tpool_wake_dispatch(mt->out_queue); // unstick the reader
    pthread_mutex_unlock(&mt->command_m);

    // Check for thread worker failure, indicated by is_shutdown returning 2
    // It's possible really late errors might be missed, but we can live with
    // that.
    ret = -(hts_tpool_process_is_shutdown(mt->out_queue) > 1);
    // Destroying the queue first forces the writer to exit.
    // mt->out_queue is reference counted, so destroy gets called in both
    // this and the IO threads.  The last to do it will clean up.
    hts_tpool_process_destroy(mt->out_queue);

    // IO thread will now exit.  Wait for it and perform final clean-up.
    // If it returned non-NULL, it was not happy.
    void *retval = NULL;
    pthread_join(mt->io_task, &retval);
    ret = retval != NULL ? -1 : ret;

    pthread_mutex_destroy(&mt->job_pool_m);
    pthread_mutex_destroy(&mt->command_m);
    pthread_mutex_destroy(&mt->idx_m);
    pthread_cond_destroy(&mt->command_c);
    if (mt->curr_job)
        pool_free(mt->job_pool, mt->curr_job);

    if (mt->own_pool)
        hts_tpool_destroy(mt->pool);

    pool_destroy(mt->job_pool);

    if (mt->idx_cache.e)
        free(mt->idx_cache.e);

    free(mt);
    //fflush(stderr);

    return ret;
}

static int mt_queue(BGZF *fp)
{
    mtaux_t *mt = fp->mt;

    mt->block_number++;

    // Also updated by writer thread
    pthread_mutex_lock(&mt->job_pool_m);
    bgzf_job *j = pool_alloc(mt->job_pool);
    if (j) mt->jobs_pending++;
    pthread_mutex_unlock(&mt->job_pool_m);
    if (!j) return -1;

    j->fp = fp;
    j->errcode = 0;
    j->uncomp_len  = fp->block_offset;
    if (fp->compress_level == 0) {
        memcpy(j->comp_data + BLOCK_HEADER_LENGTH + 5, fp->uncompressed_block,
               j->uncomp_len);
        if (hts_tpool_dispatch3(mt->pool, mt->out_queue,
                                bgzf_encode_level0_func, j,
                                job_cleanup, job_cleanup, 0) < 0) {
            goto fail;
        }
    } else {
        memcpy(j->uncomp_data, fp->uncompressed_block, j->uncomp_len);

        // Need non-block vers & job_pending?
        if (hts_tpool_dispatch3(mt->pool, mt->out_queue, bgzf_encode_func, j,
                                job_cleanup, job_cleanup, 0) < 0) {
            goto fail;
        }
    }

    fp->block_offset = 0;
    return 0;

 fail:
    job_cleanup(j);
    pthread_mutex_lock(&mt->job_pool_m);
    mt->jobs_pending--;
    pthread_mutex_unlock(&mt->job_pool_m);
    return -1;
}

static int mt_flush_queue(BGZF *fp)
{
    mtaux_t *mt = fp->mt;

    // Drain the encoder jobs.
    // We cannot use hts_tpool_flush here as it can cause deadlock if
    // the queue is full up of decoder tasks.  The best solution would
    // be to have one input queue per type of job, but we don't right now.
    //hts_tpool_flush(mt->pool);
    pthread_mutex_lock(&mt->job_pool_m);
    int shutdown = 0;
    while (mt->jobs_pending != 0) {
        if ((shutdown = hts_tpool_process_is_shutdown(mt->out_queue)))
            break;
        pthread_mutex_unlock(&mt->job_pool_m);
        usleep(10000); // FIXME: replace by condition variable
        pthread_mutex_lock(&mt->job_pool_m);
    }
    pthread_mutex_unlock(&mt->job_pool_m);

    if (shutdown)
        return -1;

    // Wait on bgzf_mt_writer to drain the queue
    if (hts_tpool_process_flush(mt->out_queue) != 0)
        return -1;

    return (fp->errcode == 0)? 0 : -1;
}

static int lazy_flush(BGZF *fp)
{
    if (fp->mt)
        return fp->block_offset ? mt_queue(fp) : 0;
    else
        return bgzf_flush(fp);
}

#else  // ~ #ifdef BGZF_MT

int bgzf_mt(BGZF *fp, int n_threads, int n_sub_blks)
{
    return 0;
}

static inline int lazy_flush(BGZF *fp)
{
    return bgzf_flush(fp);
}

#endif // ~ #ifdef BGZF_MT

int bgzf_flush(BGZF *fp)
{
    if (!fp->is_write) return 0;
#ifdef BGZF_MT
    if (fp->mt) {
        int ret = 0;
        if (fp->block_offset) ret = mt_queue(fp);
        if (!ret) ret = mt_flush_queue(fp);

        // We maintain mt->block_address when threading as the
        // main code can call bgzf_tell without any locks.
        // (The result from tell are wrong, but we only care about the last
        // 16-bits worth except for the final flush process.
        pthread_mutex_lock(&fp->mt->idx_m);
        fp->block_address = fp->mt->block_address;
        pthread_mutex_unlock(&fp->mt->idx_m);

        return ret;
    }
#endif
    while (fp->block_offset > 0) {
        int block_length;
        if ( fp->idx_build_otf )
        {
            bgzf_index_add_block(fp);
            fp->idx->ublock_addr += fp->block_offset;
        }
        block_length = deflate_block(fp, fp->block_offset);
        if (block_length < 0) {
            hts_log_debug("Deflate block operation failed: %s", bgzf_zerr(block_length, NULL));
            return -1;
        }
        if (hwrite(fp->fp, fp->compressed_block, block_length) != block_length) {
            hts_log_error("File write failed (wrong size)");
            fp->errcode |= BGZF_ERR_IO; // possibly truncated file
            return -1;
        }
        fp->block_address += block_length;
    }
    return 0;
}

int bgzf_flush_try(BGZF *fp, ssize_t size)
{
    if (fp->block_offset + size > BGZF_BLOCK_SIZE) return lazy_flush(fp);
    return 0;
}

ssize_t bgzf_write(BGZF *fp, const void *data, size_t length)
{
    if ( !fp->is_compressed ) {
        size_t push = length + (size_t) fp->block_offset;
        fp->block_offset = push % BGZF_MAX_BLOCK_SIZE;
        fp->block_address += (push - fp->block_offset);
        return hwrite(fp->fp, data, length);
    }

    const uint8_t *input = (const uint8_t*)data;
    ssize_t remaining = length;
    assert(fp->is_write);
    while (remaining > 0) {
        uint8_t* buffer = (uint8_t*)fp->uncompressed_block;
        int copy_length = BGZF_BLOCK_SIZE - fp->block_offset;
        if (copy_length > remaining) copy_length = remaining;
        memcpy(buffer + fp->block_offset, input, copy_length);
        fp->block_offset += copy_length;
        input += copy_length;
        remaining -= copy_length;
        if (fp->block_offset == BGZF_BLOCK_SIZE) {
            if (lazy_flush(fp) != 0) return -1;
        }
    }
    return length - remaining;
}

// Helper function for tidying up fp->mt and setting errcode
static void bgzf_close_mt(BGZF *fp) {
    if (fp->mt) {
        if (!fp->mt->free_block)
            fp->uncompressed_block = NULL;
        if (mt_destroy(fp->mt) < 0)
            fp->errcode = BGZF_ERR_IO;
    }
}

HTSLIB_EXPORT
int bgzf_close(BGZF* fp)
{
    int ret, block_length;
    if (fp == 0) return -1;
    if (fp->is_write && fp->is_compressed) {
        if (bgzf_flush(fp) != 0) {
            bgzf_close_mt(fp);
            return -1;
        }
        fp->compress_level = -1;
        block_length = deflate_block(fp, 0); // write an empty block
        if (block_length < 0) {
            hts_log_debug("Deflate block operation failed: %s", bgzf_zerr(block_length, NULL));
            bgzf_close_mt(fp);
            return -1;
        }
        if (hwrite(fp->fp, fp->compressed_block, block_length) < 0
            || hflush(fp->fp) != 0) {
            hts_log_error("File write failed");
            fp->errcode |= BGZF_ERR_IO;
            return -1;
        }
    }

    bgzf_close_mt(fp);

    if ( fp->is_gzip )
    {
        if (fp->gz_stream == NULL) ret = Z_OK;
        else if (!fp->is_write) ret = inflateEnd(fp->gz_stream);
        else ret = deflateEnd(fp->gz_stream);
        if (ret != Z_OK) {
            hts_log_error("Call to inflateEnd/deflateEnd failed: %s", bgzf_zerr(ret, NULL));
        }
        free(fp->gz_stream);
    }
    ret = hclose(fp->fp);
    if (ret != 0) return -1;
    bgzf_index_destroy(fp);
    free(fp->uncompressed_block);
    free_cache(fp);
    ret = fp->errcode ? -1 : 0;
    free(fp);
    return ret;
}

void bgzf_set_cache_size(BGZF *fp, int cache_size)
{
    if (fp && fp->mt) return; // Not appropriate when multi-threading
    if (fp && fp->cache) fp->cache_size = cache_size;
}

int bgzf_check_EOF(BGZF *fp) {
    int has_eof;

    if (fp->mt) {
        pthread_mutex_lock(&fp->mt->command_m);
        // fp->mt->command state transitions should be:
        // NONE -> HAS_EOF -> HAS_EOF_DONE -> NONE
        // (HAS_EOF -> HAS_EOF_DONE happens in bgzf_mt_reader thread)
        if (fp->mt->command != CLOSE)
            fp->mt->command = HAS_EOF;
        pthread_cond_signal(&fp->mt->command_c);
        hts_tpool_wake_dispatch(fp->mt->out_queue);
        do {
            if (fp->mt->command == CLOSE) {
                // possible error in bgzf_mt_reader
                pthread_mutex_unlock(&fp->mt->command_m);
                return 0;
            }
            pthread_cond_wait(&fp->mt->command_c, &fp->mt->command_m);
            switch (fp->mt->command) {
            case HAS_EOF_DONE: break;
            case HAS_EOF:
                // Resend signal intended for bgzf_mt_reader()
                pthread_cond_signal(&fp->mt->command_c);
                break;
            case CLOSE:
                continue;
            default: 
               // commented to comply with cran's warning
               // abort();  // Should not get to any other state
               return -1;
            }
        } while (fp->mt->command != HAS_EOF_DONE);
        fp->mt->command = NONE;
        has_eof = fp->mt->eof;
        pthread_mutex_unlock(&fp->mt->command_m);
    } else {
        has_eof = bgzf_check_EOF_common(fp);
    }

    fp->no_eof_block = (has_eof == 0);

    return has_eof;
}

static inline int64_t bgzf_seek_common(BGZF* fp,
                                       int64_t block_address, int block_offset)
{
    if (fp->mt) {
        // The reader runs asynchronous and does loops of:
        //    Read block
        //    Check & process command
        //    Dispatch decode job
        //
        // Once at EOF it then switches to loops of
        //    Wait for command
        //    Process command (possibly switching back to above loop).
        //
        // To seek we therefore send the reader thread a SEEK command,
        // waking it up if blocked in dispatch and signalling if
        // waiting for a command.  We then wait for the response so we
        // know the seek succeeded.
        pthread_mutex_lock(&fp->mt->command_m);
        fp->mt->hit_eof = 0;
        // fp->mt->command state transitions should be:
        // NONE -> SEEK -> SEEK_DONE -> NONE
        // (SEEK -> SEEK_DONE happens in bgzf_mt_reader thread)
        fp->mt->command = SEEK;
        fp->mt->block_address = block_address;
        pthread_cond_signal(&fp->mt->command_c);
        hts_tpool_wake_dispatch(fp->mt->out_queue);
        do {
            pthread_cond_wait(&fp->mt->command_c, &fp->mt->command_m);
            switch (fp->mt->command) {
            case SEEK_DONE: break;
            case SEEK:
                // Resend signal intended for bgzf_mt_reader()
                pthread_cond_signal(&fp->mt->command_c);
                break;
            default:
                // commented to comply with cran's warning
                // abort();  // Should not get to any other state
                return -1;
            }
        } while (fp->mt->command != SEEK_DONE);
        fp->mt->command = NONE;

        fp->block_length = 0;  // indicates current block has not been loaded
        fp->block_address = block_address;
        fp->block_offset = block_offset;

        pthread_mutex_unlock(&fp->mt->command_m);
    } else {
        if (hseek(fp->fp, block_address, SEEK_SET) < 0) {
            fp->errcode |= BGZF_ERR_IO;
            return -1;
        }
        fp->block_length = 0;  // indicates current block has not been loaded
        fp->block_address = block_address;
        fp->block_offset = block_offset;
    }

    return 0;
}

int64_t bgzf_seek(BGZF* fp, int64_t pos, int where)
{
    if (fp->is_write || where != SEEK_SET || fp->is_gzip) {
        fp->errcode |= BGZF_ERR_MISUSE;
        return -1;
    }

    // This is a flag to indicate we've jumped elsewhere in the stream, to act
    // as a hint to any other code which is wrapping up bgzf for its own
    // purposes.  We may not be able to tell when seek happens as it can be
    // done on our behalf, eg by the iterator.
    //
    // This is never cleared here.  Any tool that needs to handle it is also
    // responsible for clearing it.
    fp->seeked = pos;

    return bgzf_seek_common(fp, pos >> 16, pos & 0xFFFF);
}

int bgzf_is_bgzf(const char *fn)
{
    uint8_t buf[16];
    int n;
    hFILE *fp;
    if ((fp = hopen(fn, "r")) == 0) return 0;
    n = hread(fp, buf, 16);
    if (hclose(fp) < 0) return 0;
    if (n != 16) return 0;
    return check_header(buf) == 0? 1 : 0;
}

int bgzf_compression(BGZF *fp)
{
    return (!fp->is_compressed)? no_compression : (fp->is_gzip)? gzip : bgzf;
}

int bgzf_getc(BGZF *fp)
{
    if (fp->block_offset+1 < fp->block_length) {
        fp->uncompressed_address++;
        return ((unsigned char*)fp->uncompressed_block)[fp->block_offset++];
    }

    int c;
    if (fp->block_offset >= fp->block_length) {
        if (bgzf_read_block(fp) != 0) return -2; /* error */
        if (fp->block_length == 0) return -1; /* end-of-file */
    }
    c = ((unsigned char*)fp->uncompressed_block)[fp->block_offset++];
    if (fp->block_offset == fp->block_length) {
        fp->block_address = bgzf_htell(fp);
        fp->block_offset = 0;
        fp->block_length = 0;
    }
    fp->uncompressed_address++;
    return c;
}

int bgzf_getline(BGZF *fp, int delim, kstring_t *str)
{
    int l, state = 0;
    str->l = 0;
    do {
        if (fp->block_offset >= fp->block_length) {
            if (bgzf_read_block(fp) != 0) { state = -2; break; }
            if (fp->block_length == 0) { state = -1; break; }
        }
        unsigned char *buf = fp->uncompressed_block;

        // Equivalent to a naive byte by byte search from
        // buf + block_offset to buf + block_length.
        void *e = memchr(&buf[fp->block_offset], delim,
                         fp->block_length - fp->block_offset);
        l = e ? (unsigned char *)e - buf : fp->block_length;

        if (l < fp->block_length) state = 1;
        l -= fp->block_offset;
        if (ks_expand(str, l + 2) < 0) { state = -3; break; }
        memcpy(str->s + str->l, buf + fp->block_offset, l);
        str->l += l;
        fp->block_offset += l + 1;
        if (fp->block_offset >= fp->block_length) {
            fp->block_address = bgzf_htell(fp);
            fp->block_offset = 0;
            fp->block_length = 0;
        }
    } while (state == 0);
    if (state < -1) return state;
    if (str->l == 0 && state < 0) return state;
    fp->uncompressed_address += str->l + 1;
    if ( delim=='\n' && str->l>0 && str->s[str->l-1]=='\r' ) str->l--;
    str->s[str->l] = 0;
    return str->l <= INT_MAX ? (int) str->l : INT_MAX;
}

void bgzf_index_destroy(BGZF *fp)
{
    if ( !fp->idx ) return;
    free(fp->idx->offs);
    free(fp->idx);
    fp->idx = NULL;
    fp->idx_build_otf = 0;
}

int bgzf_index_build_init(BGZF *fp)
{
    bgzf_index_destroy(fp);
    fp->idx = (bgzidx_t*) calloc(1,sizeof(bgzidx_t));
    if ( !fp->idx ) return -1;
    fp->idx_build_otf = 1;  // build index on the fly
    return 0;
}

int bgzf_index_add_block(BGZF *fp)
{
    fp->idx->noffs++;
    if ( fp->idx->noffs > fp->idx->moffs )
    {
        fp->idx->moffs = fp->idx->noffs;
        kroundup32(fp->idx->moffs);
        fp->idx->offs = (bgzidx1_t*) realloc(fp->idx->offs, fp->idx->moffs*sizeof(bgzidx1_t));
        if ( !fp->idx->offs ) return -1;
    }
    fp->idx->offs[ fp->idx->noffs-1 ].uaddr = fp->idx->ublock_addr;
    fp->idx->offs[ fp->idx->noffs-1 ].caddr = fp->block_address;
    return 0;
}

static inline int hwrite_uint64(uint64_t x, hFILE *f)
{
    if (ed_is_big()) x = ed_swap_8(x);
    if (hwrite(f, &x, sizeof(x)) != sizeof(x)) return -1;
    return 0;
}

static char * get_name_suffix(const char *bname, const char *suffix)
{
    size_t len = strlen(bname) + strlen(suffix) + 1;
    char *buff = malloc(len);
    if (!buff) return NULL;
    snprintf(buff, len, "%s%s", bname, suffix);
    return buff;
}

int bgzf_index_dump_hfile(BGZF *fp, struct hFILE *idx, const char *name)
{
    // Note that the index contains one extra record when indexing files opened
    // for reading. The terminating record is not present when opened for writing.
    // This is not a bug.

    int i;

    if (!fp->idx) {
        hts_log_error("Called for BGZF handle with no index");
        errno = EINVAL;
        return -1;
    }

    if (bgzf_flush(fp) != 0) return -1;

    // discard the entry marking the end of the file
    if (fp->mt && fp->idx)
        fp->idx->noffs--;

    if (hwrite_uint64(fp->idx->noffs - 1, idx) < 0) goto fail;
    for (i=1; i<fp->idx->noffs; i++)
    {
        if (hwrite_uint64(fp->idx->offs[i].caddr, idx) < 0) goto fail;
        if (hwrite_uint64(fp->idx->offs[i].uaddr, idx) < 0) goto fail;
    }
    return 0;

 fail:
    hts_log_error("Error writing to %s : %s", name ? name : "index", strerror(errno));
    return -1;
}

int bgzf_index_dump(BGZF *fp, const char *bname, const char *suffix)
{
    const char *name = bname, *msg = NULL;
    char *tmp = NULL;
    hFILE *idx = NULL;

    if (!fp->idx) {
        hts_log_error("Called for BGZF handle with no index");
        errno = EINVAL;
        return -1;
    }

    if ( suffix )
    {
        tmp = get_name_suffix(bname, suffix);
        if ( !tmp ) return -1;
        name = tmp;
    }

    idx = hopen(name, "wb");
    if ( !idx ) {
        msg = "Error opening";
        goto fail;
    }

    if (bgzf_index_dump_hfile(fp, idx, name) != 0) goto fail;

    if (hclose(idx) < 0)
    {
        idx = NULL;
        msg = "Error on closing";
        goto fail;
    }

    free(tmp);
    return 0;

 fail:
    if (msg != NULL) {
        hts_log_error("%s %s : %s", msg, name, strerror(errno));
    }
    if (idx) hclose_abruptly(idx);
    free(tmp);
    return -1;
}

static inline int hread_uint64(uint64_t *xptr, hFILE *f)
{
    if (hread(f, xptr, sizeof(*xptr)) != sizeof(*xptr)) return -1;
    if (ed_is_big()) ed_swap_8p(xptr);
    return 0;
}

int bgzf_index_load_hfile(BGZF *fp, struct hFILE *idx, const char *name)
{
    fp->idx = (bgzidx_t*) calloc(1,sizeof(bgzidx_t));
    if (fp->idx == NULL) goto fail;
    uint64_t x;
    if (hread_uint64(&x, idx) < 0) goto fail;

    fp->idx->noffs = fp->idx->moffs = x + 1;
    fp->idx->offs  = (bgzidx1_t*) malloc(fp->idx->moffs*sizeof(bgzidx1_t));
    if (fp->idx->offs == NULL) goto fail;
    fp->idx->offs[0].caddr = fp->idx->offs[0].uaddr = 0;

    int i;
    for (i=1; i<fp->idx->noffs; i++)
    {
        if (hread_uint64(&fp->idx->offs[i].caddr, idx) < 0) goto fail;
        if (hread_uint64(&fp->idx->offs[i].uaddr, idx) < 0) goto fail;
    }

    return 0;

 fail:
    hts_log_error("Error reading %s : %s", name ? name : "index", strerror(errno));
    if (fp->idx) {
        free(fp->idx->offs);
        free(fp->idx);
        fp->idx = NULL;
    }
    return -1;
}

int bgzf_index_load(BGZF *fp, const char *bname, const char *suffix)
{
    const char *name = bname, *msg = NULL;
    char *tmp = NULL;
    hFILE *idx = NULL;
    if ( suffix )
    {
        tmp = get_name_suffix(bname, suffix);
        if ( !tmp ) return -1;
        name = tmp;
    }

    idx = hopen(name, "rb");
    if ( !idx ) {
        msg = "Error opening";
        goto fail;
    }

    if (bgzf_index_load_hfile(fp, idx, name) != 0) goto fail;

    if (hclose(idx) != 0) {
        idx = NULL;
        msg = "Error closing";
        goto fail;
    }

    free(tmp);
    return 0;

 fail:
    if (msg != NULL) {
        hts_log_error("%s %s : %s", msg, name, strerror(errno));
    }
    if (idx) hclose_abruptly(idx);
    free(tmp);
    return -1;
}

int bgzf_useek(BGZF *fp, off_t uoffset, int where)
{
    if (fp->is_write || where != SEEK_SET || fp->is_gzip) {
        fp->errcode |= BGZF_ERR_MISUSE;
        return -1;
    }
    if (uoffset >= fp->uncompressed_address - fp->block_offset &&
        uoffset < fp->uncompressed_address + fp->block_length - fp->block_offset) {
        // Can seek into existing data
        fp->block_offset += uoffset - fp->uncompressed_address;
        fp->uncompressed_address = uoffset;
        return 0;
    }
    if ( !fp->is_compressed )
    {
        if (hseek(fp->fp, uoffset, SEEK_SET) < 0)
        {
            fp->errcode |= BGZF_ERR_IO;
            return -1;
        }
        fp->block_length = 0;  // indicates current block has not been loaded
        fp->block_address = uoffset;
        fp->block_offset = 0;
        if (bgzf_read_block(fp) < 0) {
            fp->errcode |= BGZF_ERR_IO;
            return -1;
        }
        fp->uncompressed_address = uoffset;
        return 0;
    }

    if ( !fp->idx )
    {
        fp->errcode |= BGZF_ERR_IO;
        return -1;
    }

    // binary search
    int ilo = 0, ihi = fp->idx->noffs - 1;
    while ( ilo<=ihi )
    {
        int i = (ilo+ihi)*0.5;
        if ( uoffset < fp->idx->offs[i].uaddr ) ihi = i - 1;
        else if ( uoffset >= fp->idx->offs[i].uaddr ) ilo = i + 1;
        else break;
    }
    int i = ilo-1;
    off_t offset = 0;
    if (bgzf_seek_common(fp, fp->idx->offs[i].caddr, 0) < 0)
        return -1;

    if ( bgzf_read_block(fp) < 0 ) {
        fp->errcode |= BGZF_ERR_IO;
        return -1;
    }
    offset = uoffset - fp->idx->offs[i].uaddr;
    if ( offset > 0 )
    {
        if (offset > fp->block_length) {
            fp->errcode |= BGZF_ERR_IO;
            return -1;                                      //offset outside the available data
        }
        fp->block_offset = offset;
        assert( fp->block_offset <= fp->block_length );     // todo: skipped, unindexed, blocks
    }
    fp->uncompressed_address = uoffset;
    return 0;
}

off_t bgzf_utell(BGZF *fp)
{
    return fp->uncompressed_address;    // currently maintained only when reading
}

/* prototype is in hfile_internal.h */
struct hFILE *bgzf_hfile(struct BGZF *fp) {
    return fp->fp;
}
