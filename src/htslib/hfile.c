/*  hfile.c -- buffered low-level input/output streams.

    Copyright (C) 2013-2021 Genome Research Ltd.

    Author: John Marshall <jm18@sanger.ac.uk>

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

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <errno.h>
#include <limits.h>

#include <pthread.h>

#ifdef ENABLE_PLUGINS
#if defined(_WIN32) || defined(__CYGWIN__) || defined(__MSYS__)
#define USING_WINDOWS_PLUGIN_DLLS
#include <dlfcn.h>
#endif
#endif

#include "htslib/hfile.h"
#include "hfile_internal.h"
#include "htslib/kstring.h"

#ifndef ENOTSUP
#define ENOTSUP EINVAL
#endif
#ifndef EOVERFLOW
#define EOVERFLOW ERANGE
#endif
#ifndef EPROTONOSUPPORT
#define EPROTONOSUPPORT ENOSYS
#endif

#ifndef SSIZE_MAX /* SSIZE_MAX is POSIX 1 */
#define SSIZE_MAX LONG_MAX
#endif

/* hFILE fields are used as follows:

   char *buffer;     // Pointer to the start of the I/O buffer
   char *begin;      // First not-yet-read character / unused position
   char *end;        // First unfilled/unfillable position
   char *limit;      // Pointer to the first position past the buffer

   const hFILE_backend *backend;  // Methods to refill/flush I/O buffer

   off_t offset;     // Offset within the stream of buffer position 0
   unsigned at_eof:1;// For reading, whether EOF has been seen
   unsigned mobile:1;// Buffer is a mobile window or fixed full contents
   unsigned readonly:1;// Whether opened as "r" rather than "r+"/"w"/"a"
   int has_errno;    // Error number from the last failure on this stream

For reading, begin is the first unread character in the buffer and end is the
first unfilled position:

   -----------ABCDEFGHIJKLMNO---------------
   ^buffer    ^begin         ^end           ^limit

For writing, begin is the first unused position and end is unused so remains
equal to buffer:

   ABCDEFGHIJKLMNOPQRSTUVWXYZ---------------
   ^buffer                   ^begin         ^limit
   ^end

Thus if begin > end then there is a non-empty write buffer, if begin < end
then there is a non-empty read buffer, and if begin == end then both buffers
are empty.  In all cases, the stream's file position indicator corresponds
to the position pointed to by begin.

The above is the normal scenario of a mobile window.  For in-memory
streams (eg via hfile_init_fixed) the buffer can be used as the full
contents without any separate backend behind it.  These always have at_eof
set, offset set to 0, need no read() method, and should just return EINVAL
for seek():

   abcdefghijkLMNOPQRSTUVWXYZ------
   ^buffer    ^begin         ^end  ^limit
*/
HTSLIB_EXPORT
hFILE *hfile_init(size_t struct_size, const char *mode, size_t capacity)
{
    hFILE *fp = (hFILE *) malloc(struct_size);
    if (fp == NULL) goto error;

    if (capacity == 0) capacity = 32768;
    // FIXME For now, clamp input buffer sizes so mpileup doesn't eat memory
    if (strchr(mode, 'r') && capacity > 32768) capacity = 32768;

    fp->buffer = (char *) malloc(capacity);
    if (fp->buffer == NULL) goto error;

    fp->begin = fp->end = fp->buffer;
    fp->limit = &fp->buffer[capacity];

    fp->offset = 0;
    fp->at_eof = 0;
    fp->mobile = 1;
    fp->readonly = (strchr(mode, 'r') && ! strchr(mode, '+'));
    fp->has_errno = 0;
    return fp;

error:
    hfile_destroy(fp);
    return NULL;
}

hFILE *hfile_init_fixed(size_t struct_size, const char *mode,
                        char *buffer, size_t buf_filled, size_t buf_size)
{
    hFILE *fp = (hFILE *) malloc(struct_size);
    if (fp == NULL) return NULL;

    fp->buffer = fp->begin = buffer;
    fp->end = &fp->buffer[buf_filled];
    fp->limit = &fp->buffer[buf_size];

    fp->offset = 0;
    fp->at_eof = 1;
    fp->mobile = 0;
    fp->readonly = (strchr(mode, 'r') && ! strchr(mode, '+'));
    fp->has_errno = 0;
    return fp;
}

static const struct hFILE_backend mem_backend;

HTSLIB_EXPORT
void hfile_destroy(hFILE *fp)
{
    int save = errno;
    if (fp) free(fp->buffer);
    free(fp);
    errno = save;
}

static inline int writebuffer_is_nonempty(hFILE *fp)
{
    return fp->begin > fp->end;
}

/* Refills the read buffer from the backend (once, so may only partially
   fill the buffer), returning the number of additional characters read
   (which might be 0), or negative when an error occurred.  */
static ssize_t refill_buffer(hFILE *fp)
{
    ssize_t n;

    // Move any unread characters to the start of the buffer
    if (fp->mobile && fp->begin > fp->buffer) {
        fp->offset += fp->begin - fp->buffer;
        memmove(fp->buffer, fp->begin, fp->end - fp->begin);
        fp->end = &fp->buffer[fp->end - fp->begin];
        fp->begin = fp->buffer;
    }

    // Read into the available buffer space at fp->[end,limit)
    if (fp->at_eof || fp->end == fp->limit) n = 0;
    else {
        n = fp->backend->read(fp, fp->end, fp->limit - fp->end);
        if (n < 0) { fp->has_errno = errno; return n; }
        else if (n == 0) fp->at_eof = 1;
    }

    fp->end += n;
    return n;
}

/*
 * Changes the buffer size for an hFILE.  Ideally this is done
 * immediately after opening.  If performed later, this function may
 * fail if we are reducing the buffer size and the current offset into
 * the buffer is beyond the new capacity.
 *
 * Returns 0 on success;
 *        -1 on failure.
 */
HTSLIB_EXPORT
int hfile_set_blksize(hFILE *fp, size_t bufsiz) {
    char *buffer;
    ptrdiff_t curr_used;
    if (!fp) return -1;
    curr_used = (fp->begin > fp->end ? fp->begin : fp->end) - fp->buffer;
    if (bufsiz == 0) bufsiz = 32768;

    // Ensure buffer resize will not erase live data
    if (bufsiz < curr_used)
        return -1;

    if (!(buffer = (char *) realloc(fp->buffer, bufsiz))) return -1;

    fp->begin  = buffer + (fp->begin - fp->buffer);
    fp->end    = buffer + (fp->end   - fp->buffer);
    fp->buffer = buffer;
    fp->limit  = &fp->buffer[bufsiz];

    return 0;
}

/* Called only from hgetc(), when our buffer is empty.  */
HTSLIB_EXPORT
int hgetc2(hFILE *fp)
{
    return (refill_buffer(fp) > 0)? (unsigned char) *(fp->begin++) : EOF;
}

ssize_t hgetdelim(char *buffer, size_t size, int delim, hFILE *fp)
{
    char *found;
    size_t n, copied = 0;
    ssize_t got;

    if (size < 1 || size > SSIZE_MAX) {
        fp->has_errno = errno = EINVAL;
        return -1;
    }
    if (writebuffer_is_nonempty(fp)) {
        fp->has_errno = errno = EBADF;
        return -1;
    }

    --size; /* to allow space for the NUL terminator */

    do {
        n = fp->end - fp->begin;
        if (n > size - copied) n = size - copied;

        /* Look in the hFILE buffer for the delimiter */
        found = memchr(fp->begin, delim, n);
        if (found != NULL) {
            n = found - fp->begin + 1;
            memcpy(buffer + copied, fp->begin, n);
            buffer[n + copied] = '\0';
            fp->begin += n;
            return n + copied;
        }

        /* No delimiter yet, copy as much as we can and refill if necessary */
        memcpy(buffer + copied, fp->begin, n);
        fp->begin += n;
        copied += n;

        if (copied == size) { /* Output buffer full */
            buffer[copied] = '\0';
            return copied;
        }

        got = refill_buffer(fp);
    } while (got > 0);

    if (got < 0) return -1; /* Error on refill. */

    buffer[copied] = '\0';  /* EOF, return anything that was copied. */
    return copied;
}

char *hgets(char *buffer, int size, hFILE *fp)
{
    if (size < 1) {
        fp->has_errno = errno = EINVAL;
        return NULL;
    }
    return hgetln(buffer, size, fp) > 0 ? buffer : NULL;
}

ssize_t hpeek(hFILE *fp, void *buffer, size_t nbytes)
{
    size_t n = fp->end - fp->begin;
    while (n < nbytes) {
        ssize_t ret = refill_buffer(fp);
        if (ret < 0) return ret;
        else if (ret == 0) break;
        else n += ret;
    }

    if (n > nbytes) n = nbytes;
    memcpy(buffer, fp->begin, n);
    return n;
}

/* Called only from hread(); when called, our buffer is empty and nread bytes
   have already been placed in the destination buffer.  */
HTSLIB_EXPORT
ssize_t hread2(hFILE *fp, void *destv, size_t nbytes, size_t nread)
{
    const size_t capacity = fp->limit - fp->buffer;
    int buffer_invalidated = 0;
    char *dest = (char *) destv;
    dest += nread, nbytes -= nread;

    // Read large requests directly into the destination buffer
    while (nbytes * 2 >= capacity && !fp->at_eof) {
        ssize_t n = fp->backend->read(fp, dest, nbytes);
        if (n < 0) { fp->has_errno = errno; return n; }
        else if (n == 0) fp->at_eof = 1;
        else buffer_invalidated = 1;
        fp->offset += n;
        dest += n, nbytes -= n;
        nread += n;
    }

    if (buffer_invalidated) {
        // Our unread buffer is empty, so begin == end, but our already-read
        // buffer [buffer,begin) is likely non-empty and is no longer valid as
        // its contents are no longer adjacent to the file position indicator.
        // Discard it so that hseek() can't try to take advantage of it.
        fp->offset += fp->begin - fp->buffer;
        fp->begin = fp->end = fp->buffer;
    }

    while (nbytes > 0 && !fp->at_eof) {
        size_t n;
        ssize_t ret = refill_buffer(fp);
        if (ret < 0) return ret;

        n = fp->end - fp->begin;
        if (n > nbytes) n = nbytes;
        memcpy(dest, fp->begin, n);
        fp->begin += n;
        dest += n, nbytes -= n;
        nread += n;
    }

    return nread;
}

/* Flushes the write buffer, fp->[buffer,begin), out through the backend
   returning 0 on success or negative if an error occurred.  */
static ssize_t flush_buffer(hFILE *fp)
{
    const char *buffer = fp->buffer;
    while (buffer < fp->begin) {
        ssize_t n = fp->backend->write(fp, buffer, fp->begin - buffer);
        if (n < 0) { fp->has_errno = errno; return n; }
        buffer += n;
        fp->offset += n;
    }

    fp->begin = fp->buffer;  // Leave the buffer empty
    return 0;
}

int hflush(hFILE *fp)
{
    if (flush_buffer(fp) < 0) return EOF;
    if (fp->backend->flush) {
        if (fp->backend->flush(fp) < 0) { fp->has_errno = errno; return EOF; }
    }
    return 0;
}

/* Called only from hputc(), when our buffer is already full.  */
HTSLIB_EXPORT
int hputc2(int c, hFILE *fp)
{
    if (flush_buffer(fp) < 0) return EOF;
    *(fp->begin++) = c;
    return c;
}

/* Called only from hwrite() and hputs2(); when called, our buffer is either
   full and ncopied bytes from the source have already been copied to our
   buffer; or completely empty, ncopied is zero and totalbytes is greater than
   the buffer size.  */
HTSLIB_EXPORT
ssize_t hwrite2(hFILE *fp, const void *srcv, size_t totalbytes, size_t ncopied)
{
    const char *src = (const char *) srcv;
    ssize_t ret;
    const size_t capacity = fp->limit - fp->buffer;
    size_t remaining = totalbytes - ncopied;
    src += ncopied;

    ret = flush_buffer(fp);
    if (ret < 0) return ret;

    // Write large blocks out directly from the source buffer
    while (remaining * 2 >= capacity) {
        ssize_t n = fp->backend->write(fp, src, remaining);
        if (n < 0) { fp->has_errno = errno; return n; }
        fp->offset += n;
        src += n, remaining -= n;
    }

    // Just buffer any remaining characters
    memcpy(fp->begin, src, remaining);
    fp->begin += remaining;

    return totalbytes;
}

/* Called only from hputs(), when our buffer is already full.  */
HTSLIB_EXPORT
int hputs2(const char *text, size_t totalbytes, size_t ncopied, hFILE *fp)
{
    return (hwrite2(fp, text, totalbytes, ncopied) >= 0)? 0 : EOF;
}

off_t hseek(hFILE *fp, off_t offset, int whence)
{
    off_t curpos, pos;

    if (writebuffer_is_nonempty(fp) && fp->mobile) {
        int ret = flush_buffer(fp);
        if (ret < 0) return ret;
    }

    curpos = htell(fp);

    // Relative offsets are given relative to the hFILE's stream position,
    // which may differ from the backend's physical position due to buffering
    // read-ahead.  Correct for this by converting to an absolute position.
    if (whence == SEEK_CUR) {
        if (curpos + offset < 0) {
            // Either a negative offset resulted in a position before the
            // start of the file, or we overflowed when given a positive offset
            fp->has_errno = errno = (offset < 0)? EINVAL : EOVERFLOW;
            return -1;
        }

        whence = SEEK_SET;
        offset = curpos + offset;
    }
    // For fixed immobile buffers, convert everything else to SEEK_SET too
    // so that seeking can be avoided for all (within range) requests.
    else if (! fp->mobile && whence == SEEK_END) {
        size_t length = fp->end - fp->buffer;
        if (offset > 0 || -offset > length) {
            fp->has_errno = errno = EINVAL;
            return -1;
        }

        whence = SEEK_SET;
        offset = length + offset;
    }

    // Avoid seeking if the desired position is within our read buffer.
    // (But not when the next operation may be a write on a mobile buffer.)
    if (whence == SEEK_SET && (! fp->mobile || fp->readonly) &&
        offset >= fp->offset && offset - fp->offset <= fp->end - fp->buffer) {
        fp->begin = &fp->buffer[offset - fp->offset];
        return offset;
    }

    pos = fp->backend->seek(fp, offset, whence);
    if (pos < 0) { fp->has_errno = errno; return pos; }

    // Seeking succeeded, so discard any non-empty read buffer
    fp->begin = fp->end = fp->buffer;
    fp->at_eof = 0;

    fp->offset = pos;
    return pos;
}

int hclose(hFILE *fp)
{
    int err = fp->has_errno;

    if (writebuffer_is_nonempty(fp) && hflush(fp) < 0) err = fp->has_errno;
    if (fp->backend->close(fp) < 0) err = errno;
    hfile_destroy(fp);

    if (err) {
        errno = err;
        return EOF;
    }
    else return 0;
}

void hclose_abruptly(hFILE *fp)
{
    int save = errno;
    if (fp->backend->close(fp) < 0) { /* Ignore subsequent errors */ }
    hfile_destroy(fp);
    errno = save;
}


/***************************
 * File descriptor backend *
 ***************************/

#ifndef _WIN32
#include <sys/socket.h>
#include <sys/stat.h>
#define HAVE_STRUCT_STAT_ST_BLKSIZE
#else
#include <winsock2.h>
#define HAVE_CLOSESOCKET
#define HAVE_SETMODE
#endif
#include <fcntl.h>
#include <unistd.h>

/* For Unix, it doesn't matter whether a file descriptor is a socket.
   However Windows insists on send()/recv() and its own closesocket()
   being used when fd happens to be a socket.  */

typedef struct {
    hFILE base;
    int fd;
    unsigned is_socket:1, is_shared:1;
} hFILE_fd;

static ssize_t fd_read(hFILE *fpv, void *buffer, size_t nbytes)
{
    hFILE_fd *fp = (hFILE_fd *) fpv;
    ssize_t n;
    do {
        n = fp->is_socket? recv(fp->fd, buffer, nbytes, 0)
                         : read(fp->fd, buffer, nbytes);
    } while (n < 0 && errno == EINTR);
    return n;
}

static ssize_t fd_write(hFILE *fpv, const void *buffer, size_t nbytes)
{
    hFILE_fd *fp = (hFILE_fd *) fpv;
    ssize_t n;
    do {
        n = fp->is_socket?  send(fp->fd, buffer, nbytes, 0)
                         : write(fp->fd, buffer, nbytes);
    } while (n < 0 && errno == EINTR);
#ifdef _WIN32
        // On windows we have no SIGPIPE.  Instead write returns
        // EINVAL.  We check for this and our fd being a pipe.
        // If so, we raise SIGTERM instead of SIGPIPE.  It's not
        // ideal, but I think the only alternative is extra checking
        // in every single piece of code.
        if (n < 0 && errno == EINVAL &&
            GetLastError() == ERROR_NO_DATA &&
            GetFileType((HANDLE)_get_osfhandle(fp->fd)) == FILE_TYPE_PIPE) {
            raise(SIGTERM);
        }
#endif
    return n;
}

static off_t fd_seek(hFILE *fpv, off_t offset, int whence)
{
    hFILE_fd *fp = (hFILE_fd *) fpv;
#ifdef _WIN32
    // On windows lseek can return non-zero values even on a pipe.  Instead
    // it's likely to seek somewhere within the pipe memory buffer.
    // This breaks bgzf_check_EOF among other things.
    if (GetFileType((HANDLE)_get_osfhandle(fp->fd)) == FILE_TYPE_PIPE) {
        errno = ESPIPE;
        return -1;
    }
#endif

    return lseek(fp->fd, offset, whence);
}

static int fd_flush(hFILE *fpv)
{
    int ret = 0;
    do {
#ifdef HAVE_FDATASYNC
        hFILE_fd *fp = (hFILE_fd *) fpv;
        ret = fdatasync(fp->fd);
#elif defined(HAVE_FSYNC)
        hFILE_fd *fp = (hFILE_fd *) fpv;
        ret = fsync(fp->fd);
#endif
        // Ignore invalid-for-fsync(2) errors due to being, e.g., a pipe,
        // and operation-not-supported errors (Mac OS X)
        if (ret < 0 && (errno == EINVAL || errno == ENOTSUP)) ret = 0;
    } while (ret < 0 && errno == EINTR);
    return ret;
}

static int fd_close(hFILE *fpv)
{
    hFILE_fd *fp = (hFILE_fd *) fpv;
    int ret;

    // If we don't own the fd, return successfully without actually closing it
    if (fp->is_shared) return 0;

    do {
#ifdef HAVE_CLOSESOCKET
        ret = fp->is_socket? closesocket(fp->fd) : close(fp->fd);
#else
        ret = close(fp->fd);
#endif
    } while (ret < 0 && errno == EINTR);
    return ret;
}

static const struct hFILE_backend fd_backend =
{
    fd_read, fd_write, fd_seek, fd_flush, fd_close
};

static size_t blksize(int fd)
{
#ifdef HAVE_STRUCT_STAT_ST_BLKSIZE
    struct stat sbuf;
    if (fstat(fd, &sbuf) != 0) return 0;
    return sbuf.st_blksize;
#else
    return 0;
#endif
}

static hFILE *hopen_fd(const char *filename, const char *mode)
{
    hFILE_fd *fp = NULL;
    int fd = open(filename, hfile_oflags(mode), 0666);
    if (fd < 0) goto error;

    fp = (hFILE_fd *) hfile_init(sizeof (hFILE_fd), mode, blksize(fd));
    if (fp == NULL) goto error;

    fp->fd = fd;
    fp->is_socket = 0;
    fp->is_shared = 0;
    fp->base.backend = &fd_backend;
    return &fp->base;

error:
    if (fd >= 0) { int save = errno; (void) close(fd); errno = save; }
    hfile_destroy((hFILE *) fp);
    return NULL;
}

// Loads the contents of filename to produced a read-only, in memory,
// immobile hfile.  fp is the already opened file.  We always close this
// input fp, irrespective of whether we error or whether we return a new
// immobile hfile.
static hFILE *hpreload(hFILE *fp) {
    hFILE *mem_fp;
    char *buf = NULL;
    off_t buf_sz = 0, buf_a = 0, buf_inc = 8192, len;

    for (;;) {
        if (buf_a - buf_sz < 5000) {
            buf_a += buf_inc;
            char *t = realloc(buf, buf_a);
            if (!t) goto err;
            buf = t;
            if (buf_inc < 1000000) buf_inc *= 1.3;
        }
        len = hread(fp, buf+buf_sz, buf_a-buf_sz);
        if (len > 0)
            buf_sz += len;
        else
            break;
    }

    if (len < 0) goto err;
    mem_fp = hfile_init_fixed(sizeof(hFILE), "r", buf, buf_sz, buf_a);
    if (!mem_fp) goto err;
    mem_fp->backend = &mem_backend;

    if (hclose(fp) < 0) {
        hclose_abruptly(mem_fp);
        goto err;
    }
    return mem_fp;

 err:
    free(buf);
    hclose_abruptly(fp);
    return NULL;
}

static int is_preload_url_remote(const char *url){
    return hisremote(url + 8); // len("preload:") = 8
}

static hFILE *hopen_preload(const char *url, const char *mode){
    hFILE* fp = hopen(url + 8, mode);
    return hpreload(fp);
}

hFILE *hdopen(int fd, const char *mode)
{
    hFILE_fd *fp = (hFILE_fd*) hfile_init(sizeof (hFILE_fd), mode, blksize(fd));
    if (fp == NULL) return NULL;

    fp->fd = fd;
    fp->is_socket = (strchr(mode, 's') != NULL);
    fp->is_shared = (strchr(mode, 'S') != NULL);
    fp->base.backend = &fd_backend;
    return &fp->base;
}

static hFILE *hopen_fd_fileuri(const char *url, const char *mode)
{
    if (strncmp(url, "file://localhost/", 17) == 0) url += 16;
    else if (strncmp(url, "file:///", 8) == 0) url += 7;
    else { errno = EPROTONOSUPPORT; return NULL; }

#if defined(_WIN32) || defined(__MSYS__)
    // For cases like C:/foo
    if (url[0] == '/' && url[1] && url[2] == ':' && url[3] == '/') url++;
#endif

    return hopen_fd(url, mode);
}

static hFILE *hopen_fd_stdinout(const char *mode)
{
    int fd = (strchr(mode, 'r') != NULL)? STDIN_FILENO : STDOUT_FILENO;
    char mode_shared[101];
    snprintf(mode_shared, sizeof mode_shared, "S%s", mode);
#if defined HAVE_SETMODE && defined O_BINARY
    if (setmode(fd, O_BINARY) < 0) return NULL;
#endif
    return hdopen(fd, mode_shared);
}

HTSLIB_EXPORT
int hfile_oflags(const char *mode)
{
    int rdwr = 0, flags = 0;
    const char *s;
    for (s = mode; *s; s++)
        switch (*s) {
        case 'r': rdwr = O_RDONLY;  break;
        case 'w': rdwr = O_WRONLY; flags |= O_CREAT | O_TRUNC;  break;
        case 'a': rdwr = O_WRONLY; flags |= O_CREAT | O_APPEND;  break;
        case '+': rdwr = O_RDWR;  break;
#ifdef O_CLOEXEC
        case 'e': flags |= O_CLOEXEC;  break;
#endif
#ifdef O_EXCL
        case 'x': flags |= O_EXCL;  break;
#endif
        default:  break;
        }

#ifdef O_BINARY
    flags |= O_BINARY;
#endif

    return rdwr | flags;
}


/*********************
 * In-memory backend *
 *********************/

#include "hts_internal.h"

typedef struct {
    hFILE base;
} hFILE_mem;

static off_t mem_seek(hFILE *fpv, off_t offset, int whence)
{
    errno = EINVAL;
    return -1;
}

static int mem_close(hFILE *fpv)
{
    return 0;
}

static const struct hFILE_backend mem_backend =
{
    NULL, NULL, mem_seek, NULL, mem_close
};

static int cmp_prefix(const char *key, const char *s)
{
    while (*key)
        if (tolower_c(*s) != *key) return +1;
        else s++, key++;

    return 0;
}

static hFILE *create_hfile_mem(char* buffer, const char* mode, size_t buf_filled, size_t buf_size)
{
    hFILE_mem *fp = (hFILE_mem *) hfile_init_fixed(sizeof(hFILE_mem), mode, buffer, buf_filled, buf_size);
    if (fp == NULL)
        return NULL;

    fp->base.backend = &mem_backend;
    return &fp->base;
}

static hFILE *hopen_mem(const char *url, const char *mode)
{
    size_t length, size;
    char *buffer;
    const char *data, *comma = strchr(url, ',');
    if (comma == NULL) { errno = EINVAL; return NULL; }
    data = comma+1;

    // TODO Implement write modes
    if (strchr(mode, 'r') == NULL) { errno = EROFS; return NULL; }

    if (comma - url >= 7 && cmp_prefix(";base64", &comma[-7]) == 0) {
        size = hts_base64_decoded_length(strlen(data));
        buffer = malloc(size);
        if (buffer == NULL) return NULL;
        hts_decode_base64(buffer, &length, data);
    }
    else {
        size = strlen(data) + 1;
        buffer = malloc(size);
        if (buffer == NULL) return NULL;
        hts_decode_percent(buffer, &length, data);
    }
    hFILE* hf;

    if(!(hf = create_hfile_mem(buffer, mode, length, size))){
        free(buffer);
        return NULL;
    }

    return hf;
}

static hFILE *hopenv_mem(const char *filename, const char *mode, va_list args)
{
    char* buffer = va_arg(args, char*);
    size_t sz = va_arg(args, size_t);
    va_end(args);

    hFILE* hf;

    if(!(hf = create_hfile_mem(buffer, mode, sz, sz))){
        free(buffer);
        return NULL;
    }

    return hf;
}

char *hfile_mem_get_buffer(hFILE *file, size_t *length) {
    if (file->backend != &mem_backend) {
        errno = EINVAL;
        return NULL;
    }

    if (length)
        *length = file->buffer - file->limit;

    return file->buffer;
}

char *hfile_mem_steal_buffer(hFILE *file, size_t *length) {
    char *buf = hfile_mem_get_buffer(file, length);
    if (buf)
        file->buffer = NULL;
    return buf;
}

int hfile_plugin_init_mem(struct hFILE_plugin *self)
{
    // mem files are declared remote so they work with a tabix index
    static const struct hFILE_scheme_handler handler =
            {NULL, hfile_always_remote, "mem", 2000 + 50, hopenv_mem};
    self->name = "mem";
    hfile_add_scheme_handler("mem", &handler);
    return 0;
}

/**********************************************************************
 * Dummy crypt4gh plug-in.  Does nothing apart from advise how to get *
 * the real one.  It will be overridden by the actual plug-in.        *
 **********************************************************************/

static hFILE *crypt4gh_needed(const char *url, const char *mode)
{
    const char *u = strncmp(url, "crypt4gh:", 9) == 0 ? url + 9 : url;
#if defined(ENABLE_PLUGINS)
    const char *enable_plugins = "";
#else
    const char *enable_plugins = "You also need to rebuild HTSlib with plug-ins enabled.\n";
#endif

    hts_log_error("Accessing \"%s\" needs the crypt4gh plug-in.\n"
                  "It can be found at "
                  "https://github.com/samtools/htslib-crypt4gh\n"
                  "%s"
                  "If you have the plug-in, please ensure it can be "
                  "found on your HTS_PATH.",
                  u, enable_plugins);

    errno = EPROTONOSUPPORT;
    return NULL;
}

int hfile_plugin_init_crypt4gh_needed(struct hFILE_plugin *self)
{
    static const struct hFILE_scheme_handler handler =
        { crypt4gh_needed, NULL, "crypt4gh-needed", 0, NULL };
    self->name = "crypt4gh-needed";
    hfile_add_scheme_handler("crypt4gh", &handler);
    return 0;
}


/*****************************************
 * Plugin and hopen() backend dispatcher *
 *****************************************/

#include "htslib/khash.h"

KHASH_MAP_INIT_STR(scheme_string, const struct hFILE_scheme_handler *)
static khash_t(scheme_string) *schemes = NULL;

struct hFILE_plugin_list {
    struct hFILE_plugin plugin;
    struct hFILE_plugin_list *next;
};

static struct hFILE_plugin_list *plugins = NULL;
static pthread_mutex_t plugins_lock = PTHREAD_MUTEX_INITIALIZER;

void hfile_shutdown(int do_close_plugin)
{
    pthread_mutex_lock(&plugins_lock);

    if (schemes) {
        kh_destroy(scheme_string, schemes);
        schemes = NULL;
    }

    while (plugins != NULL) {
        struct hFILE_plugin_list *p = plugins;
        if (p->plugin.destroy) p->plugin.destroy();
#ifdef ENABLE_PLUGINS
        if (p->plugin.obj && do_close_plugin) close_plugin(p->plugin.obj);
#endif
        plugins = p->next;
        free(p);
    }

    pthread_mutex_unlock(&plugins_lock);
}

static void hfile_exit()
{
    hfile_shutdown(0);
    pthread_mutex_destroy(&plugins_lock);
}

static inline int priority(const struct hFILE_scheme_handler *handler)
{
    return handler->priority % 1000;
}

#ifdef USING_WINDOWS_PLUGIN_DLLS
/*
 * Work-around for Windows plug-in dlls where the plug-in could be
 * using a different HTSlib library to the executable (for example
 * because the latter was build against a static libhts.a).  When this
 * happens, the plug-in can call the wrong copy of hfile_add_scheme_handler().
 * If this is detected, it calls this function which attempts to fix the
 * problem by redirecting to the hfile_add_scheme_handler() in the main
 * executable.
 */
static int try_exe_add_scheme_handler(const char *scheme,
                                      const struct hFILE_scheme_handler *handler)
{
    static void (*add_scheme_handler)(const char *scheme,
                                      const struct hFILE_scheme_handler *handler);
    if (!add_scheme_handler) {
        // dlopen the main executable and resolve hfile_add_scheme_handler
        void *exe_handle = dlopen(NULL, RTLD_LAZY);
        if (!exe_handle) return -1;
        *(void **) (&add_scheme_handler) = dlsym(exe_handle, "hfile_add_scheme_handler");
        dlclose(exe_handle);
    }
    // Check that the symbol was obtained and isn't the one in this copy
    // of the library (to avoid infinite recursion)
    if (!add_scheme_handler || add_scheme_handler == hfile_add_scheme_handler)
        return -1;
    add_scheme_handler(scheme, handler);
    return 0;
}
#else
static int try_exe_add_scheme_handler(const char *scheme,
                                      const struct hFILE_scheme_handler *handler)
{
    return -1;
}
#endif

HTSLIB_EXPORT
void hfile_add_scheme_handler(const char *scheme,
                              const struct hFILE_scheme_handler *handler)
{
    int absent;
    if (!schemes) {
        if (try_exe_add_scheme_handler(scheme, handler) != 0) {
            hts_log_warning("Couldn't register scheme handler for %s", scheme);
        }
        return;
    }
    khint_t k = kh_put(scheme_string, schemes, scheme, &absent);
    if (absent < 0) {
        hts_log_warning("Couldn't register scheme handler for %s : %s",
                        scheme, strerror(errno));
        return;
    }
    if (absent || priority(handler) > priority(kh_value(schemes, k))) {
        kh_value(schemes, k) = handler;
    }
}

static int init_add_plugin(void *obj, int (*init)(struct hFILE_plugin *),
                           const char *pluginname)
{
    struct hFILE_plugin_list *p = malloc (sizeof (struct hFILE_plugin_list));
    if (p == NULL) {
        hts_log_debug("Failed to allocate memory for plugin \"%s\"", pluginname);
        return -1;
    }

    p->plugin.api_version = 1;
    p->plugin.obj = obj;
    p->plugin.name = NULL;
    p->plugin.destroy = NULL;

    int ret = (*init)(&p->plugin);

    if (ret != 0) {
        hts_log_debug("Initialisation failed for plugin \"%s\": %d", pluginname, ret);
        free(p);
        return ret;
    }

    hts_log_debug("Loaded \"%s\"", pluginname);

    p->next = plugins, plugins = p;
    return 0;
}

/*
 * Returns 0 on success,
 *        <0 on failure
 */
static int load_hfile_plugins()
{
    static const struct hFILE_scheme_handler
        data = { hopen_mem, hfile_always_local, "built-in", 80 },
        file = { hopen_fd_fileuri, hfile_always_local, "built-in", 80 },
        preload = { hopen_preload, is_preload_url_remote, "built-in", 80 };

    schemes = kh_init(scheme_string);
    if (schemes == NULL)
        return -1;

    hfile_add_scheme_handler("data", &data);
    hfile_add_scheme_handler("file", &file);
    hfile_add_scheme_handler("preload", &preload);
    init_add_plugin(NULL, hfile_plugin_init_mem, "mem");
    init_add_plugin(NULL, hfile_plugin_init_crypt4gh_needed, "crypt4gh-needed");

#ifdef ENABLE_PLUGINS
    struct hts_path_itr path;
    const char *pluginname;
    hts_path_itr_setup(&path, NULL, NULL, "hfile_", 6, NULL, 0);
    while ((pluginname = hts_path_itr_next(&path)) != NULL) {
        void *obj;
        int (*init)(struct hFILE_plugin *) = (int (*)(struct hFILE_plugin *))
            load_plugin(&obj, pluginname, "hfile_plugin_init");

        if (init) {
            if (init_add_plugin(obj, init, pluginname) != 0)
                close_plugin(obj);
        }
    }
#else

#ifdef ENABLE_GCS
    init_add_plugin(NULL, hfile_plugin_init_gcs, "gcs");
#endif
#ifdef ENABLE_S3
    init_add_plugin(NULL, hfile_plugin_init_s3, "s3");
    init_add_plugin(NULL, hfile_plugin_init_s3_write, "s3w");
#endif

#endif

    // In the unlikely event atexit() fails, it's better to succeed here and
    // carry on; then eventually when the program exits, we'll merely close
    // down the plugins uncleanly, as if we had aborted.
    (void) atexit(hfile_exit);

    return 0;
}

/* A filename like "foo:bar" in which we don't recognise the scheme is
   either an ordinary file or an indication of a missing or broken plugin.
   Try to open it as an ordinary file; but if there's no such file, set
   errno distinctively to make the plugin issue apparent.  */
static hFILE *hopen_unknown_scheme(const char *fname, const char *mode)
{
    hFILE *fp = hopen_fd(fname, mode);
    if (fp == NULL && errno == ENOENT) errno = EPROTONOSUPPORT;
    return fp;
}

/* Returns the appropriate handler, or NULL if the string isn't an URL.  */
static const struct hFILE_scheme_handler *find_scheme_handler(const char *s)
{
    static const struct hFILE_scheme_handler unknown_scheme =
        { hopen_unknown_scheme, hfile_always_local, "built-in", 0 };

    char scheme[12];
    int i;

    for (i = 0; i < sizeof scheme; i++)
        if (isalnum_c(s[i]) || s[i] == '+' || s[i] == '-' || s[i] == '.')
            scheme[i] = tolower_c(s[i]);
        else if (s[i] == ':') break;
        else return NULL;

    // 1 byte schemes are likely windows C:/foo pathnames
    if (i <= 1 || i >= sizeof scheme) return NULL;
    scheme[i] = '\0';

    pthread_mutex_lock(&plugins_lock);
    if (!schemes && load_hfile_plugins() < 0) {
        pthread_mutex_unlock(&plugins_lock);
        return NULL;
    }
    pthread_mutex_unlock(&plugins_lock);

    khint_t k = kh_get(scheme_string, schemes, scheme);
    return (k != kh_end(schemes))? kh_value(schemes, k) : &unknown_scheme;
}


/***************************
 * Library introspection functions
 ***************************/

/*
 * Fills out sc_list[] with the list of known URL schemes.
 * This can be restricted to just ones from a specific plugin,
 * or all (plugin == NULL).
 *
 * Returns number of schemes found on success;
 *        -1 on failure.
 */
HTSLIB_EXPORT
int hfile_list_schemes(const char *plugin, const char *sc_list[], int *nschemes)
{
    pthread_mutex_lock(&plugins_lock);
    if (!schemes && load_hfile_plugins() < 0) {
        pthread_mutex_unlock(&plugins_lock);
        return -1;
    }
    pthread_mutex_unlock(&plugins_lock);

    khiter_t k;
    int ns = 0;

    for (k = kh_begin(schemes); k != kh_end(schemes); k++) {
        if (!kh_exist(schemes, k))
            continue;

        const struct hFILE_scheme_handler *s = kh_value(schemes, k);
        if (plugin && strcmp(s->provider, plugin) != 0)
            continue;

        if (ns < *nschemes)
            sc_list[ns] = kh_key(schemes, k);
        ns++;
    }

    if (*nschemes > ns)
        *nschemes = ns;

    return ns;
}


/*
 * Fills out plist[] with the list of known hFILE plugins.
 *
 * Returns number of schemes found on success;
 *        -1 on failure
 */
HTSLIB_EXPORT
int hfile_list_plugins(const char *plist[], int *nplugins)
{
    pthread_mutex_lock(&plugins_lock);
    if (!schemes && load_hfile_plugins() < 0) {
        pthread_mutex_unlock(&plugins_lock);
        return -1;
    }
    pthread_mutex_unlock(&plugins_lock);

    int np = 0;
    if (*nplugins)
        plist[np++] = "built-in";

    struct hFILE_plugin_list *p = plugins;
    while (p) {
        if (np < *nplugins)
            plist[np] = p->plugin.name;

        p = p->next;
        np++;
    }

    if (*nplugins > np)
        *nplugins = np;

    return np;
}


/*
 * Tests for the presence of a specific hFILE plugin.
 *
 * Returns 1 if true
 *         0 otherwise
 */
HTSLIB_EXPORT
int hfile_has_plugin(const char *name)
{
    pthread_mutex_lock(&plugins_lock);
    if (!schemes && load_hfile_plugins() < 0) {
        pthread_mutex_unlock(&plugins_lock);
        return -1;
    }
    pthread_mutex_unlock(&plugins_lock);

    struct hFILE_plugin_list *p = plugins;
    while (p) {
        if (strcmp(p->plugin.name, name) == 0)
            return 1;
        p = p->next;
    }

    return 0;
}

/***************************
 * hFILE interface proper
 ***************************/

hFILE *hopen(const char *fname, const char *mode, ...)
{
    const struct hFILE_scheme_handler *handler = find_scheme_handler(fname);
    if (handler) {
        if (strchr(mode, ':') == NULL
            || handler->priority < 2000
            || handler->vopen == NULL) {
            return handler->open(fname, mode);
        }
        else {
            hFILE *fp;
            va_list arg;
            va_start(arg, mode);
            fp = handler->vopen(fname, mode, arg);
            va_end(arg);
            return fp;
        }
    }
    else if (strcmp(fname, "-") == 0) return hopen_fd_stdinout(mode);
    else return hopen_fd(fname, mode);
}

HTSLIB_EXPORT
int hfile_always_local (const char *fname) { return 0; }

HTSLIB_EXPORT
int hfile_always_remote(const char *fname) { return 1; }

int hisremote(const char *fname)
{
    const struct hFILE_scheme_handler *handler = find_scheme_handler(fname);
    return handler? handler->isremote(fname) : 0;
}

// Remove an extension, if any, from the basename part of [start,limit).
// Note: Doesn't notice percent-encoded '.' and '/' characters. Don't do that.
static const char *strip_extension(const char *start, const char *limit)
{
    const char *s = limit;
    while (s > start) {
        --s;
        if (*s == '.') return s;
        else if (*s == '/') break;
    }
    return limit;
}

char *haddextension(struct kstring_t *buffer, const char *filename,
                    int replace, const char *new_extension)
{
    const char *trailing, *end;

    if (find_scheme_handler(filename)) {
        // URL, so alter extensions before any trailing query or fragment parts
        // Allow # symbols in s3 URLs
        trailing = filename + ((strncmp(filename, "s3://", 5) && strncmp(filename, "s3+http://", 10) && strncmp(filename, "s3+https://", 11))  ? strcspn(filename, "?#") : strcspn(filename, "?"));
    }
    else {
        // Local path, so alter extensions at the end of the filename
        trailing = strchr(filename, '\0');
    }

    end = replace? strip_extension(filename, trailing) : trailing;

    buffer->l = 0;
    if (kputsn(filename, end - filename, buffer) >= 0 &&
        kputs(new_extension, buffer) >= 0 &&
        kputs(trailing, buffer) >= 0) return buffer->s;
    else return NULL;
}


/*
 * ----------------------------------------------------------------------
 * Minimal stub functions for knet, added after the removal of
 * hfile_net.c and knetfile.c.
 *
 * They exist purely for ABI compatibility, but are simply wrappers to
 * hFILE.  API should be compatible except knet_fileno (unused?).
 *
 * CULL THESE and knetfile.h at the next .so version bump.
 */
typedef struct knetFile_s {
    // As per htslib/knetfile.h.  Duplicated here as we don't wish to
    // have any dependence on the deprecated knetfile.h interface, plus
    // it's hopefully only temporary.
    int type, fd;
    int64_t offset;
    char *host, *port;
    int ctrl_fd, pasv_ip[4], pasv_port, max_response, no_reconnect, is_ready;
    char *response, *retr, *size_cmd;
    int64_t seek_offset;
    int64_t file_size;
    char *path, *http_host;

    // Our local addition
    hFILE *hf;
} knetFile;

HTSLIB_EXPORT
knetFile *knet_open(const char *fn, const char *mode) {
    knetFile *fp = calloc(1, sizeof(*fp));
    if (!fp) return NULL;
    if (!(fp->hf = hopen(fn, mode))) {
        free(fp);
        return NULL;
    }

    // FD backend is the only one implementing knet_fileno
    fp->fd = fp->hf->backend == &fd_backend
        ? ((hFILE_fd *)fp->hf)->fd
        : -1;

    return fp;
}

HTSLIB_EXPORT
knetFile *knet_dopen(int fd, const char *mode) {
    knetFile *fp = calloc(1, sizeof(*fp));
    if (!fp) return NULL;
    if (!(fp->hf = hdopen(fd, mode))) {
        free(fp);
        return NULL;
    }
    fp->fd = fd;
    return fp;
}

HTSLIB_EXPORT
ssize_t knet_read(knetFile *fp, void *buf, size_t len) {
    ssize_t r = hread(fp->hf, buf, len);
    fp->offset += r>0?r:0;
    return r;
}

HTSLIB_EXPORT
off_t knet_seek(knetFile *fp, off_t off, int whence) {
    off_t r = hseek(fp->hf, off, whence);
    if (r >= 0)
        fp->offset = r;
    return r;
}

HTSLIB_EXPORT
int knet_close(knetFile *fp) {
    int r = hclose(fp->hf);
    free(fp);
    return r;
}
