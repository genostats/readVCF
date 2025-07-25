/// @file htslib/hfile.h
/// Buffered low-level input/output streams.
/*
    Copyright (C) 2013-2022 Genome Research Ltd.

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

#ifndef HTSLIB_HFILE_H
#define HTSLIB_HFILE_H

#include <string.h>

#include <sys/types.h>

#include "hts_defs.h"

// Ensure ssize_t exists within this header. All #includes must precede this,
// and ssize_t must be undefined again at the end of this header.
#if defined _MSC_VER && defined _INTPTR_T_DEFINED && !defined _SSIZE_T_DEFINED && !defined ssize_t
#define HTSLIB_SSIZE_T
#define ssize_t intptr_t
#endif

#ifdef __cplusplus
extern "C" {
#endif

struct hFILE_backend;
struct kstring_t;

/// Low-level input/output stream handle
/** The fields of this structure are declared here solely for the benefit
of the hFILE-related inline functions.  They may change in future releases.
User code should not use them directly; you should imagine that hFILE is an
opaque incomplete type.
*/
typedef struct hFILE {
    // @cond internal
    char *buffer, *begin, *end, *limit;
    const struct hFILE_backend *backend;
    off_t offset;
    unsigned at_eof:1, mobile:1, readonly:1;
    int has_errno;
    // @endcond
} hFILE;

/// Open the named file or URL as a stream
/** @return An hFILE pointer, or `NULL` (with _errno_ set) if an error occurred.

The usual `fopen(3)` _mode_ letters are supported: one of
`r` (read), `w` (write), `a` (append), optionally followed by any of
`+` (update), `e` (close on `exec(2)`), `x` (create exclusively),
`:` (indicates scheme-specific variable arguments follow).
*/
HTSLIB_EXPORT
hFILE *hopen(const char *filename, const char *mode, ...) HTS_RESULT_USED;

/// Associate a stream with an existing open file descriptor
/** @return An hFILE pointer, or `NULL` (with _errno_ set) if an error occurred.

Note that the file must be opened in binary mode, or else
there will be problems on platforms that make a difference
between text and binary mode.

By default, the returned hFILE "takes ownership" of the file descriptor
and _fd_ will be closed by hclose(). When _mode_ contains `S` (shared fd),
hclose() will destroy the hFILE but not close the underlying _fd_.

For socket descriptors (on Windows), _mode_ should contain `s`.
*/
HTSLIB_EXPORT
hFILE *hdopen(int fd, const char *mode) HTS_RESULT_USED;

/// Report whether the file name or URL denotes remote storage
/** @return  0 if local, 1 if remote.

"Remote" means involving e.g. explicit network access, with the implication
that callers may wish to cache such files' contents locally.
*/
HTSLIB_EXPORT
int hisremote(const char *filename) HTS_RESULT_USED;

/// Append an extension or replace an existing extension
/** @param buffer     The kstring to be used to store the modified filename
    @param filename   The filename to be (copied and) adjusted
    @param replace    If non-zero, one extension (if any) is removed first
    @param extension  The extension to be added (e.g. ".csi")
    @return  The modified filename (i.e., `buffer->s`), or NULL on error.
    @since   1.10

If _filename_ is an URL, alters extensions at the end of the `hier-part`,
leaving any trailing `?query` or `#fragment` unchanged.
*/
HTSLIB_EXPORT
char *haddextension(struct kstring_t *buffer, const char *filename,
                    int replace, const char *extension) HTS_RESULT_USED;

/// Flush (for output streams) and close the stream
/** @return  0 if successful, or `EOF` (with _errno_ set) if an error occurred.
*/
HTSLIB_EXPORT
int hclose(hFILE *fp) HTS_RESULT_USED;

/// Close the stream, without flushing or propagating errors
/** For use while cleaning up after an error only.  Preserves _errno_.
*/
HTSLIB_EXPORT
void hclose_abruptly(hFILE *fp);

/// Return the stream's error indicator
/** @return  Non-zero (in fact, an _errno_ value) if an error has occurred.

This would be called `herror()` and return true/false to parallel `ferror(3)`,
but a networking-related `herror(3)` function already exists.
*/
static inline int herrno(hFILE *fp)
{
    return fp->has_errno;
}

/// Clear the stream's error indicator
static inline void hclearerr(hFILE *fp)
{
    fp->has_errno = 0;
}

/// Reposition the read/write stream offset
/** @return  The resulting offset within the stream (as per `lseek(2)`),
    or negative if an error occurred.
*/
HTSLIB_EXPORT
off_t hseek(hFILE *fp, off_t offset, int whence) HTS_RESULT_USED;

/// Report the current stream offset
/** @return  The offset within the stream, starting from zero.
*/
static inline off_t htell(hFILE *fp)
{
    return fp->offset + (fp->begin - fp->buffer);
}

/// Read one character from the stream
/** @return  The character read, or `EOF` on end-of-file or error.
*/
static inline int hgetc(hFILE *fp)
{
    HTSLIB_EXPORT
    extern int hgetc2(hFILE *);
    return (fp->end > fp->begin)? (unsigned char) *(fp->begin++) : hgetc2(fp);
}

/// Read from the stream until the delimiter, up to a maximum length
/** @param buffer  The buffer into which bytes will be written
    @param size    The size of the buffer
    @param delim   The delimiter (interpreted as an `unsigned char`)
    @param fp      The file stream
    @return  The number of bytes read, or negative on error.
    @since   1.4

Bytes will be read into the buffer up to and including a delimiter, until
EOF is reached, or _size-1_ bytes have been written, whichever comes first.
The string will then be terminated with a NUL byte (`\0`).
*/
HTSLIB_EXPORT
ssize_t hgetdelim(char *buffer, size_t size, int delim, hFILE *fp)
    HTS_RESULT_USED;

/// Read a line from the stream, up to a maximum length
/** @param buffer  The buffer into which bytes will be written
    @param size    The size of the buffer
    @param fp      The file stream
    @return  The number of bytes read, or negative on error.
    @since   1.4

Specialization of hgetdelim() for a `\n` delimiter.
*/
static inline ssize_t HTS_RESULT_USED
hgetln(char *buffer, size_t size, hFILE *fp)
{
    return hgetdelim(buffer, size, '\n', fp);
}

/// Read a line from the stream, up to a maximum length
/** @param buffer  The buffer into which bytes will be written
    @param size    The size of the buffer (must be > 1 to be useful)
    @param fp      The file stream
    @return  _buffer_ on success, or `NULL` if an error occurred.
    @since   1.4

This function can be used as a replacement for `fgets(3)`, or together with
kstring's `kgetline()` to read arbitrarily-long lines into a _kstring_t_.
*/
HTSLIB_EXPORT
char *hgets(char *buffer, int size, hFILE *fp) HTS_RESULT_USED;

/// Peek at characters to be read without removing them from buffers
/** @param fp      The file stream
    @param buffer  The buffer to which the peeked bytes will be written
    @param nbytes  The number of bytes to peek at; limited by the size of the
                   internal buffer, which could be as small as 4K.
    @return  The number of bytes peeked, which may be less than _nbytes_
             if EOF is encountered; or negative, if there was an I/O error.

The characters peeked at remain in the stream's internal buffer, and will be
returned by later hread() etc calls.
*/
HTSLIB_EXPORT
ssize_t hpeek(hFILE *fp, void *buffer, size_t nbytes) HTS_RESULT_USED;

/// Read a block of characters from the file
/** @return  The number of bytes read, or negative if an error occurred.

The full _nbytes_ requested will be returned, except as limited by EOF
or I/O errors.
*/
static inline ssize_t HTS_RESULT_USED
hread(hFILE *fp, void *buffer, size_t nbytes)
{
    HTSLIB_EXPORT
    extern ssize_t hread2(hFILE *, void *, size_t, size_t);

    size_t n = fp->end - fp->begin;
    if (n > nbytes) n = nbytes;
    memcpy(buffer, fp->begin, n);
    fp->begin += n;
    return (n == nbytes || !fp->mobile)? (ssize_t) n : hread2(fp, buffer, nbytes, n);
}

/// Write a character to the stream
/** @return  The character written, or `EOF` if an error occurred.
*/
static inline int hputc(int c, hFILE *fp)
{
    HTSLIB_EXPORT
    extern int hputc2(int, hFILE *);
    if (fp->begin < fp->limit) *(fp->begin++) = c;
    else c = hputc2(c, fp);
    return c;
}

/// Write a string to the stream
/** @return  0 if successful, or `EOF` if an error occurred.
*/
static inline int hputs(const char *text, hFILE *fp)
{
    HTSLIB_EXPORT
    extern int hputs2(const char *, size_t, size_t, hFILE *);

    size_t nbytes = strlen(text), n = fp->limit - fp->begin;
    if (n > nbytes) n = nbytes;
    memcpy(fp->begin, text, n);
    fp->begin += n;
    return (n == nbytes)? 0 : hputs2(text, nbytes, n, fp);
}

/// Write a block of characters to the file
/** @return  Either _nbytes_, or negative if an error occurred.

In the absence of I/O errors, the full _nbytes_ will be written.
*/
static inline ssize_t HTS_RESULT_USED
hwrite(hFILE *fp, const void *buffer, size_t nbytes)
{
    HTSLIB_EXPORT
    extern ssize_t hwrite2(hFILE *, const void *, size_t, size_t);
    HTSLIB_EXPORT
    extern int hfile_set_blksize(hFILE *fp, size_t bufsiz);

    if (!fp->mobile) {
        size_t n = fp->limit - fp->begin;
        if (n < nbytes) {
            hfile_set_blksize(fp, fp->limit - fp->buffer + nbytes);
            fp->end = fp->limit;
        }
    }

    size_t n = fp->limit - fp->begin;
    if (nbytes >= n && fp->begin == fp->buffer) {
        // Go straight to hwrite2 if the buffer is empty and the request
        // won't fit.
        return hwrite2(fp, buffer, nbytes, 0);
    }

    if (n > nbytes) n = nbytes;
    memcpy(fp->begin, buffer, n);
    fp->begin += n;
    return (n==nbytes)? (ssize_t) n : hwrite2(fp, buffer, nbytes, n);
}

/// For writing streams, flush buffered output to the underlying stream
/** @return  0 if successful, or `EOF` if an error occurred.

This includes low-level flushing such as via `fdatasync(2)`.
*/
HTSLIB_EXPORT
int hflush(hFILE *fp) HTS_RESULT_USED;

/// For hfile_mem: get the internal buffer and it's size from a hfile
/** @return  buffer if successful, or NULL if an error occurred

The buffer returned should not be freed as this will happen when the
hFILE is closed.
*/
HTSLIB_EXPORT
char *hfile_mem_get_buffer(hFILE *file, size_t *length);

/// For hfile_mem: get the internal buffer and it's size from a hfile.
/** @return  buffer if successful, or NULL if an error occurred

This is similar to hfile_mem_get_buffer except that ownership of the
buffer is granted to the caller, who now has responsibility for freeing
it.  From this point onwards, the hFILE should not be used for any
purpose other than closing.
*/
HTSLIB_EXPORT
char *hfile_mem_steal_buffer(hFILE *file, size_t *length);

/// Fills out sc_list[] with the list of known URL schemes.
/**
 * @param plugin   [in]     Restricts schemes to only those from 'plugin.
 * @param sc_list  [out]    Filled out with the scheme names
 * @param nschemes [in/out] Size of sc_list (in) and number returned (out)
 *
 * Plugin may be passed in as NULL in which case all schemes are returned.
 * Use plugin "built-in" to list the built in schemes.
 * The size of sc_list is determined by the input value of *nschemes.
 * This is updated to return the output size.  It is up to the caller to
 * determine whether to call again with a larger number if this is too small.
 *
 * The return value represents the total number found matching plugin, which
 * may be larger than *nschemes if too small a value was specified.
 *
 * @return the number of schemes found on success.
 *         -1 on failure
 */
HTSLIB_EXPORT
int hfile_list_schemes(const char *plugin, const char *sc_list[], int *nschemes);

/// Fills out plist[] with the list of known hFILE plugins.
/*
 * @param plist    [out]    Filled out with the plugin names
 * @param nplugins [in/out] Size of plist (in) and number returned (out)
 *
 * The size of plist is determined by the input value of *nplugins.
 * This is updated to return the output size.  It is up to the caller to
 * determine whether to call again with a larger number if this is too small.
 *
 * The return value represents the total number found, which may be
 * larger than *nplugins if too small a value was specified.
 *
 * @return the number of plugins found on success.
 *         -1 on failure
 */
HTSLIB_EXPORT
int hfile_list_plugins(const char *plist[], int *nplugins);

/// Tests for the presence of a specific hFILE plugin.
/*
 * @param name     The name of the plugin to query.
 *
 * @return 1 if found, 0 otherwise.
 */
HTSLIB_EXPORT
int hfile_has_plugin(const char *name);

HTSLIB_EXPORT
void hfile_destroy(hFILE *);

#ifdef __cplusplus
}
#endif

#ifdef HTSLIB_SSIZE_T
#undef HTSLIB_SSIZE_T
#undef ssize_t
#endif

#endif
