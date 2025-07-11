/* The MIT License

   Copyright (C) 2011 by Attractive Chaos <attractor@live.co.uk>
   Copyright (C) 2013-2018, 2020-2021 Genome Research Ltd.

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

#define HTS_BUILDING_LIBRARY // Enables HTSLIB_EXPORT, see htslib/hts_defs.h
#include <config.h>

#include <stdarg.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include "htslib/kstring.h"

int kputd(double d, kstring_t *s) {
	int len = 0;
	char buf[21], *cp = buf+20, *ep;
	if (d == 0) {
		if (signbit(d)) {
			kputsn("-0",2,s);
			return 2;
		} else {
			kputsn("0",1,s);
			return 1;
		}
	}

	if (d < 0) {
		kputc('-',s);
		len = 1;
		d=-d;
	}
	if (!(d >= 0.0001 && d <= 999999)) {
		if (ks_resize(s, s->l + 50) < 0)
			return EOF;
		// We let stdio handle the exponent cases
		int s2 = snprintf(s->s + s->l, s->m - s->l, "%g", d);
		len += s2;
		s->l += s2;
		return len;
	}

	uint64_t i = d*10000000000LL;
	// Correction for rounding - rather ugly

	// Optimised for small numbers.
	// Better still would be __builtin_clz on hi/lo 32 and get the
	// starting point very rapidly.
	if (d<.0001)
		i+=0;
	else if (d<0.001)
		i+=5;
	else if (d < 0.01)
		i+=50;
	else if (d < 0.1)
		i+=500;
	else if (d < 1)
		i+=5000;
	else if (d < 10)
		i+=50000;
	else if (d < 100)
		i+=500000;
	else if (d < 1000)
		i+=5000000;
	else if (d < 10000)
		i+=50000000;
	else if (d < 100000)
		i+=500000000;
	else
		i+=5000000000LL;

	do {
		*--cp = '0' + i%10;
		i /= 10;
	} while (i >= 1);
	buf[20] = 0;
	int p = buf+20-cp;
	if (p <= 10) { // d < 1
		//assert(d/1);
		cp[6] = 0; ep = cp+5;// 6 precision
		while (p < 10) {
			*--cp = '0';
			p++;
		}
		*--cp = '.';
		*--cp = '0';
	} else {
		char *xp = --cp;
		while (p > 10) {
			xp[0] = xp[1];
			p--;
			xp++;
		}
		xp[0] = '.';
		cp[7] = 0; ep=cp+6;
		if (cp[6] == '.') cp[6] = 0;
	}

	// Cull trailing zeros
	while (*ep == '0' && ep > cp)
		ep--;
	char *z = ep+1;
	while (ep > cp) {
		if (*ep == '.') {
			if (z[-1] == '.')
				z[-1] = 0;
			else
				z[0] = 0;
			break;
		}
		ep--;
	}

	int sl = strlen(cp);
	len += sl;
	kputsn(cp, sl, s);
	return len;
}

int kvsprintf(kstring_t *s, const char *fmt, va_list ap)
{
	va_list args;
	int l;
	va_copy(args, ap);

	if (fmt[0] == '%' && fmt[1] == 'g' && fmt[2] == 0) {
		double d = va_arg(args, double);
		l = kputd(d, s);
		va_end(args);
		return l;
	}

	if (!s->s) {
		const size_t sz = 64;
		s->s = malloc(sz);
		if (!s->s)
			return -1;
		s->m = sz;
		s->l = 0;
	}

	l = vsnprintf(s->s + s->l, s->m - s->l, fmt, args); // This line does not work with glibc 2.0. See `man snprintf'.
	va_end(args);
	if (l + 1 > s->m - s->l) {
		if (ks_resize(s, s->l + l + 2) < 0)
			return -1;
		va_copy(args, ap);
		l = vsnprintf(s->s + s->l, s->m - s->l, fmt, args);
		va_end(args);
	}
	s->l += l;
	return l;
}

int ksprintf(kstring_t *s, const char *fmt, ...)
{
	va_list ap;
	int l;
	va_start(ap, fmt);
	l = kvsprintf(s, fmt, ap);
	va_end(ap);
	return l;
}

char *kstrtok(const char *str, const char *sep_in, ks_tokaux_t *aux)
{
	const unsigned char *p, *start, *sep = (unsigned char *) sep_in;
	if (sep) { // set up the table
		if (str == 0 && aux->finished) return 0; // no need to set up if we have finished
		aux->finished = 0;
		if (sep[0] && sep[1]) {
			aux->sep = -1;
			aux->tab[0] = aux->tab[1] = aux->tab[2] = aux->tab[3] = 0;
			for (p = sep; *p; ++p) aux->tab[*p>>6] |= 1ull<<(*p&0x3f);
		} else aux->sep = sep[0];
	}
	if (aux->finished) return 0;
	else if (str) start = (unsigned char *) str, aux->finished = 0;
	else start = (unsigned char *) aux->p + 1;
	if (aux->sep < 0) {
		for (p = start; *p; ++p)
			if (aux->tab[*p>>6]>>(*p&0x3f)&1) break;
	} else {
		// Using strchr is fast for next token, but slower for
		// last token due to extra pass from strlen.  Overall
		// on a VCF parse this func was 146% faster with // strchr.
		// Equiv to:
		// for (p = start; *p; ++p) if (*p == aux->sep) break;

		// NB: We could use strchrnul() here from glibc if detected,
		// which is ~40% faster again, but it's not so portable.
		// i.e.   p = (uint8_t *)strchrnul((char *)start, aux->sep);
		uint8_t *p2 = (uint8_t *)strchr((char *)start, aux->sep);
		p = p2 ? p2 : start + strlen((char *)start);
	}
	aux->p = (const char *) p; // end of token
	if (*p == 0) aux->finished = 1; // no more tokens
	return (char*)start;
}


int kgetline(kstring_t *s, kgets_func *fgets_fn, void *fp)
{
	size_t l0 = s->l;

	while (s->l == l0 || s->s[s->l-1] != '\n') {
		if (s->m - s->l < 200) {
			if (ks_resize(s, s->m + 200) < 0)
				return EOF;
		}
		if (fgets_fn(s->s + s->l, s->m - s->l, fp) == NULL) break;
		s->l += strlen(s->s + s->l);
	}

	if (s->l == l0) return EOF;

	if (s->l > l0 && s->s[s->l-1] == '\n') {
		s->l--;
		if (s->l > l0 && s->s[s->l-1] == '\r') s->l--;
	}
	s->s[s->l] = '\0';
	return 0;
}

int kgetline2(kstring_t *s, kgets_func2 *fgets_fn, void *fp)
{
	size_t l0 = s->l;

	while (s->l == l0 || s->s[s->l-1] != '\n') {
		if (s->m - s->l < 200) {
			// We return EOF for both EOF and error and the caller
			// needs to check for errors in fp, and we haven't
			// even got there yet.
			//
			// The only way of propagating memory errors is to
			// deliberately call something that we know triggers
			// and error so fp is also set.  This works for
			// hgets, but not for gets where reading <= 0 bytes
			// isn't an error.
			if (ks_resize(s, s->m + 200) < 0) {
				fgets_fn(s->s + s->l, 0, fp);
				return EOF;
			}
		}
		ssize_t len = fgets_fn(s->s + s->l, s->m - s->l, fp);
		if (len <= 0) break;
		s->l += len;
	}

	if (s->l == l0) return EOF;

	if (s->l > l0 && s->s[s->l-1] == '\n') {
		s->l--;
		if (s->l > l0 && s->s[s->l-1] == '\r') s->l--;
	}
	s->s[s->l] = '\0';
	return 0;
}