#include "Biostrings.h"
#include <S.h> /* for Srealloc() */


static int debug = 0;

SEXP Biostrings_debug_utils()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'utils.c'\n", debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'utils.c'\n");
#endif
	return R_NilValue;
}



/****************************************************************************
 Functions defined below are NOT .Call methods: they are low level routines
 used by .Call methods. They are almost "R independent" (i.e. except for the
 use of R_alloc() or the error()/warning() macros, they don't use/know
 anything about R internals).
 DON'T REGISTER THEM IN R_init_Biostrings.c!
 They are all prefixed with "_Biostrings_" to minimize the risk of clash with
 symbols found in libc (before "_memcmp" was renamed "_Biostrings_memcmp" it
 was clashing with "_memcmp" from libc on churchill).
 ****************************************************************************/


/* Alloc memory for a string of length n */
char * _Biostrings_alloc_string(int n)
{
	char * s;

	s = (char *) R_alloc((long) (n) + 1L, sizeof(char));
	s[n] = (char) 0;
	return s;
}


/* ==========================================================================
 * Memory comparison
 * --------------------------------------------------------------------------
 */

int _Biostrings_memcmp(const char *a, int ia, const char *b, int ib, int n, size_t size)
{
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _Biostrings_memcmp(): ");
		Rprintf("a=%p ia=%d b=%p ib=%d n=%d size=%d\n",
			a, ia, b, ib, n, size);
	}
#endif
	a += ia * size;
	b += ib * size;
	/* memcmp() doesn't try to be smart by checking if a == b */
	return a == b ? 0 : memcmp(a, b, n * size);
}


/* ==========================================================================
 * Memory copy:
 *   dest[(i-i1) % dest_nmemb] <- src[i] for i1 <= i <= i2
 * --------------------------------------------------------------------------
 * Reads a linear subset from 'src' defined by 'i1', 'i2'.
 * IMPORTANT: Assumes that i1 <= i2.
 * Writing is recycled in 'dest': it starts at its first member
 * and comes back to it after it reaches its last member.
 */
void _Biostrings_memcpy_from_i1i2(int i1, int i2,
		char *dest, size_t dest_nmemb,
		const char *src, size_t src_nmemb, size_t size)
{
	const char *b;
	int i2next, i1max, q;
	size_t dest_size;

	if (i1 < 0 || i2 >= src_nmemb)
		error("subscript out of bounds");
	if (dest_nmemb == 0)
		error("no destination to copy to");
	i2next = i2 + 1;
	i1max = i2next - dest_nmemb;
	b = src + i1 * size;
	dest_size = dest_nmemb * size;
	while (i1 <= i1max) {
		memcpy(dest, b, dest_size);
		b += dest_size;
		i1 += dest_nmemb;
	}
	q = i2next - i1;
	if (q > 0) {
		/* Safe because q is always < dest_nmemb */
		memcpy(dest, b, q * size);
		warning("number of items to replace is not a multiple "
			"of replacement length");
	}
	return;
}


/* ==========================================================================
 * Memory copy:
 *   dest[k % dest_nmemb] <- src[subset[k] - 1] for 0 <= k <= n
 * --------------------------------------------------------------------------
 * Reads from the members of 'src' that have the offsets passed in 'subset'.
 * Writing is recycled in 'dest': it starts at its first member
 * and comes back to it after it reaches its last member.
 */
void _Biostrings_memcpy_from_subset(int *subset, int n,
		char *dest, size_t dest_nmemb,
		const char *src, size_t src_nmemb, size_t size)
{
	char *a;
        const char *b;
	int i, j, k, z;

	if (dest_nmemb == 0 && n != 0)
		error("no destination to copy to");
	a = dest;
	for (i = k = 0; k < n; i++, k++) {
		j = subset[k] - 1;
		if (j < 0 || j >= src_nmemb)
			error("subscript out of bounds");
		if (i >= dest_nmemb) {
			i = 0; /* recycle */
			a = dest;
		}
		b = src + j * size;
		for (z = 0; z < size; z++) {
			*(a++) = *(b++);
		}
	}
	if (i != dest_nmemb)
		warning("number of items to replace is not a multiple "
			"of replacement length");
	return;
}


/* ==========================================================================
 * Memory copy:
 *   dest[i] <- src[(i-i1) % src_nmemb] for i1 <= i <= i2
 * --------------------------------------------------------------------------
 * Writes to a linear subset of 'dest' defined by 'i1', 'i2'.
 * IMPORTANT: Assumes that i1 <= i2.
 * Reading is recycled in 'src': it starts at its first member
 * and comes back to it after it reaches its last member.
 */
void _Biostrings_memcpy_to_i1i2(int i1, int i2,
		char *dest, size_t dest_nmemb,
		const char *src, size_t src_nmemb, size_t size)
{
	char *a;
	int i2next, i1max, q;
	size_t src_size;

	if (i1 < 0 || i2 >= dest_nmemb)
		error("subscript out of bounds");
	if (src_nmemb == 0)
		error("no value provided");
	i2next = i2 + 1;
	i1max = i2next - src_nmemb;
	a = dest + i1 * size;
	src_size = src_nmemb * size;
	while (i1 <= i1max) {
		memcpy(a, src, src_size);
		a += src_size;
		i1 += src_nmemb;
	}
	q = i2next - i1;
	if (q > 0) {
		/* Safe because q is always < src_nmemb */
		memcpy(a, src, q * size);
		warning("number of items to replace is not a multiple "
			"of replacement length");
	}
	return;
}


/* ==========================================================================
 * Memory copy:
 *   dest[subset[k] - 1] <- src[k % src_nmemb] for 0 <= k <= n
 * --------------------------------------------------------------------------
 * Writes to the members of 'dest' that have the offsets passed in 'subset'.
 * Reading is recycled in 'src': it starts at its first member
 * and comes back to it after it reaches its last member.
 */
void _Biostrings_memcpy_to_subset(int *subset, int n,
		char *dest, size_t dest_nmemb,
		const char *src, size_t src_nmemb, size_t size)
{
	char *a;
        const char *b;
	int i, j, k, z;

	if (src_nmemb == 0 && n != 0)
		error("no value provided");
	b = src;
	for (j = k = 0; k < n; j++, k++) {
		i = subset[k] - 1;
		if (i < 0 || i >= dest_nmemb)
			error("subscript out of bounds");
		if (j >= src_nmemb) {
			j = 0; /* recycle */
			b = src;
		}
		a = dest + i * size;
		for (z = 0; z < size; z++) {
			*(a++) = *(b++);
		}
	}
	if (j != src_nmemb)
		warning("number of items to replace is not a multiple "
			"of replacement length");
	return;
}


/* ==========================================================================
 * Memory copy with translation:
 *   dest[(i-i1) % dest_length] <- tr(src[i]) for i1 <= i <= i2
 * --------------------------------------------------------------------------
 * Reads a linear subset from 'src' defined by 'i1', 'i2'.
 * IMPORTANT: Assumes that i1 <= i2.
 * Writing is recycled in 'dest': it starts at its first member
 * and comes back to it after it reaches its last member.
 */
void _Biostrings_translate_charcpy_from_i1i2(int i1, int i2,
		char *dest, int dest_length,
		const char *src, int src_length,
		const int *lkup, int lkup_length)
{
	const char *b;
        char src_val;
	int i, j, lkup_key, lkup_val;

	if (i1 < 0 || i2 >= src_length)
		error("subscript out of bounds");
	if (dest_length == 0)
		error("no destination to copy to");
	b = src + i1;
	for (i = i1, j = 0; i <= i2; i++, j++) {
		if (j >= dest_length) { /* recycle */
			j = 0;
		}
		src_val = *(b++);
		lkup_key = (unsigned char) src_val;
		if (lkup_key >= lkup_length || (lkup_val = lkup[lkup_key]) == NA_INTEGER) {
			error("key %d not in lookup table", lkup_key);
		}
		dest[j] = (char) lkup_val;
	}
	if (j < dest_length)
		warning("number of items to replace is not a multiple "
			"of replacement length");
	return;
}


/* ==========================================================================
 * Memory copy with translation:
 *   dest[k % dest_length] <- tr(src[subset[k] - 1]) for 0 <= k <= n
 * --------------------------------------------------------------------------
 * Reads from the members of 'src' that have the offsets passed in 'subset'.
 * Writing is recycled in 'dest': it starts at its first member
 * and comes back to it after it reaches its last member.
 */
void _Biostrings_translate_charcpy_from_subset(int *subset, int n,
		char *dest, int dest_length,
		const char *src, int src_length,
		const int *lkup, int lkup_length)
{
	char src_val;
	int i, j, k, lkup_key, lkup_val;

	if (dest_length == 0 && n != 0)
		error("no destination to copy to");
	for (k = j = 0; k < n; k++, j++) {
		if (j >= dest_length) { /* recycle */
			j = 0;
		}
		i = subset[k] - 1;
		if (i < 0 || i >= src_length)
			error("subscript out of bounds");
		src_val = src[i];
		lkup_key = (unsigned char) src_val;
		if (lkup_key >= lkup_length || (lkup_val = lkup[lkup_key]) == NA_INTEGER) {
			error("key %d not in lookup table", lkup_key);
		}
		dest[j] = (char) lkup_val;
	}
	if (j < dest_length)
		warning("number of items to replace is not a multiple "
			"of replacement length");
	return;
}


/* ==========================================================================
 * Memory copy with translation:
 *   dest[i] <- tr(src[(i-i1) % src_length]) for i1 <= i <= i2
 * --------------------------------------------------------------------------
 * Writes to a linear subset of 'dest' defined by 'i1', 'i2'.
 * IMPORTANT: Assumes that i1 <= i2.
 * Reading is recycled in 'src': it starts at its first member
 * and comes back to it after it reaches its last member.
 */
void _Biostrings_translate_charcpy_to_i1i2(int i1, int i2,
		char *dest, int dest_length,
		const char *src, int src_length,
		const int *lkup, int lkup_length)
{
	char *a, src_val;
	int i, j, lkup_key, lkup_val;

	if (i1 < 0 || i2 >= dest_length)
		error("subscript out of bounds");
	if (src_length == 0)
		error("no value provided");
	a = dest + i1;
	for (i = i1, j = 0; i <= i2; i++, j++) {
		if (j >= src_length) { /* recycle */
			j = 0;
		}
		src_val = src[j];
		lkup_key = (unsigned char) src_val;
		if (lkup_key >= lkup_length || (lkup_val = lkup[lkup_key]) == NA_INTEGER) {
			error("key %d not in lookup table", lkup_key);
		}
		*(a++) = (char) lkup_val;
	}
	if (j < src_length)
		warning("number of items to replace is not a multiple "
			"of replacement length");
	return;
}


/* ==========================================================================
 * Memory copy with translation:
 *   dest[subset[k] - 1] <- tr(src[k % src_length]) for 0 <= k <= n
 * --------------------------------------------------------------------------
 * Writes to the members of 'dest' that have the offsets passed in 'subset'.
 * Reading is recycled in 'src': it starts at its first member
 * and comes back to it after it reaches its last member.
 */
void _Biostrings_translate_charcpy_to_subset(int *subset, int n,
		char *dest, int dest_length,
		const char *src, int src_length,
		const int *lkup, int lkup_length)
{
	char src_val;
	int i, j, k, lkup_key, lkup_val;

	if (src_length == 0 && n != 0)
		error("no value provided");
	for (k = j = 0; k < n; k++, j++) {
		if (j >= src_length) { /* recycle */
			j = 0;
		}
		i = subset[k] - 1;
		if (i < 0 || i >= dest_length)
			error("subscript out of bounds");
		src_val = src[j];
		lkup_key = (unsigned char) src_val;
		if (lkup_key >= lkup_length || (lkup_val = lkup[lkup_key]) == NA_INTEGER) {
			error("key %d not in lookup table", lkup_key);
		}
		dest[i] = (char) lkup_val;
	}
	if (j < src_length)
		warning("number of items to replace is not a multiple "
			"of replacement length");
	return;
}


/* ==========================================================================
 * Memory copy with reverse order:
 *   dest[(dest_nmemb-1-(i-i1)) % dest_nmemb] <- src[i] for i1 <= i <= i2
 * --------------------------------------------------------------------------
 * Reads a linear subset from 'src' defined by 'i1', 'i2'.
 * IMPORTANT: Assumes that i1 <= i2.
 * Writing is recycled in 'dest': it starts at its last member
 * and comes back to it after it reaches its first member.
 */
void _Biostrings_reverse_memcpy_from_i1i2(int i1, int i2,
		char *dest, size_t dest_nmemb,
		const char *src, size_t src_nmemb, size_t size)
{
	char *a;
        const char *b;
	int i, j, z;

	if (i1 < 0 || i2 >= src_nmemb)
		error("subscript out of bounds");
	if (dest_nmemb == 0)
		error("no destination to copy to");
	b = src + i1 * size;
	for (i = i1, j = dest_nmemb - 1; i <= i2; i++, j--) {
		if (j < 0) { /* recycle */
			j = dest_nmemb - 1;
		}
		a = dest + j * size;
		for (z = 0; z < size; z++) {
			*(a++) = *(b++);
		}
	}
	if (j >= 0)
		warning("number of items to replace is not a multiple "
			"of replacement length");
	return;
}


/* ==========================================================================
 * Memory copy with reverse order and translation:
 *   dest[(dest_length-1-(i-i1)) % dest_length] <- tr(src[i]) for i1 <= i <= i2
 * --------------------------------------------------------------------------
 * Reads a linear subset from 'src' defined by 'i1', 'i2'.
 * IMPORTANT: Assumes that i1 <= i2.
 * Writing is recycled in 'dest': it starts at its last member
 * and comes back to it after it reaches its first member.
 */
void _Biostrings_reverse_translate_charcpy_from_i1i2(int i1, int i2,
		char *dest, int dest_length,
		const char *src, int src_length,
		const int *lkup, int lkup_length)
{
	const char *b;
        char src_val;
	int i, j, lkup_key, lkup_val;

	if (i1 < 0 || i2 >= src_length)
		error("subscript out of bounds");
	if (dest_length == 0)
		error("no destination to copy to");
	b = src + i1;
	for (i = i1, j = dest_length - 1; i <= i2; i++, j--) {
		if (j < 0) { /* recycle */
			j = dest_length - 1;
		}
		src_val = *(b++);
		lkup_key = (unsigned char) src_val;
		if (lkup_key >= lkup_length || (lkup_val = lkup[lkup_key]) == NA_INTEGER) {
			error("key %d not in lookup table", lkup_key);
		}
		dest[j] = (char) lkup_val;
	}
	if (j >= 0)
		warning("number of items to replace is not a multiple "
			"of replacement length");
	return;
}


/* ==========================================================================
 * Memory copy with conversion to complex values:
 *   dest[(i-i1) % dest_length] <- toComplex(src[i]) for i1 <= i <= i2
 * --------------------------------------------------------------------------
 * Reads a linear subset from 'src' defined by 'i1', 'i2'.
 * IMPORTANT: Assumes that i1 <= i2.
 * Writing is recycled in 'dest': it starts at its first member
 * and comes back to it after it reaches its last member.
 */
void _Biostrings_coerce_to_complex_from_i1i2(int i1, int i2,
		Rcomplex *dest, int dest_length,
		const char *src, int src_length,
		const Rcomplex *lkup, int lkup_length)
{
	const char *b;
        char src_val;
	int i, j, lkup_key;
	Rcomplex lkup_val;

	if (i1 < 0 || i2 >= src_length)
		error("subscript out of bounds");
	if (dest_length == 0)
		error("no destination to copy to");
	b = src + i1;
	for (i = i1, j = 0; i <= i2; i++, j++) {
		if (j >= dest_length) { /* recycle */
			j = 0;
		}
		src_val = *(b++);
		lkup_key = (unsigned char) src_val;
		if (lkup_key >= lkup_length
		 || ISNA((lkup_val = lkup[lkup_key]).r)
		 || ISNA(lkup_val.i)) {
			error("key %d not in lookup table", lkup_key);
		}
		dest[j] = lkup_val;
	}
	if (j < dest_length)
		warning("number of items to replace is not a multiple "
			"of replacement length");
	return;
}


/* ==========================================================================
 * Helper functions used for storing views in a temporary buffer like the
 * matches found by the matching algos.
 * --------------------------------------------------------------------------
 */

static int *views_start, *views_end;
static int views_buffer_size, views_count;

/* Reset views buffer */
void _Biostrings_reset_views_buffer()
{
	/* No memory leak here, because we use transient storage allocation */
	views_start = views_end = NULL;
	views_buffer_size = views_count = 0;
	return;
}

/* Return the new number of views */
int _Biostrings_report_view(int start, int end)
{
	long new_size;

	if (views_count >= views_buffer_size) {
		/* Buffer is full */
		if (views_buffer_size == 0)
			new_size = 1024;
		else
			new_size = 2 * views_buffer_size;
		views_start = Srealloc((char *) views_start, new_size,
						(long) views_buffer_size, int);
		views_end = Srealloc((char *) views_end, new_size,
						(long) views_buffer_size, int);
		views_buffer_size = new_size;
	}
	views_start[views_count] = start;
	views_end[views_count] = end;
	return ++views_count;
}

/* Return the new number of views */
int _Biostrings_report_match(int Lpos, int Rpos)
{
	return _Biostrings_report_view(++Lpos, ++Rpos);
}

int *_Biostrings_get_views_start()
{
	return views_start;
}

int *_Biostrings_get_views_end()
{
	return views_end;
}

