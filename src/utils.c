#include "Biostrings.h"


static int debug = 0;

SEXP utils_debug()
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
 Functions defined below are NOT .Call methods.
 DON'T REGISTER THEM IN init.c!
 They are all prefixed with "Biostrings_" to minimize the risk of clash with
 symbols found in libc ("_memcmp" was clashing with "_memcmp" from libc on
 churchill).
 ****************************************************************************/

/* ==========================================================================
 * Low-level memory operations.
 * --------------------------------------------------------------------------
 */

int Biostrings_memcmp(char *a, int ia, char *b, int ib, int n, size_t size)
{
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] Biostrings_memcmp(): ");
		Rprintf("a=%p ia=%d b=%p ib=%d n=%d size=%d\n",
			a, ia, b, ib, n, size);
	}
#endif
	a += ia * size;
	b += ib * size;
	/* memcmp() doesn't try to be smart by checking if a == b */
	return a == b ? 0 : memcmp(a, b, n * size);
}

/*
 * Copy memory area.
 * Reads a linear subset from 'src' defined by 'i1', 'i2'.
 * Writing is recycled in 'dest': it starts at its first byte
 * and comes back to it after it reaches its last byte.
 * IMPORTANT: Assumes that i1 <= i2.
 */
void Biostrings_memcpy_from_range(int i1, int i2,
		char *dest, size_t dest_nmemb,
		char *src, size_t src_nmemb, size_t size)
{
	register char *b;
	register int i2next, i1max, q;
	register size_t dest_size;

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

/*
 * Copy memory area.
 * Reading is recycled in 'src': it starts at its first byte
 * and comes back to it after it reaches its last byte.
 * Writes to a linear subset of 'dest' defined by 'i1', 'i2'.
 * IMPORTANT: Assumes that i1 <= i2.
 */
void Biostrings_memcpy_to_range(int i1, int i2,
		char *dest, size_t dest_nmemb,
		char *src, size_t src_nmemb, size_t size)
{
	register char *a;
	register int i2next, i1max, q;
	register size_t src_size;

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

/*
 * Copy memory area.
 * Reads from the members of 'src' that have the offsets passed in 'subset'.
 * Writing is recycled in 'dest': it starts at its first byte
 * and comes back to it after it reaches its last byte.
 */
void Biostrings_memcpy_from_subset(int *subset, int n,
		char *dest, size_t dest_nmemb,
		char *src, size_t src_nmemb, size_t size)
{
	register char *a, *b;
	register int i, j, k, z;

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

/*
 * Copy memory area.
 * Reading is recycled in 'src': it starts at its first byte
 * and comes back to it after it reaches its last byte.
 * Writes to the members of 'dest' that have the offsets passed in 'subset'.
 */
void Biostrings_memcpy_to_subset(int *subset, int n,
		char *dest, size_t dest_nmemb,
		char *src, size_t src_nmemb, size_t size)
{
	register char *a, *b;
	register int i, j, k, z;

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
