#include "Biostrings.h"

static int debug = 0;

SEXP debug_copy_seq()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'copy_seq.c'\n", debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'copy_seq.c'\n");
#endif
	return R_NilValue;
}


/* --------------------------------------------------------------------------
 * Linear copy and linear reverse copy (eventually with translation).
 */

void _copy_seq(char *dest, const char *src, size_t n, const int *chrtrtable)
{
	int i, trkey, trval;

	if (chrtrtable == NULL) {
		memcpy(dest, src, n);
	} else {
		for (i = 0; i < n; i++, dest++, src++) {
			trkey = (unsigned char) *src;
			if ((trval = chrtrtable[trkey]) == -1)
				error("sequence contains invalid code %d",
				      trkey);
			*dest = (char) trval;
		}
	}
	return;
}

void _revcopy_seq(char *dest, const char *src, size_t n, const int *chrtrtable)
{
	int i, trkey, trval;

	src += n - 1;
	if (chrtrtable == NULL) {
		for (i = 0; i < n; i++, dest++, src--)
			*dest = *src;
	} else {
		for (i = 0; i < n; i++, dest++, src--) {
			trkey = (unsigned char) *src;
			if ((trval = chrtrtable[trkey]) == -1)
				error("sequence contains invalid code %d",
				      trkey);
			*dest = (char) trval;
		}
	}
	return;
}


/* --------------------------------------------------------------------------
 * Reads a linear subset from 'src' defined by 'i1', 'i2'.
 * Writing is recycled in 'dest': it starts at its first member
 * and comes back to it after it reaches its last member.
 * Don't do anything if i1 > i2.
 *
 *   dest[(i-i1) % dest_length] <- tr(src[i]) for i1 <= i <= i2
 */
void _copy_seq_from_i1i2(int i1, int i2,
		char *dest, int dest_length,
		const char *src, int src_length,
		const int *chrtrtable)
{
	int nic; // nb of (remaining) items to copy

	if (i1 > i2)
		return;
	if (i1 < 0 || i2 >= src_length)
		error("subscript out of bounds");
	if (dest_length == 0)
		error("no destination to copy to");
	src += i1;
	nic = i2 - i1 + 1;
	while (nic >= dest_length) {
		_copy_seq(dest, src, dest_length, chrtrtable);
		src += dest_length;
		nic -= dest_length;
	}
	if (nic > 0) {
		_copy_seq(dest, src, nic, chrtrtable);
		warning("number of items to replace is not a multiple "
			"of replacement length");
	}
	return;
}


/* --------------------------------------------------------------------------
 * Writes to a linear subset of 'dest' defined by 'i1', 'i2'.
 * Reading is recycled in 'src': it starts at its first member
 * and comes back to it after it reaches its last member.
 * Don't do anything if i1 > i2.
 *
 *   dest[i] <- tr(src[(i-i1) % src_length]) for i1 <= i <= i2
 */
void _copy_seq_to_i1i2(int i1, int i2,
		char *dest, int dest_length,
		const char *src, int src_length,
		const int *chrtrtable)
{
	int nic; // nb of (remaining) items to copy

	if (i1 > i2)
		return;
	if (i1 < 0 || i2 >= dest_length)
		error("subscript out of bounds");
	if (src_length == 0)
		error("no value provided");
	dest += i1;
	nic = i2 - i1 + 1;
	while (nic >= src_length) {
		_copy_seq(dest, src, src_length, chrtrtable);
		dest += src_length;
		nic -= src_length;
	}
	if (nic > 0) {
		_copy_seq(dest, src, nic, chrtrtable);
		warning("number of items to replace is not a multiple "
			"of replacement length");
	}
	return;
}


/* --------------------------------------------------------------------------
 * Reads from the members of 'src' that have the offsets passed in 'subset'.
 * Writing is recycled in 'dest': it starts at its first member
 * and comes back to it after it reaches its last member.
 *
 *   dest[k % dest_length] <- tr(src[subset[k] - 1]) for 0 <= k < n
 */
void _copy_seq_from_subset(int *subset, int n,
		char *dest, int dest_length,
		const char *src, int src_length,
		const int *chrtrtable)
{
	int i, j, k, trkey, trval;

	if (dest_length == 0 && n != 0)
		error("no destination to copy to");
	if (chrtrtable == NULL) {
		for (k = i = 0; k < n; k++, i++) {
			j = subset[k] - 1;
			if (j < 0 || j >= src_length)
				error("subscript out of bounds");
			if (i >= dest_length)
				i = 0; /* recycle */
			dest[i] = src[j];
		}
	} else {
		for (k = i = 0; k < n; k++, i++) {
			j = subset[k] - 1;
			if (j < 0 || j >= src_length)
				error("subscript out of bounds");
			trkey = (unsigned char) src[j];
			if ((trval = chrtrtable[trkey]) == -1)
				error("sequence contains invalid code %d",
				      trkey);
			if (i >= dest_length)
				i = 0; /* recycle */
			dest[i] = (char) trval;
		}
	}
	if (i < dest_length)
		warning("number of items to replace is not a multiple "
			"of replacement length");
	return;
}


/* --------------------------------------------------------------------------
 * Writes to the members of 'dest' that have the offsets passed in 'subset'.
 * Reading is recycled in 'src': it starts at its first member
 * and comes back to it after it reaches its last member.
 *
 *   dest[subset[k] - 1] <- tr(src[k % src_length]) for 0 <= k < n
 */
void _copy_seq_to_subset(int *subset, int n,
		char *dest, int dest_length,
		const char *src, int src_length,
		const int *chrtrtable)
{
	int i, j, k, trkey, trval;

	if (src_length == 0 && n != 0)
		error("no value provided");
	if (chrtrtable == NULL) {
		for (k = j = 0; k < n; k++, j++) {
			i = subset[k] - 1;
			if (i < 0 || i >= dest_length)
				error("subscript out of bounds");
			if (j >= src_length)
				j = 0; /* recycle */
			dest[i] = src[j];
		}
	} else {
		for (k = j = 0; k < n; k++, j++) {
			i = subset[k] - 1;
			if (i < 0 || i >= dest_length)
				error("subscript out of bounds");
			trkey = (unsigned char) src[i];
			if ((trval = chrtrtable[trkey]) == -1)
				error("sequence contains invalid code %d",
				      trkey);
			if (j >= src_length)
				j = 0; /* recycle */
			dest[i] = (char) trval;
		}
	}
	if (j < src_length)
		warning("number of items to replace is not a multiple "
			"of replacement length");
	return;
}

