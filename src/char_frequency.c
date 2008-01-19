#include "Biostrings.h"
#include <string.h>


/*
 * Arguments:
 *   'x_xp': x@data@xp
 *   'x_offset': x@offset
 *   'x_length': x@length
 * Returns a vector of 256 integers.
 */
SEXP char_frequency(SEXP x_xp, SEXP x_offset, SEXP x_length)
{
	SEXP ans;
	int o, l, i;
	const Rbyte *x;
	Rbyte c;

	PROTECT(ans = NEW_INTEGER(256));
	memset(INTEGER(ans), 0, 256 * sizeof(int));
	o = INTEGER(x_offset)[0];
	x = RAW(R_ExternalPtrTag(x_xp)) + o;
	l = INTEGER(x_length)[0];
	for (i = 0; i < l; i++) {
		c = *(x++);
		INTEGER(ans)[(unsigned char) c]++;
	}
	UNPROTECT(1);
	return ans;
}

SEXP multinucleotide_frequency(SEXP x_xp, SEXP x_offset, SEXP x_length,
		SEXP base_codes, SEXP width)
{
	SEXP ans;
	int ans_len, i;

	if (INTEGER(width)[0] > 8)
		error("'width' must be <= 8");
	ans_len = 1;
	for (i = 0; i < INTEGER(width)[0]; i++)
		ans_len *= 4;
	PROTECT(ans = NEW_INTEGER(ans_len));
	warning("THIS FUNCTION IS NOT READY YET");
	UNPROTECT(1);
	return ans;
}

