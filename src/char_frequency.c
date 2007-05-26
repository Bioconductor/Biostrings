#include "Biostrings.h"
#include <string.h>


/*
 * Arguments:
 *   'x_xp': x@data@xp
 *   'x_offset': x@offset
 *   'x_length': x@length
 * Returns a vector of 256 integers.
 */
SEXP XRaw_char_frequency(SEXP x_xp, SEXP x_offset, SEXP x_length)
{
	SEXP ans;
	int offset, length, i;
	char *x, c;

	PROTECT(ans = allocVector(INTSXP, 256));
	memset(INTEGER(ans), 0, 256 * sizeof(int));
	offset = INTEGER(x_offset)[0];
	x = CHAR(R_ExternalPtrTag(x_xp)) + offset;
	length = INTEGER(x_length)[0];
	for (i = 0; i < length; i++) {
		c = *(x++);
		INTEGER(ans)[(unsigned char) c]++;
	}
	UNPROTECT(1);
	return ans;
}
