#include "Biostrings.h"
#include <S.h> /* for Salloc() */

/*
 * --- .Call ENTRY POINT ---
 * Return an XString object.
 */
SEXP inject_code(SEXP x, SEXP safe_starts, SEXP safe_widths, SEXP code)
{
	const char *x_class;
	RoSeq x_seq;
	int nranges, i, *start, *width;
	SEXP tag, xdata, ans;

	x_class = _get_class(x);
	x_seq = _get_XString_asRoSeq(x);
	nranges = LENGTH(safe_starts); /* must be == LENGTH(safe_widths) */
	PROTECT(tag = NEW_RAW(x_seq.nelt));
	memcpy(RAW(tag), x_seq.elts, x_seq.nelt);
	for (i = 0, start = INTEGER(safe_starts), width = INTEGER(safe_widths);
	     i < nranges;
	     i++, start++, width++)
	{
		memset(RAW(tag) + *start - 1, INTEGER(code)[0], *width);
	}
	PROTECT(xdata = _new_XRaw(tag));
	PROTECT(ans = _new_XString(x_class, xdata, 0, LENGTH(tag)));
	UNPROTECT(3);
	return ans;
}

