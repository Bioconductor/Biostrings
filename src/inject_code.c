#include "Biostrings.h"
#include "XVector_interface.h"
#include "IRanges_interface.h"
#include "S4Vectors_interface.h"

/*
 * --- .Call ENTRY POINT ---
 * Return an XString object.
 */
SEXP XString_inject_code(SEXP x, SEXP start, SEXP width, SEXP code)
{
	const char *x_classname;
	Chars_holder X;
	int nranges, i, s, w;
	const int *s_p, *w_p;
	SEXP tag, ans;

	x_classname = get_classname(x);
	X = hold_XRaw(x);
	nranges = LENGTH(start); /* must be == LENGTH(width) */
	PROTECT(tag = NEW_RAW(X.length));
	memcpy(RAW(tag), X.ptr, X.length);
	for (i = 0, s_p = INTEGER(start),  w_p = INTEGER(width);
	     i < nranges;
	     i++, s_p++, w_p++)
	{
		s = *s_p;
		w = *w_p;
		if (s == NA_INTEGER || w == NA_INTEGER)
			error("Biostrings internal error in XString_inject_code():"
			      "NAs in 'start' or 'width' are not supported");
		s--; // 0-based start (offset)
		if (s < 0 || w < 0 || s + w > X.length)
			error("Biostrings internal error in XString_inject_code():"
			      "invalid start/width values");
		memset(RAW(tag) + s, INTEGER(code)[0], w);
	}
	PROTECT(ans = new_XRaw_from_tag(x_classname, tag));
	UNPROTECT(2);
	return ans;
}

