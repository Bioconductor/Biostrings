#include "Biostrings.h"
#include "IRanges_interface.h"
#include <S.h> /* for Salloc() */

/*
 * --- .Call ENTRY POINT ---
 * Return an XString object, not an XStringSet object!
 */
SEXP XStringSet_char_translate(SEXP x, SEXP lkup, SEXP reverse)
{
	int x_length, x_ncharsum, x_ncharmax, i, write_at;
	cachedXStringSet cached_x;
	RoSeq xx, yy;
	const char *x_baseClass;
	SEXP ans;
	char *buf;

	x_length = _get_XStringSet_length(x);
	x_ncharsum = x_ncharmax = 0;
	cached_x = _cache_XStringSet(x);
	for (i = 0; i < x_length; i++) {
		xx = _get_cachedXStringSet_elt(&cached_x, i);
		x_ncharsum += xx.nelt;
		if (xx.nelt > x_ncharmax)
			x_ncharmax = xx.nelt;
	}
	if (x_ncharmax == 0)
		return x;
	x_baseClass = _get_XStringSet_baseClass(x);
	PROTECT(ans = _alloc_XString(x_baseClass, x_ncharsum));
	buf = Salloc((long) x_ncharmax, char);
	yy.elts = buf;
	write_at = 1;
	for (i = 0; i < x_length; i++) {
		xx = _get_cachedXStringSet_elt(&cached_x, i);
		IRanges_charcpy_from_i1i2_with_lkup(0, xx.nelt - 1,
			buf, xx.nelt,
			xx.elts, xx.nelt,
			INTEGER(lkup), LENGTH(lkup));
		yy.nelt = xx.nelt;
		_write_RoSeq_to_XString(ans, write_at, &yy, 0);
		write_at += yy.nelt;
	}
	UNPROTECT(1);
	return ans;
}

