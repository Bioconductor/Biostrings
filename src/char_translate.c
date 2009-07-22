#include "Biostrings.h"
#include "IRanges_interface.h"
#include <S.h> /* for Salloc() */

/*
 * --- .Call ENTRY POINT ---
 * Return an XString object, not an XStringSet object!
 */
SEXP XStringSet_char_translate(SEXP x, SEXP lkup, SEXP reverse)
{
	cachedXStringSet cached_x;
	int x_length, x_ncharsum, x_ncharmax, i, write_at;
	cachedCharSeq xx, yy;
	const char *x_xsbaseclassname;
	SEXP ans;
	char *buf;

	cached_x = _cache_XStringSet(x);
	x_length = _get_cachedXStringSet_length(&cached_x);
	x_ncharsum = x_ncharmax = 0;
	for (i = 0; i < x_length; i++) {
		xx = _get_cachedXStringSet_elt(&cached_x, i);
		x_ncharsum += xx.length;
		if (xx.length > x_ncharmax)
			x_ncharmax = xx.length;
	}
	if (x_ncharmax == 0)
		return x;
	x_xsbaseclassname = _get_XStringSet_xsbaseclassname(x);
	PROTECT(ans = _alloc_XString(x_xsbaseclassname, x_ncharsum));
	buf = Salloc((long) x_ncharmax, char);
	yy.seq = buf;
	write_at = 1;
	for (i = 0; i < x_length; i++) {
		xx = _get_cachedXStringSet_elt(&cached_x, i);
		IRanges_charcpy_from_i1i2_with_lkup(0, xx.length - 1,
			buf, xx.length,
			xx.seq, xx.length,
			INTEGER(lkup), LENGTH(lkup));
		yy.length = xx.length;
		_write_RoSeq_to_XString(ans, write_at, &yy, 0);
		write_at += yy.length;
	}
	UNPROTECT(1);
	return ans;
}

