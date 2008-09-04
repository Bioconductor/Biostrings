#include "Biostrings.h"
#include "IRanges_interface.h"
#include <S.h> /* for Salloc() */

/*
 * See XRaw_copy_from_i1i2() in XRaw_utils.c for a description of the first 4 arguments.
 */
SEXP XRaw_translate_copy_from_i1i2(SEXP dest_xp, SEXP src_xp, SEXP imin, SEXP imax, SEXP lkup)
{
	SEXP dest, src;
	int i1, i2;

	dest = R_ExternalPtrTag(dest_xp);
	i1 = INTEGER(imin)[0] - 1;
	i2 = INTEGER(imax)[0] - 1;
	src = R_ExternalPtrTag(src_xp);
	IRanges_charcpy_from_i1i2_with_lkup(i1, i2,
			(char *) RAW(dest), LENGTH(dest),
			(char *) RAW(src), LENGTH(src),
			INTEGER(lkup), LENGTH(lkup));
	return dest_xp;
}

/*
 * See XRaw_copy_from_subset() in XRaw_utils.c for a description of the first 3 arguments.
 */
SEXP XRaw_translate_copy_from_subset(SEXP dest_xp, SEXP src_xp, SEXP subset, SEXP lkup)
{
	SEXP dest, src;

	dest = R_ExternalPtrTag(dest_xp);
	src = R_ExternalPtrTag(src_xp);
	IRanges_charcpy_from_subset_with_lkup(INTEGER(subset), LENGTH(subset),
			(char *) RAW(dest), LENGTH(dest),
			(char *) RAW(src), LENGTH(src),
			INTEGER(lkup), LENGTH(lkup));
	return dest_xp;
}

/*
 * See XRaw_copy_from_i1i2() in XRaw_utils.c for a description of the arguments.
 */
SEXP XRaw_reverse_copy_from_i1i2(SEXP dest_xp, SEXP src_xp, SEXP imin, SEXP imax)
{
	SEXP dest, src;
	int i1, i2;

	dest = R_ExternalPtrTag(dest_xp);
	i1 = INTEGER(imin)[0] - 1;
	i2 = INTEGER(imax)[0] - 1;
	src = R_ExternalPtrTag(src_xp);
	IRanges_reverse_memcpy_from_i1i2(i1, i2,
			(char *) RAW(dest), LENGTH(dest),
			(char *) RAW(src), LENGTH(src), sizeof(char));
	return dest_xp;
}

/*
 * See XRaw_copy_from_i1i2() in XRaw_utils.c for a description of the first 4 arguments.
 */
SEXP XRaw_reverse_translate_copy_from_i1i2(SEXP dest_xp, SEXP src_xp, SEXP imin, SEXP imax, SEXP lkup)
{
	SEXP dest, src;
	int i1, i2;

	dest = R_ExternalPtrTag(dest_xp);
	i1 = INTEGER(imin)[0] - 1;
	i2 = INTEGER(imax)[0] - 1;
	src = R_ExternalPtrTag(src_xp);
	IRanges_reverse_charcpy_from_i1i2_with_lkup(i1, i2,
			(char *) RAW(dest), LENGTH(dest),
			(char *) RAW(src), LENGTH(src),
			INTEGER(lkup), LENGTH(lkup));
	return dest_xp;
}

/*
 * --- .Call ENTRY POINT ---
 * Return an XString object, not an XStringSet object!
 */
SEXP XStringSet_char_translate(SEXP x, SEXP lkup, SEXP reverse)
{
	int x_length, x_ncharsum, x_ncharmax, i, write_at;
	CachedXStringSet cached_x;
	RoSeq xx, yy;
	const char *x_baseClass;
	SEXP ans;
	char *buf;

	x_length = _get_XStringSet_length(x);
	x_ncharsum = x_ncharmax = 0;
	cached_x = _new_CachedXStringSet(x);
	for (i = 0; i < x_length; i++) {
		xx = _get_CachedXStringSet_elt_asRoSeq(&cached_x, i);
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
		xx = _get_CachedXStringSet_elt_asRoSeq(&cached_x, i);
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

