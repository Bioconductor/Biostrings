#include "Biostrings.h"


/****************************************************************************
 * Memory allocation for a "XInteger" object
 * -----------------------------------------
 * An "XInteger" object stores its data in an "external" integer vector
 * (INTSXP vector).
 */

SEXP XInteger_alloc(SEXP ib_xp, SEXP length)
{
	SEXP tag;
	int tag_length;

	tag_length = INTEGER(length)[0];

	PROTECT(tag = allocVector(INTSXP, tag_length));
	R_SetExternalPtrTag(ib_xp, tag);
	UNPROTECT(1);
	return ib_xp;
}

SEXP XInteger_show(SEXP ib_xp)
{
	SEXP tag;
	int tag_length;

	tag = R_ExternalPtrTag(ib_xp);
	tag_length = LENGTH(tag);
	Rprintf("%d-integer buffer (starting at address %p)\n",
		tag_length, INTEGER(tag));
	return R_NilValue;
}

SEXP XInteger_length(SEXP ib_xp)
{
	SEXP tag, ans;
	int tag_length;

	tag = R_ExternalPtrTag(ib_xp);
	tag_length = LENGTH(tag);

	PROTECT(ans = allocVector(INTSXP, 1));
	INTEGER(ans)[0] = tag_length;
	UNPROTECT(1);
	return ans;
}

/*
 * From R:
 *   ib <- XInteger(30)
 *   .Call("XInteger_memcmp", ib@xp, 1:1, ib@xp, 10:10, 21:21, PACKAGE="Biostrings")
 */
SEXP XInteger_memcmp(SEXP ib1_xp, SEXP start1,
		 SEXP ib2_xp, SEXP start2, SEXP width)
{
	SEXP tag1, tag2, ans;
	int i1, i2, n;

	tag1 = R_ExternalPtrTag(ib1_xp);
	i1 = INTEGER(start1)[0] - 1;
	tag2 = R_ExternalPtrTag(ib2_xp);
	i2 = INTEGER(start2)[0] - 1;
	n = INTEGER(width)[0];

	PROTECT(ans = allocVector(INTSXP, 1));
	INTEGER(ans)[0] = Biostrings_memcmp((char *) INTEGER(tag1), i1,
					(char *) INTEGER(tag2), i2,
					n, sizeof(int));
	UNPROTECT(1);
	return ans;
}

/* ==========================================================================
 * Read/write integers to a "XInteger" object
 * --------------------------------------------------------------------------
 */

SEXP XInteger_read_ints_from_i1i2(SEXP src_xp, SEXP imin, SEXP imax)
{
	SEXP src, ans;
	int i1, i2, n;

	src = R_ExternalPtrTag(src_xp);
	i1 = INTEGER(imin)[0] - 1;
	i2 = INTEGER(imax)[0] - 1;
	n = i2 - i1 + 1;
	PROTECT(ans = allocVector(INTSXP, n));
	Biostrings_memcpy_from_i1i2(i1, i2,
			(char *) INTEGER(ans), LENGTH(ans),
			(char *) INTEGER(src), LENGTH(src), sizeof(int));
	UNPROTECT(1);
	return ans;
}

SEXP XInteger_read_ints_from_subset(SEXP src_xp, SEXP subset)
{
	SEXP src, ans;
	int n;

	src = R_ExternalPtrTag(src_xp);
	n = LENGTH(subset);
	PROTECT(ans = allocVector(INTSXP, n));
	Biostrings_memcpy_from_subset(INTEGER(subset), n,
			(char *) INTEGER(ans), n,
			(char *) INTEGER(src), LENGTH(src), sizeof(int));
	UNPROTECT(1);
	return ans;
}

/*
 * 'val' must be an integer vector.
 */
SEXP XInteger_write_ints_to_i1i2(SEXP dest_xp, SEXP imin, SEXP imax, SEXP val)
{
	SEXP dest;
	int i1, i2;

	dest = R_ExternalPtrTag(dest_xp);
	i1 = INTEGER(imin)[0] - 1;
	i2 = INTEGER(imax)[0] - 1;
	Biostrings_memcpy_to_i1i2(i1, i2,
			(char *) INTEGER(dest), LENGTH(dest),
			(char *) INTEGER(val), LENGTH(val), sizeof(int));
	return dest_xp;
}

SEXP XInteger_write_ints_to_subset(SEXP dest_xp, SEXP subset, SEXP val)
{
	SEXP dest;

	dest = R_ExternalPtrTag(dest_xp);
	Biostrings_memcpy_to_subset(INTEGER(subset), LENGTH(subset),
			(char *) INTEGER(dest), LENGTH(dest),
			(char *) INTEGER(val), LENGTH(val), sizeof(int));
	return dest_xp;
}
