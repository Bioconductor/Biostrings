#include "Biostrings.h"


/****************************************************************************
 * Memory allocation for a "ibuf" object
 * -------------------------------------
 * An "ibuf" object stores its data in an "external" integer vector
 * (INTSXP vector).
 */

SEXP ibuf_alloc(SEXP ib_xp, SEXP length)
{
	SEXP tag;
	int tag_length;

	tag_length = INTEGER(length)[0];

	PROTECT(tag = allocVector(INTSXP, tag_length));
	R_SetExternalPtrTag(ib_xp, tag);
	UNPROTECT(1);
	return ib_xp;
}

SEXP ibuf_show(SEXP ib_xp)
{
	SEXP tag;
	int tag_length;

	tag = R_ExternalPtrTag(ib_xp);
	tag_length = LENGTH(tag);
	Rprintf("%d-integer buffer (starting at address %p)\n",
		tag_length, INTEGER(tag));
	return R_NilValue;
}

SEXP ibuf_length(SEXP ib_xp)
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
 *   ib <- ibuf(30)
 *   .Call("ibuf_memcmp", ib@xp, 1:1, ib@xp, 10:10, 21:21, PACKAGE="Biostrings")
 */
SEXP ibuf_memcmp(SEXP ib1_xp, SEXP first1,
		 SEXP ib2_xp, SEXP first2, SEXP width)
{
	SEXP tag1, tag2, ans;
	int i1, i2, n;

	tag1 = R_ExternalPtrTag(ib1_xp);
	i1 = INTEGER(first1)[0] - 1;
	tag2 = R_ExternalPtrTag(ib2_xp);
	i2 = INTEGER(first2)[0] - 1;
	n = INTEGER(width)[0];

	PROTECT(ans = allocVector(INTSXP, 1));
	INTEGER(ans)[0] = _memcmp((char *) INTEGER(tag1), i1,
				  (char *) INTEGER(tag2), i2,
				  n, sizeof(int));
	UNPROTECT(1);
	return ans;
}

/* ==========================================================================
 * Read/write integers to a "ibuf" object
 * --------------------------------------------------------------------------
 */

SEXP ibuf_read_ints(SEXP ib_xp, SEXP imin, SEXP imax)
{
	SEXP tag, ans;
	int i1, i2, n;

	tag = R_ExternalPtrTag(ib_xp);
	i1 = INTEGER(imin)[0] - 1;
	i2 = INTEGER(imax)[0] - 1;
	n = i2 - i1 + 1;

	PROTECT(ans = allocVector(INTSXP, n));
	_memcpy_from_range(i1, i2,
			(char *) INTEGER(ans), LENGTH(ans),
			(char *) INTEGER(tag), LENGTH(tag), sizeof(int));
	Rprintf("Il est passe par ici\n");
	printf("Il repassera par la\n");
	UNPROTECT(1);
	return ans;
}

SEXP ibuf_readii_ints(SEXP ib_xp, SEXP ii)
{
	SEXP tag, ans;	
	int n;

	tag = R_ExternalPtrTag(ib_xp);
	n = LENGTH(ii);

	PROTECT(ans = allocVector(INTSXP, n));
	_memcpy_from_subset(INTEGER(ii), n,
			(char *) INTEGER(ans), n,
			(char *) INTEGER(tag), LENGTH(tag), sizeof(int));
	UNPROTECT(1);
	return ans;
}

/*
 * 'val' must be an integer vector.
 */
SEXP ibuf_write_ints(SEXP ib_xp, SEXP imin, SEXP imax, SEXP val)
{
	SEXP tag;
	int i1, i2;

	tag = R_ExternalPtrTag(ib_xp);
	i1 = INTEGER(imin)[0] - 1;
	i2 = INTEGER(imax)[0] - 1;
	_memcpy_to_range(i1, i2,
			(char *) INTEGER(tag), LENGTH(tag),
			(char *) INTEGER(val), LENGTH(val), sizeof(int));
	return ib_xp;
}

SEXP ibuf_writeii_ints(SEXP ib_xp, SEXP ii, SEXP val)
{
	SEXP tag;

	tag = R_ExternalPtrTag(ib_xp);
	_memcpy_to_subset(INTEGER(ii), LENGTH(ii),
			(char *) INTEGER(tag), LENGTH(tag),
			(char *) INTEGER(val), LENGTH(val), sizeof(int));
	return ib_xp;
}
