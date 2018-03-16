/****************************************************************************
 *                 Basic manipulation of XStringSet objects                 *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "Biostrings.h"
#include "XVector_interface.h"
#include "IRanges_interface.h"
#include "S4Vectors_interface.h"

#include <limits.h>  /* for INT_MAX */


/****************************************************************************
 * C-level slot getters.
 *
 * Be careful that these functions do NOT duplicate the returned slot.
 * Thus they cannot be made .Call() entry points!
 */

/* Not strict "slot getters" but very much like. */

int _get_XStringSet_length(SEXP x)
{
	return get_XVectorList_length(x);
}

SEXP _get_XStringSet_width(SEXP x)
{
	return get_XVectorList_width(x);
}


/****************************************************************************
 * C-level abstract getters.
 */

XStringSet_holder _hold_XStringSet(SEXP x)
{
	return hold_XVectorList(x);
}

int _get_length_from_XStringSet_holder(const XStringSet_holder *x_holder)
{
	return get_length_from_XVectorList_holder(x_holder);
}

Chars_holder _get_elt_from_XStringSet_holder(
		const XStringSet_holder *x_holder, int i)
{
	return get_elt_from_XRawList_holder(x_holder, i);
}

XStringSet_holder _get_linear_subset_from_XStringSet_holder(
		const XStringSet_holder *x_holder, int offset, int length)
{
	return get_linear_subset_from_XVectorList_holder(x_holder,
							 offset, length);
}


/****************************************************************************
 * C-level slot setters.
 *
 */

/* WARNING: x@ranges@NAMES is modified in-place! */
void _set_XStringSet_names(SEXP x, SEXP names)
{
	set_XVectorList_names(x, names);
	return;
}


/****************************************************************************
 * _alloc_XStringSet()
 */

SEXP _alloc_XStringSet(const char *element_type, SEXP width)
{
	char classname[40];  /* longest string should be "DNAStringSet" */

	if (snprintf(classname, sizeof(classname), "%sSet", element_type) >=
	    sizeof(classname))
	{
		error("Biostrings internal error in _alloc_XStringSet(): "
		      "'classname' buffer too small");
	}
	return alloc_XRawList(classname, element_type, width);
}


/****************************************************************************
 * From CHARACTER to XStringSet and vice-versa.
 */

/* --- .Call ENTRY POINT --- */
SEXP new_XStringSet_from_CHARACTER(SEXP classname, SEXP elementType,
		SEXP x, SEXP start, SEXP width, SEXP lkup)
{
	SEXP ans, x_elt;
	XVectorList_holder ans_holder;
	const int *lkup0;
	int ans_len, lkup_len, i;
	Chars_holder ans_elt_holder;

	PROTECT(ans = alloc_XRawList(CHAR(STRING_ELT(classname, 0)),
				     CHAR(STRING_ELT(elementType, 0)),
				     width));
	ans_holder = hold_XVectorList(ans);
	ans_len = get_length_from_XVectorList_holder(&ans_holder);
	if (lkup == R_NilValue) {
		lkup0 = NULL;
	} else {
		lkup0 = INTEGER(lkup);
		lkup_len = LENGTH(lkup);
	}
	for (i = 0; i < ans_len; i++) {
		ans_elt_holder = get_elt_from_XRawList_holder(&ans_holder, i);
		x_elt = STRING_ELT(x, i);
		if (x_elt == NA_STRING) {
			UNPROTECT(1);
			error("input sequence %d is NA", i + 1);
		}
		_copy_CHARSXP_to_Chars_holder(&ans_elt_holder, x_elt,
				INTEGER(start)[i], lkup0, lkup_len);
	}
	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP new_CHARACTER_from_XStringSet(SEXP x, SEXP lkup)
{
	XStringSet_holder x_holder;
	int x_len, i;
	SEXP ans, ans_elt;
	Chars_holder x_elt_holder;

	x_holder = hold_XVectorList(x);
	x_len = get_length_from_XVectorList_holder(&x_holder);
	PROTECT(ans = NEW_CHARACTER(x_len));
	for (i = 0; i < x_len; i++) {
		x_elt_holder = get_elt_from_XRawList_holder(&x_holder, i);
		PROTECT(ans_elt = _new_CHARSXP_from_Chars_holder(
					&x_elt_holder, lkup));
		SET_STRING_ELT(ans, i, ans_elt);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * Creating a set of sequences (RoSeqs struct) from an XStringSet object.
 */

RoSeqs _new_RoSeqs_from_XStringSet(int nelt, SEXP x)
{
	RoSeqs seqs;
	XStringSet_holder x_holder;
	Chars_holder *elt1;
	int i;

	if (nelt > _get_XStringSet_length(x))
		error("_new_RoSeqs_from_XStringSet(): "
		      "'nelt' must be <= '_get_XStringSet_length(x)'");
	seqs = _alloc_RoSeqs(nelt);
	x_holder = _hold_XStringSet(x);
	for (i = 0, elt1 = seqs.elts; i < nelt; i++, elt1++)
		*elt1 = _get_elt_from_XStringSet_holder(&x_holder, i);
	return seqs;
}


/****************************************************************************
 * unlist().
 */

/* Note that XStringSet_unlist() is VERY similar to XString_xscat().
   Maybe both could be unified under a fast c() for XRaw objects. */
SEXP XStringSet_unlist(SEXP x)
{
	SEXP ans_tag, ans;
	unsigned int ans_len;
	int x_len, tag_offset, i;
	XStringSet_holder x_holder;
	Chars_holder xx;

	x_holder = _hold_XStringSet(x);
	x_len = _get_length_from_XStringSet_holder(&x_holder);

	/* 1st pass: determine 'ans_len' */
	ans_len = 0;
	for (i = 0; i < x_len; i++) {
		xx = _get_elt_from_XStringSet_holder(&x_holder, i);
		ans_len += xx.length;
		if (ans_len > INT_MAX)
			error("XStringSet object is too big to be "
			      "unlisted (would result in an XString\n"
			      "  object of length 2^31 or more)");
	}
	PROTECT(ans_tag = NEW_RAW(ans_len));

	/* 2nd pass: fill 'ans' */
	tag_offset = 0;
	for (i = 0; i < x_len; i++) {
		xx = _get_elt_from_XStringSet_holder(&x_holder, i);
		Ocopy_bytes_to_i1i2_with_lkup(tag_offset,
				tag_offset + xx.length - 1,
				(char *) RAW(ans_tag), LENGTH(ans_tag),
				xx.ptr, xx.length,
				NULL, 0);
		tag_offset += xx.length;
	}

	/* Make 'ans' */
	PROTECT(ans = new_XRaw_from_tag(get_List_elementType(x), ans_tag));
	UNPROTECT(2);
	return ans;
}

