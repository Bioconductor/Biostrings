/****************************************************************************
 *                   Basic manipulation of MIndex objects                   *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "Biostrings.h"
#include "IRanges_interface.h"
#include "S4Vectors_interface.h"


/****************************************************************************
 * C-level slot getters.
 *
 * Be careful that these functions do NOT duplicate the returned slot.
 * Thus they cannot be made .Call() entry points!
 */

static SEXP
	width0_symbol = NULL,
	NAMES_symbol = NULL,
	dups0_symbol = NULL,
	ends_symbol = NULL;

static SEXP get_MIndex_width0(SEXP x)
{
	INIT_STATIC_SYMBOL(width0)
	return GET_SLOT(x, width0_symbol);
}

static SEXP get_MIndex_names(SEXP x)
{
	INIT_STATIC_SYMBOL(NAMES)
	return GET_SLOT(x, NAMES_symbol);
}

static SEXP get_MIndex_dups0(SEXP x)
{
	INIT_STATIC_SYMBOL(dups0)
	return GET_SLOT(x, dups0_symbol);
}

static SEXP get_MIndex_ends(SEXP x)
{
	INIT_STATIC_SYMBOL(ends)
	return GET_SLOT(x, ends_symbol);
}


/****************************************************************************
 * C-level abstract getters.
 */

MIndex_holder _hold_MIndex(SEXP x)
{
	MIndex_holder x_holder;
	SEXP dups0;

	x_holder.classname = get_classname(x);
	x_holder.width0 = get_MIndex_width0(x);
	x_holder.names = get_MIndex_names(x);
	x_holder.length = LENGTH(x_holder.width0);
	x_holder.ends = get_MIndex_ends(x);
	dups0 = get_MIndex_dups0(x);
	if (dups0 == R_NilValue) {
		x_holder.dups0_high2low = R_NilValue;
		x_holder.dups0_low2high = R_NilValue;
	} else {
		x_holder.dups0_high2low = get_H2LGrouping_high2low(dups0);
		x_holder.dups0_low2high = get_H2LGrouping_low2high(dups0);
	}
	return x_holder;
}

int _get_length_from_MIndex_holder(const MIndex_holder *x_holder)
{
	return x_holder->length;
}

int _get_width0_elt_from_MIndex_holder(const MIndex_holder *x_holder, int i)
{
	return INTEGER(x_holder->width0)[i];
}

IRanges_holder _get_elt_from_MIndex_holder(const MIndex_holder *x_holder, int i)
{
	IRanges_holder iranges_holder;
	int low;
	SEXP ends_elt;

	if (x_holder->dups0_high2low != R_NilValue
	 && LENGTH(x_holder->dups0_high2low) != 0
	 && (low = INTEGER(x_holder->dups0_high2low)[i]) != NA_INTEGER)
		i = low - 1;
	iranges_holder.classname = "IRanges";
	iranges_holder.is_constant_width = 1;
	iranges_holder.width = INTEGER(x_holder->width0) + i;
	iranges_holder.start = NULL;
	iranges_holder.SEXP_offset = 0;
	iranges_holder.names = R_NilValue;
	ends_elt = VECTOR_ELT(x_holder->ends, i);
	if (ends_elt == R_NilValue) {
		/* No need to initialize iranges_holder.end */
		iranges_holder.length = 0;
	} else {
		iranges_holder.length = LENGTH(ends_elt);
		iranges_holder.end = INTEGER(ends_elt);
	}
	return iranges_holder;
}


/****************************************************************************
 * Other MIndex utilities.
 */

/*
 * Does *inplace* addition of 'val' to all the elements in 'x' (INTSXP).
 * Never use it if 'x' is coming from the user space!
 */
static void add_val_to_INTEGER(SEXP x, int val)
{
	int i, *x_elt;

	for (i = 0, x_elt = INTEGER(x); i < LENGTH(x); i++, x_elt++)
		*x_elt += val;
	return;
}

/*
 * --- .Call ENTRY POINT ---
 * If 'x_width0' is NULL => returns the endIndex (list).
 * Otherwise 'x_width0' must be an integer vector of same length as the
 * 'x_ends' list and the startIndex is returned.
 */
SEXP ByPos_MIndex_endIndex(SEXP x_high2low, SEXP x_ends, SEXP x_width0)
{
	SEXP ans, ans_elt;
	int i, low;

	PROTECT(ans = duplicate(x_ends));
	for (i = 0; i < LENGTH(ans); i++) {
		if (x_high2low != R_NilValue
		 && LENGTH(x_high2low) != 0
		 && (low = INTEGER(x_high2low)[i]) != NA_INTEGER) {
			PROTECT(ans_elt = duplicate(VECTOR_ELT(ans, low - 1)));
			SET_ELEMENT(ans, i, ans_elt);
			UNPROTECT(1);
			continue;
		}
		if (x_width0 == R_NilValue)
			continue;
		ans_elt = VECTOR_ELT(ans, i);
		if (!IS_INTEGER(ans_elt)) // could be NULL
			continue;
		add_val_to_INTEGER(ans_elt, 1 - INTEGER(x_width0)[i]);
	}
	UNPROTECT(1);
	return ans;
}

/*
 * --- .Call ENTRY POINT ---
 */
SEXP ByPos_MIndex_combine(SEXP ends_listlist)
{
	int NTB, ans_length, i, j;
	SEXP ans, ans_elt, ends;
	IntAE *ends_buf;

	NTB = LENGTH(ends_listlist);
	if (NTB == 0)
		error("nothing to combine");
	ans_length = LENGTH(VECTOR_ELT(ends_listlist, 0));
	for (j = 1; j < NTB; j++)
		if (LENGTH(VECTOR_ELT(ends_listlist, j)) != ans_length)
			error("cannot combine MIndex objects of different lengths");
	ends_buf = new_IntAE(0, 0, 0);
	PROTECT(ans = NEW_LIST(ans_length));
	for (i = 0; i < ans_length; i++) {
		IntAE_set_nelt(ends_buf, 0);
		for (j = 0; j < NTB; j++) {
			ends = VECTOR_ELT(VECTOR_ELT(ends_listlist, j), i);
			if (ends == R_NilValue)
				continue;
			IntAE_append(ends_buf, INTEGER(ends), LENGTH(ends));
		}
		if (IntAE_get_nelt(ends_buf) == 0)
			continue;
		IntAE_qsort(ends_buf, 0, 0);
		IntAE_uniq(ends_buf, 0);
		PROTECT(ans_elt = new_INTEGER_from_IntAE(ends_buf));
		SET_ELEMENT(ans, i, ans_elt);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	return ans;
}

