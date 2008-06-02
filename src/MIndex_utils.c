/****************************************************************************
 *                          Fast MIndex utilities                           *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"

static int debug = 0;

SEXP debug_MIndex_utils()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'MIndex_utils.c'\n",
		debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'MIndex_utils.c'\n");
#endif
	return R_NilValue;
}

/*
 * 'e1' must be an INTSXP.
 * addInt() must ALWAYS duplicate 'e1', even when e2 = 0!
 */
static SEXP addInt(SEXP e1, int e2)
{
	SEXP ans;
	int i, *val;

	PROTECT(ans = duplicate(e1));
	for (i = 0, val = INTEGER(ans); i < LENGTH(ans); i++, val++)
		*val += e2;
	UNPROTECT(1);
	return ans;
}

/*
 * --- .Call ENTRY POINT ---
 * Only elements in 'x' that are integer vectors are shifted.
 */
SEXP shiftListOfInts(SEXP x, SEXP shift)
{
	SEXP ans, ans_elt;
	int shiftval, i, j, *val;

	PROTECT(ans = duplicate(x));
	shiftval = INTEGER(shift)[0];
	for (i = 0; i < LENGTH(ans); i++) {
		ans_elt = VECTOR_ELT(ans, i);
		if (!IS_INTEGER(ans_elt))
			continue;
		for (j = 0, val = INTEGER(ans_elt); j < LENGTH(ans_elt); j++, val++)
			*val += shiftval;
	}
	UNPROTECT(1);
	return ans;
}

/*
 * --- .Call ENTRY POINT ---
 * All the keys in ends_envir must be representing integers left-padded with 0s
 * so they have the same length. This works properly:
     library(Biostrings)
     ends_envir <- new.env(parent = emptyenv())
     ends_envir[['0000000010']] <- -2:1
     ends_envir[['0000000004']] <- 9:6
     .Call("extract_endIndex", ends_envir, 0L, letters[1:10], TRUE, PACKAGE="Biostrings")
     .Call("extract_endIndex", ends_envir, 0L, letters[1:10], FALSE, PACKAGE="Biostrings")
 * but this doesn't:
     ends_envir[['3']] <- 33L
     .Call("extract_endIndex", ends_envir, 0L, letters[1:10], FALSE, PACKAGE="Biostrings")
 */
SEXP extract_endIndex(SEXP ends_envir, SEXP shift, SEXP names, SEXP all_names)
{
	SEXP ans, ans_elt, ans_names, symbols, end;
	int i, j;
	IntBuf poffsets, poffsets_order;

	PROTECT(symbols = R_lsInternal(ends_envir, 1));
	poffsets = _CHARACTER_asIntBuf(symbols, -1);
	if (LOGICAL(all_names)[0]) {
		PROTECT(ans = NEW_LIST(LENGTH(names)));
		for (i = 0; i < poffsets.nelt; i++) {
			end = _get_val_from_env(STRING_ELT(symbols, i), ends_envir, 1);
			PROTECT(ans_elt = addInt(end, INTEGER(shift)[0]));
			SET_ELEMENT(ans, poffsets.elts[i], ans_elt);
			UNPROTECT(1);
		}
		SET_NAMES(ans, duplicate(names));
		UNPROTECT(1);
	} else {
		//poffsets_order = _new_IntBuf(poffsets.nelt, 0, 0);
		//get_intorder(poffsets.nelt, poffsets.elts, poffsets_order.elts);
		//poffsets_order.nelt = poffsets.nelt; /* = poffsets_order.buflength */
		PROTECT(ans = NEW_LIST(poffsets.nelt));
		PROTECT(ans_names = NEW_CHARACTER(poffsets.nelt));
		for (i = 0; i < poffsets.nelt; i++) {
			//j = poffsets_order.elts[i];
			j = i;
			end = _get_val_from_env(STRING_ELT(symbols, j), ends_envir, 1);
			SET_ELEMENT(ans, i, addInt(end, INTEGER(shift)[0]));
			SET_STRING_ELT(ans_names, i, duplicate(STRING_ELT(names, poffsets.elts[j])));
		}
		SET_NAMES(ans, ans_names);
		UNPROTECT(2);
	}
	UNPROTECT(1);
	return ans;
}

static void add_val_to_ints(int *x, int x_nelt, int i1, int i2, int val)
{
	int i;

	if (i1 < 0)
		i1 = 0;
	if (i2 >= x_nelt)
		i2 = x_nelt - 1;
	for (i = i1; i <= i2; i++)
		x[i] += val;
	return;
}

static void add_coverages(int *x, int x_nelt, int *ends, int ends_nelt, int width, int shift)
{
	int i, *end, j_min, j_max;

	for (i = 0, end = ends; i < ends_nelt; i++, end++) {
		j_max = *end - shift;
		j_min = j_max - width + 1;
		add_val_to_ints(x, x_nelt, j_min, j_max, 1);
	}
	return;
}

/*
 * --- .Call ENTRY POINT ---
 */
SEXP ByPos_MIndex_coverage(SEXP ends_list, SEXP mindex_width, SEXP start, SEXP end)
{
	SEXP ans, ends;
	int mwidth, start0, end0, ans_length, i;

	mwidth = INTEGER(mindex_width)[0];
	start0 = INTEGER(start)[0];
	end0 = INTEGER(end)[0];
	ans_length = end0 - start0 + 1;
	PROTECT(ans = NEW_INTEGER(ans_length));
	memset(INTEGER(ans), 0, ans_length * sizeof(int));
	for (i = 0; i < LENGTH(ends_list); i++) {
		ends = VECTOR_ELT(ends_list, i);
		if (ends == R_NilValue)
			continue;
		add_coverages(INTEGER(ans), ans_length, INTEGER(ends), LENGTH(ends), mwidth, start0);
	}
	UNPROTECT(1);
	return ans;
}

/*
 * --- .Call ENTRY POINT ---
 */
SEXP ByName_MIndex_coverage(SEXP ends_envir, SEXP mindex_width, SEXP start, SEXP end)
{
	SEXP ans, symbols, ends;
	int mwidth, start0, end0, ans_length, i;

	mwidth = INTEGER(mindex_width)[0];
	start0 = INTEGER(start)[0];
	end0 = INTEGER(end)[0];
	ans_length = end0 - start0 + 1;
	PROTECT(ans = NEW_INTEGER(ans_length));
	memset(INTEGER(ans), 0, ans_length * sizeof(int));
	PROTECT(symbols = R_lsInternal(ends_envir, 1));
	for (i = 0; i < LENGTH(symbols); i++) {
		ends = _get_val_from_env(STRING_ELT(symbols, i), ends_envir, 1);
		add_coverages(INTEGER(ans), ans_length, INTEGER(ends), LENGTH(ends), mwidth, start0);
	}
	UNPROTECT(2);
	return ans;
}

