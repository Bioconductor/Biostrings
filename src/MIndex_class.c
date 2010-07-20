/****************************************************************************
 *                   Basic manipulation of MIndex objects                   *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"
#include "IRanges_interface.h"

static int debug = 0;

SEXP debug_MIndex_class()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in file %s\n",
		debug ? "on" : "off", __FILE__);
#else
	Rprintf("Debug mode not available in file %s\n", __FILE__);
#endif
	return R_NilValue;
}


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

cachedMIndex _cache_MIndex(SEXP x)
{
	cachedMIndex cached_x;
	SEXP dups0;

	cached_x.classname = get_classname(x);
	cached_x.width0 = get_MIndex_width0(x);
	cached_x.names = get_MIndex_names(x);
	cached_x.length = LENGTH(cached_x.width0);
	cached_x.ends = get_MIndex_ends(x);
	dups0 = get_MIndex_dups0(x);
	cached_x.dups0_high2low = get_H2LGrouping_high2low(dups0);
	cached_x.dups0_low2high = get_H2LGrouping_low2high(dups0);
	return cached_x;
}

int _get_cachedMIndex_length(const cachedMIndex *cached_x)
{
	return cached_x->length;
}

int _get_cachedMIndex_elt_width0(const cachedMIndex *cached_x, int i)
{
	return INTEGER(cached_x->width0)[i];
}

cachedIRanges _get_cachedMIndex_elt(const cachedMIndex *cached_x, int i)
{
	cachedIRanges cached_iranges;
	int low;
	SEXP ends_elt;

	if (cached_x->dups0_high2low != R_NilValue
	 && LENGTH(cached_x->dups0_high2low) != 0
	 && (low = INTEGER(cached_x->dups0_high2low)[i]) != NA_INTEGER)
		i = low - 1;
	cached_iranges.classname = "IRanges";
	cached_iranges.is_constant_width = 1;
	cached_iranges.offset = 0;
	cached_iranges.width = INTEGER(cached_x->width0) + i;
	cached_iranges.start = NULL;
	cached_iranges.names = R_NilValue;
	ends_elt = VECTOR_ELT(cached_x->ends, i);
	if (ends_elt == R_NilValue) {
		/* No need to initialize cached_iranges.end */
		cached_iranges.length = 0;
	} else {
		cached_iranges.length = LENGTH(ends_elt);
		cached_iranges.end = INTEGER(ends_elt);
	}
	return cached_iranges;
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
 * All the keys in 'x_ends_envir' must be representing integers left-padded with 0s
 * so they have the same length. This works properly:
     library(Biostrings)
     ends_envir <- new.env(parent=emptyenv())
     ends_envir[['0000000010']] <- -2:1
     ends_envir[['0000000004']] <- 9:6
     .Call("SparseMIndex_endIndex",
           ends_envir, NULL, letters[1:10], TRUE,
           PACKAGE="Biostrings")
     .Call("SparseMIndex_endIndex",
           ends_envir, NULL, letters[1:10], FALSE,
           PACKAGE="Biostrings")
 * but this doesn't:
     ends_envir[['3']] <- 33L
     .Call("SparseMIndex_endIndex",
           ends_envir, NULL, letters[1:10], FALSE,
           PACKAGE="Biostrings")
 */
SEXP SparseMIndex_endIndex(SEXP x_ends_envir, SEXP x_width0, SEXP x_names, SEXP all_names)
{
	SEXP ans, ans_elt, ans_names, symbols, end;
	int i, j;
	IntAE poffsets, poffsets_order;

	PROTECT(symbols = R_lsInternal(x_ends_envir, 1));
	poffsets = CHARACTER_asIntAE(symbols, -1);
	if (LOGICAL(all_names)[0]) {
		PROTECT(ans = NEW_LIST(LENGTH(x_names)));
		for (i = 0; i < poffsets.nelt; i++) {
			j = poffsets.elts[i];
			end = _get_val_from_env(STRING_ELT(symbols, i), x_ends_envir, 1);
			PROTECT(ans_elt = duplicate(end));
			if (x_width0 != R_NilValue)
				add_val_to_INTEGER(ans_elt, 1 - INTEGER(x_width0)[j]);
			SET_ELEMENT(ans, j, ans_elt);
			UNPROTECT(1);
		}
		SET_NAMES(ans, duplicate(x_names));
		UNPROTECT(1);
	} else {
		//poffsets_order = new_IntAE(poffsets.nelt, 0, 0);
		//get_order_of_int_array(poffsets.elts, poffsets.nelt, 0, poffsets_order.elts, 0);
		//poffsets_order.nelt = poffsets.nelt; /* = poffsets_order.buflength */
		PROTECT(ans = NEW_LIST(poffsets.nelt));
		PROTECT(ans_names = NEW_CHARACTER(poffsets.nelt));
		for (i = 0; i < poffsets.nelt; i++) {
			//j = poffsets_order.elts[i];
			j = i;
			end = _get_val_from_env(STRING_ELT(symbols, j), x_ends_envir, 1);
			PROTECT(ans_elt = duplicate(end));
			if (x_width0 != R_NilValue)
				add_val_to_INTEGER(ans_elt, 1 - INTEGER(x_width0)[i]);
			SET_ELEMENT(ans, i, ans_elt);
			UNPROTECT(1);
			SET_STRING_ELT(ans_names, i, duplicate(STRING_ELT(x_names, poffsets.elts[j])));
		}
		SET_NAMES(ans, ans_names);
		UNPROTECT(2);
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
	IntAE ends_buf;

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
		ends_buf.nelt = 0;
		for (j = 0; j < NTB; j++) {
			ends = VECTOR_ELT(VECTOR_ELT(ends_listlist, j), i);
			if (ends == R_NilValue)
				continue;
			IntAE_append(&ends_buf, INTEGER(ends), LENGTH(ends));
		}
		if (ends_buf.nelt == 0)
			continue;
		IntAE_qsort(&ends_buf, 0);
		IntAE_delete_adjdups(&ends_buf);
		PROTECT(ans_elt = new_INTEGER_from_IntAE(&ends_buf));
		SET_ELEMENT(ans, i, ans_elt);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	return ans;
}

