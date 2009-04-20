/****************************************************************************
 *                          Fast MIndex utilities                           *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"
#include "IRanges_interface.h"

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



/****************************************************************************
 * MATCH REPORTING FACILITIES FOR PATTERN MATCHING FUNCTIONS RETURNING AN
 * MIndex OBJECT
 */

static int match_reporting_mode; // 0, 1 or 2
static int what_to_return; // 0: all matches; 1: match count; 2: matching keys
static IntAE match_count; // used when mode == 0 and initialized when mode == 2
static IntAEAE match_ends;  // used when mode >= 1
static IntAE matching_keys;

void _MIndex_init_match_reporting(int is_count_only, int with_matching_keys,
		int pdict_L)
{
	if (is_count_only == NA_LOGICAL) {
		what_to_return = 2;
		is_count_only = 1;
	} else if (is_count_only) {
		what_to_return = 1;
	} else {
		what_to_return = 0;
	}
	match_reporting_mode = is_count_only ? (with_matching_keys ? 2 : 0) : 1;
	if (match_reporting_mode == 0 || match_reporting_mode == 2)
		match_count = new_IntAE(pdict_L, pdict_L, 0);
	if (match_reporting_mode >= 1)
		match_ends = new_IntAEAE(pdict_L, pdict_L);
	matching_keys = new_IntAE(0, 0, 0);
	return;
}

void _MIndex_drop_reported_matches()
{
	int i;

	if (match_reporting_mode == 0 || match_reporting_mode == 2)
		IntAE_set_val(&match_count, 0);
	if (match_reporting_mode >= 1)
		for (i = 0; i < match_ends.nelt; i++)
			match_ends.elts[i].nelt = 0;
	matching_keys.nelt = 0;
	return;
}

int _MIndex_get_match_reporting_mode()
{
	return match_reporting_mode;
}

IntAE *_MIndex_get_match_count()
{
	return &match_count;
}

IntAE *_MIndex_get_match_ends(int key)
{
	return match_ends.elts + key;
}

IntAE *_MIndex_get_matching_keys()
{
	return &matching_keys;
}

void _MIndex_report_match(int key, int end)
{
	int is_new_matching_key;
	IntAE *ends_buf;

	if (match_reporting_mode == 0) {
		is_new_matching_key = match_count.elts[key]++ == 0;
	} else {
		ends_buf = match_ends.elts + key;
		is_new_matching_key = ends_buf->nelt == 0;
		IntAE_insert_at(ends_buf, ends_buf->nelt, end);
	}
	if (is_new_matching_key)
		IntAE_insert_at(&matching_keys,
				matching_keys.nelt, key);
	return;
}

void _MIndex_merge_matches(IntAE *global_match_count,
		const IntAEAE *global_match_ends, int view_offset)
{
	int i;
	const int *key;
	IntAE *ends_buf, *global_ends_buf;

	for (i = 0, key = matching_keys.elts;
	     i < matching_keys.nelt;
	     i++, key++)
	{
		if (match_reporting_mode == 0 || match_reporting_mode == 2) {
			global_match_count->elts[*key] += match_count.elts[*key];
			match_count.elts[*key] = 0;
		} else {
			ends_buf = match_ends.elts + *key;
			global_ends_buf = global_match_ends->elts + *key;
			IntAE_append_shifted_vals(global_ends_buf,
				ends_buf->elts, ends_buf->nelt, view_offset);
		}
		if (match_reporting_mode >= 1)
			match_ends.elts[*key].nelt = 0;
	}
	matching_keys.nelt = 0;
	return;
}

SEXP _MIndex_reported_matches_asSEXP(SEXP env)
{
	if (what_to_return == 2) {
		IntAE_sum_val(&matching_keys, 1);
		return IntAE_asINTEGER(&matching_keys);
	}
	if (what_to_return == 1)
		return IntAE_asINTEGER(&match_count);
	if (env == R_NilValue)
		return IntAEAE_asLIST(&match_ends, 1);
	return IntAEAE_toEnvir(&match_ends, env, 1);
}



/****************************************************************************
 * OTHER MIndex FACILITIES
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
 * If 'x_width' is NULL => returns the endIndex (list).
 * Otherwise 'x_width' must be an integer vector of same length as the
 * 'x_ends' list and the startIndex is returned.
 */
SEXP ByPos_MIndex_endIndex(SEXP x_dup2unq, SEXP x_ends, SEXP x_width)
{
	SEXP ans, ans_elt;
	int i, k1;

	PROTECT(ans = duplicate(x_ends));
	for (i = 0; i < LENGTH(ans); i++) {
		if (LENGTH(x_dup2unq) != 0
		 && (k1 = INTEGER(x_dup2unq)[i]) != NA_INTEGER) {
			PROTECT(ans_elt = duplicate(VECTOR_ELT(ans, k1 - 1)));
			SET_ELEMENT(ans, i, ans_elt);
			UNPROTECT(1);
			continue;
		}
		if (x_width == R_NilValue)
			continue;
		ans_elt = VECTOR_ELT(ans, i);
		if (!IS_INTEGER(ans_elt)) // could be NULL
			continue;
		add_val_to_INTEGER(ans_elt, 1 - INTEGER(x_width)[i]);
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
SEXP SparseMIndex_endIndex(SEXP x_ends_envir, SEXP x_width, SEXP x_names, SEXP all_names)
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
			if (x_width != R_NilValue)
				add_val_to_INTEGER(ans_elt, 1 - INTEGER(x_width)[j]);
			SET_ELEMENT(ans, j, ans_elt);
			UNPROTECT(1);
		}
		SET_NAMES(ans, duplicate(x_names));
		UNPROTECT(1);
	} else {
		//poffsets_order = new_IntAE(poffsets.nelt, 0, 0);
		//get_int_array_order(poffsets.elts, poffsets.nelt, poffsets_order.elts);
		//poffsets_order.nelt = poffsets.nelt; /* = poffsets_order.buflength */
		PROTECT(ans = NEW_LIST(poffsets.nelt));
		PROTECT(ans_names = NEW_CHARACTER(poffsets.nelt));
		for (i = 0; i < poffsets.nelt; i++) {
			//j = poffsets_order.elts[i];
			j = i;
			end = _get_val_from_env(STRING_ELT(symbols, j), x_ends_envir, 1);
			PROTECT(ans_elt = duplicate(end));
			if (x_width != R_NilValue)
				add_val_to_INTEGER(ans_elt, 1 - INTEGER(x_width)[i]);
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
		IntAE_qsort(&ends_buf);
		IntAE_delete_adjdups(&ends_buf);
		PROTECT(ans_elt = IntAE_asINTEGER(&ends_buf));
		SET_ELEMENT(ans, i, ans_elt);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	return ans;
}

