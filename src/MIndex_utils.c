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
static int what_to_return; // 0: all matches; 1: match count; 2: matching elts
static IntAE match_count; // used when mode == 0 and initialized when mode == 2
static IntAEAE match_ends;  // used when mode >= 1
static IntAE matching_keys;

void _MIndex_init_match_reporting(int is_count_only, int with_headtail,
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
	match_reporting_mode = is_count_only ? (with_headtail ? 2 : 0) : 1;
	if (match_reporting_mode == 0 || match_reporting_mode == 2)
		match_count = new_IntAE(pdict_L, pdict_L, 0);
	if (match_reporting_mode >= 1)
		match_ends = new_IntAEAE(pdict_L, pdict_L);
	matching_keys = new_IntAE(0, 0, 0);
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
 * Only elements in 'x_ends' that are integer vectors are shifted.
 */
SEXP ByPos_MIndex_endIndex(SEXP x_dup2unq, SEXP x_ends, SEXP shift)
{
	SEXP ans, ans_elt;
	int shift0, i, k1, j, *val;

	shift0 = INTEGER(shift)[0];
	PROTECT(ans = duplicate(x_ends));
	for (i = 0; i < LENGTH(ans); i++) {
		if (LENGTH(x_dup2unq) != 0
		 && (k1 = INTEGER(x_dup2unq)[i]) != NA_INTEGER) {
			PROTECT(ans_elt = duplicate(VECTOR_ELT(ans, k1 - 1)));
			SET_ELEMENT(ans, i, ans_elt);
			UNPROTECT(1);
			continue;
		}
		ans_elt = VECTOR_ELT(ans, i);
		if (!IS_INTEGER(ans_elt))
			continue;
		for (j = 0, val = INTEGER(ans_elt);
		     j < LENGTH(ans_elt);
		     j++, val++)
		{
			*val += shift0;
		}
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
     .Call("ByName_MIndex_endIndex", ends_envir, 0L, letters[1:10], TRUE, PACKAGE="Biostrings")
     .Call("ByName_MIndex_endIndex", ends_envir, 0L, letters[1:10], FALSE, PACKAGE="Biostrings")
 * but this doesn't:
     ends_envir[['3']] <- 33L
     .Call("ByName_MIndex_endIndex", ends_envir, 0L, letters[1:10], FALSE, PACKAGE="Biostrings")
 */
SEXP ByName_MIndex_endIndex(SEXP ends_envir, SEXP shift, SEXP names, SEXP all_names)
{
	SEXP ans, ans_elt, ans_names, symbols, end;
	int i, j;
	IntAE poffsets, poffsets_order;

	PROTECT(symbols = R_lsInternal(ends_envir, 1));
	poffsets = CHARACTER_asIntAE(symbols, -1);
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
		//poffsets_order = new_IntAE(poffsets.nelt, 0, 0);
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

