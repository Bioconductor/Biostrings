/****************************************************************************
 *                                                                          *
 *        Inexact matching of a DNA dictionary using a Trusted Band         *
 *                           Author: Herve Pages                            *
 *                                                                          *
 ****************************************************************************/
#include "Biostrings.h"
#include "IRanges_interface.h"

static int debug = 0;

SEXP debug_match_pdict()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'match_pdict.c'\n",
		debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'match_pdict.c'\n");
#endif
	return R_NilValue;
}

static CachedXStringSet *get_CachedXStringSet_ptr(SEXP x)
{
	CachedXStringSet *ptr;

	if (x == R_NilValue)
		return NULL;
	ptr = (CachedXStringSet *) R_alloc((long) 1, sizeof(CachedXStringSet));
	*ptr = _new_CachedXStringSet(x);
	return ptr;
}


/****************************************************************************
 * Inexact matching on the heads and tails of a TB_PDict object
 * ============================================================
 */

// Return the number of mismatches in the pattern head and tail.
static int nmismatch_in_headtail(const RoSeq *H, const RoSeq *T,
		const RoSeq *S, int Hshift, int Tshift, int max_mm)
{
	int nmismatch;

	nmismatch = _selected_nmismatch_at_Pshift_fun(H, S, Hshift, max_mm);
	if (nmismatch > max_mm)
		return nmismatch;
	max_mm -= nmismatch;
	nmismatch += _selected_nmismatch_at_Pshift_fun(T, S, Tshift, max_mm);
	return nmismatch;
}

/* k1 must be < k2 */
static int match_dup_headtail(int k1, int k2,
		int tb_width, const RoSeq *dup_head, const RoSeq *dup_tail,
		const RoSeq *S,
		int max_mm, int is_count_only)
{
	IntAE *match_count, *ends_buf1, *ends_buf2;
	int HTdeltashift, i, Tshift, nmismatch, end2;

	match_count = _MIndex_get_match_count();
	ends_buf1 = _MIndex_get_match_ends(k1);
	ends_buf2 = _MIndex_get_match_ends(k2);
	HTdeltashift = tb_width;
	if (dup_head != NULL)
		HTdeltashift += dup_head->nelt;
	for (i = 0; i < ends_buf1->nelt; i++) {
		Tshift = ends_buf1->elts[i];
		nmismatch = nmismatch_in_headtail(dup_head, dup_tail,
				S, Tshift - HTdeltashift, Tshift, max_mm);
		if (nmismatch > max_mm)
			continue;
		if (is_count_only) {
			match_count->elts[k2]++;
			if (_MIndex_get_match_reporting_mode() == 0)
				continue;
		}
		end2 = Tshift;
		if (dup_tail != NULL)
			end2 += dup_tail->nelt;
		IntAE_insert_at(ends_buf2, ends_buf2->nelt, end2);
	}
	return is_count_only ? match_count->elts[k2] : ends_buf2->nelt;
}

static int match_unq_headtail(int k1, int tb_width,
		const RoSeq *unq_head, const RoSeq *unq_tail,
		const RoSeq *S,
		int max_mm, int is_count_only)
{
	IntAE *match_count, *ends_buf1;
	int HTdeltashift, i, Tshift, nmismatch;

	match_count = _MIndex_get_match_count();
	ends_buf1 = _MIndex_get_match_ends(k1);
	HTdeltashift = tb_width;
	if (unq_head != NULL)
		HTdeltashift += unq_head->nelt;
	for (i = 0; i < ends_buf1->nelt; i++) {
		Tshift = ends_buf1->elts[i];
		nmismatch = nmismatch_in_headtail(unq_head, unq_tail,
				S, Tshift - HTdeltashift, Tshift, max_mm);
		if (nmismatch > max_mm) {
			if (_MIndex_get_match_reporting_mode() == 0)
				continue;
			// We need to shrink the buffer we are walking on!
			// This is safe because shrinking an IntAE object
			// should never trigger reallocation.
			IntAE_delete_at(ends_buf1, i--);
			continue;
		}
		if (is_count_only) {
			match_count->elts[k1]++;
			if (_MIndex_get_match_reporting_mode() == 0)
				continue;
		}
		if (unq_tail != NULL)
			ends_buf1->elts[i] += unq_tail->nelt;
	}
	return is_count_only ? match_count->elts[k1] : ends_buf1->nelt;
}

static void match_pdict_headtail(SEXP unq2dup, int tb_width,
		CachedXStringSet *cached_head, CachedXStringSet *cached_tail,
		const RoSeq *S,
		int max_mm, int is_count_only)
{
	IntAE *matching_keys;
	int n1, i, j, *dup, k1, k2, nmatches;
	SEXP dups;
	RoSeq Phead, Ptail;
	const RoSeq *H, *T;

#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] ENTERING match_pdict_headtail()\n");
#endif
	matching_keys = _MIndex_get_matching_keys();
	n1 = matching_keys->nelt;
	// The number of elements in matching_keys can increase or decrease
	// during the for loop!
	for (i = 0; i < n1; i++) {
		k1 = matching_keys->elts[i];
		dups = VECTOR_ELT(unq2dup, k1);
		if (dups != R_NilValue) {
			for (j = 0, dup = INTEGER(dups);
			     j < LENGTH(dups);
			     j++, dup++)
			{
				k2 = *dup - 1;
				H = T = NULL;
				if (cached_head != NULL) {
				    Phead = _get_CachedXStringSet_elt_asRoSeq(
						cached_head, k2);
				    H = &Phead;
				}
				if (cached_tail != NULL) {
				    Ptail = _get_CachedXStringSet_elt_asRoSeq(
						cached_tail, k2);
				    T = &Ptail;
				}
				nmatches = match_dup_headtail(k1, k2,
						tb_width, H, T,
						S, max_mm, is_count_only);
				if (nmatches != 0)
					IntAE_insert_at(matching_keys,
							matching_keys->nelt,
							k2);
			}
		}
		H = T = NULL;
		if (cached_head != NULL) {
			Phead = _get_CachedXStringSet_elt_asRoSeq(
					cached_head, k1);
			H = &Phead;
		}
		if (cached_tail != NULL) {
			Ptail = _get_CachedXStringSet_elt_asRoSeq(
					cached_tail, k1);
			T = &Ptail;
		}
		nmatches = match_unq_headtail(k1,
				tb_width, H, T,
				S, max_mm, is_count_only);
		if (nmatches == 0) {
			IntAE_delete_at(matching_keys, i--);
			n1--;
		}
	}
#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] LEAVING match_pdict_headtail()\n");
#endif
	return;
}

static void match_pdict(SEXP pptb,
		CachedXStringSet *cached_head, CachedXStringSet *cached_tail,
		const RoSeq *S,
		SEXP max_mismatch, SEXP fixed, int is_count_only)
{
	int tb_width, max_mm, fixedP, fixedS;
	SEXP unq2dup;
	const char *type;

#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] ENTERING match_pdict()\n");
#endif
	tb_width = _get_PreprocessedTB_width(pptb);
	unq2dup = _get_PreprocessedTB_unq2dup(pptb);
	type = get_classname(pptb);
	max_mm = INTEGER(max_mismatch)[0];
	fixedP = LOGICAL(fixed)[0];
	fixedS = LOGICAL(fixed)[1];
	if (strcmp(type, "Twobit") == 0)
		_match_Twobit(pptb, S, fixedS);
	else if (strcmp(type, "ACtree") == 0)
		_match_ACtree(pptb, S, fixedS);
	else
		error("%s: unsupported Trusted Band type in 'pdict'", type);
	if (cached_head != NULL || cached_tail != NULL) {
		_select_nmismatch_at_Pshift_fun(fixedP, fixedS);
		match_pdict_headtail(unq2dup, tb_width,
			cached_head, cached_tail,
			S, max_mm, is_count_only);
	}
#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] LEAVING match_pdict()\n");
#endif
	return;
}


/****************************************************************************
 * .Call entry point: XString_match_pdict
 *                    XStringViews_match_pdict
 *
 * Arguments:
 *   pptb: a PreprocessedTB object;
 *   pdict_head: head(pdict) (XStringSet or NULL);
 *   pdict_tail: tail(pdict) (XStringSet or NULL);
 *   subject: subject;
 *   max_mismatch: max.mismatch (max nb of mismatches out of the TB);
 *   fixed: logical vector of length 2;
 *   count_only: TRUE, FALSE or NA;
 *   envir: environment to be populated with the matches.
 *
 * Return an R object that will be assigned to the 'ends' slot of the
 * MIndex object returned by matchPDict(). Refer to the description
 * of this slot in the matchPDict.R file for the details.
 *
 ****************************************************************************/

SEXP XString_match_pdict(SEXP pptb, SEXP pdict_head, SEXP pdict_tail,
		SEXP subject,
		SEXP max_mismatch, SEXP fixed, SEXP count_only,
		SEXP envir)
{
	int tb_length, is_count_only;
	CachedXStringSet *cached_head, *cached_tail;
	RoSeq S;

#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] ENTERING XString_match_pdict()\n");
#endif
	tb_length = _get_PreprocessedTB_length(pptb);
	cached_head = get_CachedXStringSet_ptr(pdict_head);
	cached_tail = get_CachedXStringSet_ptr(pdict_tail);
	S = _get_XString_asRoSeq(subject);
	is_count_only = LOGICAL(count_only)[0];

	_MIndex_init_match_reporting(is_count_only,
		pdict_head != R_NilValue || pdict_tail != R_NilValue,
		tb_length);
	if (is_count_only == NA_LOGICAL)
		is_count_only = 1;
	match_pdict(pptb, cached_head, cached_tail,
		&S, max_mismatch, fixed, is_count_only);
#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] LEAVING XString_match_pdict()\n");
#endif
	return _MIndex_reported_matches_asSEXP(envir);
}

SEXP XStringViews_match_pdict(SEXP pptb, SEXP pdict_head, SEXP pdict_tail,
		SEXP subject, SEXP views_start, SEXP views_width,
		SEXP max_mismatch, SEXP fixed,
		SEXP count_only, SEXP envir)
{
	int tb_length, is_count_only;
	CachedXStringSet *cached_head, *cached_tail;
	RoSeq S, S_view;
	IntAE global_match_count;
	IntAEAE global_match_ends;
	int nviews, i, *view_start, *view_width, view_offset;

#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] ENTERING XStringViews_match_pdict()\n");
#endif
	tb_length = _get_PreprocessedTB_length(pptb);
	cached_head = get_CachedXStringSet_ptr(pdict_head);
	cached_tail = get_CachedXStringSet_ptr(pdict_tail);
	S = _get_XString_asRoSeq(subject);
	is_count_only = LOGICAL(count_only)[0];

	if (is_count_only)
		global_match_count = new_IntAE(tb_length, tb_length, 0);
	else
		global_match_ends = new_IntAEAE(tb_length, tb_length);
	_MIndex_init_match_reporting(is_count_only,
		pdict_head != R_NilValue || pdict_tail != R_NilValue,
		tb_length);
	if (is_count_only == NA_LOGICAL)
		error("Biostrings internal error in XStringViews_match_pdict(): "
		      "'count_only=NA' not supported");
	nviews = LENGTH(views_start);
	for (i = 0, view_start = INTEGER(views_start), view_width = INTEGER(views_width);
	     i < nviews;
	     i++, view_start++, view_width++)
	{
		view_offset = *view_start - 1;
		if (view_offset < 0 || view_offset + *view_width > S.nelt)
			error("'subject' has out of limits views");
		S_view.elts = S.elts + view_offset;
		S_view.nelt = *view_width;
		match_pdict(pptb, cached_head, cached_tail,
			&S_view, max_mismatch, fixed, is_count_only);
		_MIndex_merge_matches(&global_match_count, &global_match_ends, view_offset);
	}

#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] LEAVING XStringViews_match_pdict()\n");
#endif
	if (is_count_only)
		return IntAE_asINTEGER(&global_match_count);
	if (envir == R_NilValue)
		return IntAEAE_asLIST(&global_match_ends, 1);
	return IntAEAE_toEnvir(&global_match_ends, envir, 1);
}

SEXP XStringSet_vmatch_pdict(SEXP pptb, SEXP pdict_head, SEXP pdict_tail,
		SEXP subject,
		SEXP max_mismatch, SEXP fixed,
		SEXP count_only, SEXP envir)
{
	int tb_length, S_length, is_count_only;
	CachedXStringSet *cached_head, *cached_tail, S;
	SEXP ans;
	RoSeq S_elt;
	int i, *current_col;
	IntAE *match_count;

#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] ENTERING XStringSet_vmatch_pdict()\n");
#endif
	tb_length = _get_PreprocessedTB_length(pptb);
	cached_head = get_CachedXStringSet_ptr(pdict_head);
	cached_tail = get_CachedXStringSet_ptr(pdict_tail);
	S = _new_CachedXStringSet(subject);
	S_length = _get_XStringSet_length(subject);
	is_count_only = LOGICAL(count_only)[0];

	if (is_count_only)
		PROTECT(ans = allocMatrix(INTSXP, tb_length, S_length));
	else
		error("only vcountPDict() is supported for now, sorry");

	_MIndex_init_match_reporting(is_count_only,
		pdict_head != R_NilValue || pdict_tail != R_NilValue,
		tb_length);
	for (i = 0, current_col = INTEGER(ans);
	     i < S_length;
	     i++, current_col += tb_length)
	{
		S_elt = _get_CachedXStringSet_elt_asRoSeq(&S, i);
		match_pdict(pptb, cached_head, cached_tail,
			&S_elt, max_mismatch, fixed, is_count_only);
		match_count = _MIndex_get_match_count();
		memcpy(current_col, match_count->elts, sizeof(int) * match_count->nelt);
		_MIndex_drop_reported_matches();
	}

#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] LEAVING XStringSet_vmatch_pdict()\n");
#endif
	UNPROTECT(1);
	return ans;
}

