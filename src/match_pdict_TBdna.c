/****************************************************************************
 *                                                                          *
 *        Inexact matching of a DNA dictionary using a Trusted Band         *
 *                           Author: Herve Pages                            *
 *                                                                          *
 ****************************************************************************/
#include "Biostrings.h"

static int debug = 0;

SEXP debug_match_pdict_TBdna()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'match_pdict_TBdna.c'\n",
		debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'match_pdict_TBdna.c'\n");
#endif
	return R_NilValue;
}


/****************************************************************************
 * Inexact matching on the tails of the TBdna_PDict object
 * =======================================================
 */

/* k1 must be < k2 */
static int match_dup_tail(int k1, int k2,
		const RoSeq *dup_tail,
		const RoSeq *S,
		int max_mm, int fixedP, int fixedS, int is_count_only)
{
	int i, end1, end2, OK;
	IntBuf *match_count, *ends_buf1, *ends_buf2;

#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] ENTERING match_dup_tail()\n");
		Rprintf("[DEBUG] k1=%d k2=%d\n", k1, k2);
	}
#endif
	match_count = _MIndex_get_match_count();
	ends_buf1 = _MIndex_get_match_ends(k1);
	ends_buf2 = _MIndex_get_match_ends(k2);
	for (i = 0; i < ends_buf1->nelt; i++) {
		end1 = ends_buf1->elts[i];
		OK = _is_matching_at_Pshift(dup_tail, S, end1,
				max_mm, fixedP, fixedS);
		if (!OK)
			continue;
		if (is_count_only) {
			match_count->elts[k2]++;
			if (_MIndex_get_match_reporting_mode() == 0)
				continue;
		}
		end2 = end1 + dup_tail->nelt;
		_IntBuf_insert_at(ends_buf2, ends_buf2->nelt, end2);
	}
#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] LEAVING match_dup_tail()\n");
#endif
	return is_count_only ? match_count->elts[k2] : ends_buf2->nelt;
}

/* k1 must be < k2 */
static int match_dup_headtail(int k1, int k2,
		int pdict_W, const RoSeq *dup_head, const RoSeq *dup_tail,
		const RoSeq *S,
		int max_mm, int fixedP, int fixedS, int is_count_only)
{
	int i, end1, end2, nmismatch;
	IntBuf *match_count, *ends_buf1, *ends_buf2;

	match_count = _MIndex_get_match_count();
	ends_buf1 = _MIndex_get_match_ends(k1);
	ends_buf2 = _MIndex_get_match_ends(k2);
	for (i = 0; i < ends_buf1->nelt; i++) {
		end1 = ends_buf1->elts[i];
		nmismatch = 0;
		if (dup_head != NULL)
			nmismatch += _nmismatch_at_Pshift(dup_head, S,
					end1 - pdict_W - dup_head->nelt,
					fixedP, fixedS);
		if (dup_tail != NULL)
			nmismatch += _nmismatch_at_Pshift(dup_tail, S,
					end1,
					fixedP, fixedS);
		if (nmismatch > max_mm)
			continue;
		if (is_count_only) {
			match_count->elts[k2]++;
			if (_MIndex_get_match_reporting_mode() == 0)
				continue;
		}
		end2 = end1;
		if (dup_tail != NULL)
			end2 += dup_tail->nelt;
		_IntBuf_insert_at(ends_buf2, ends_buf2->nelt, end2);
	}
	return is_count_only ? match_count->elts[k2] : ends_buf2->nelt;
}

static int match_unq_tail(int k1,
		const RoSeq *unq_tail,
		const RoSeq *S,
		int max_mm, int fixedP, int fixedS, int is_count_only)
{
	int i, end1, OK;
	IntBuf *match_count, *ends_buf1;

#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] ENTERING match_unq_tail()\n");
		Rprintf("[DEBUG] k1=%d\n", k1);
	}
#endif
	match_count = _MIndex_get_match_count();
	ends_buf1 = _MIndex_get_match_ends(k1);
	for (i = 0; i < ends_buf1->nelt; i++) {
		end1 = ends_buf1->elts[i];
#ifdef DEBUG_BIOSTRINGS
		if (debug)
			Rprintf("[DEBUG] i=%d end1=%d\n", i, end1);
#endif
		OK = _is_matching_at_Pshift(unq_tail, S, end1,
				max_mm, fixedP, fixedS);
		if (!OK) {
#ifdef DEBUG_BIOSTRINGS
			if (debug)
				Rprintf("[DEBUG] mismatch\n");
#endif
			if (_MIndex_get_match_reporting_mode() == 0)
				continue;
			// We need to shrink the buffer we are walking on!
			// This is safe because shrinking an IntBuf object
			// should never trigger reallocation.
			_IntBuf_delete_at(ends_buf1, i--);
			continue;
		}
#ifdef DEBUG_BIOSTRINGS
		if (debug)
			Rprintf("[DEBUG] match\n");
#endif
		if (is_count_only) {
			match_count->elts[k1]++;
			if (_MIndex_get_match_reporting_mode() == 0)
				continue;
		}
		ends_buf1->elts[i] += unq_tail->nelt;
	}
#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] LEAVING match_unq_tail()\n");
#endif
	return is_count_only ? match_count->elts[k1] : ends_buf1->nelt;
}

static int match_unq_headtail(int k1, int pdict_W,
		const RoSeq *unq_head, const RoSeq *unq_tail,
		const RoSeq *S,
		int max_mm, int fixedP, int fixedS, int is_count_only)
{
	int i, end1, nmismatch;
	IntBuf *match_count, *ends_buf1;

	match_count = _MIndex_get_match_count();
	ends_buf1 = _MIndex_get_match_ends(k1);
	for (i = 0; i < ends_buf1->nelt; i++) {
		end1 = ends_buf1->elts[i];
		nmismatch = 0;
		if (unq_head != NULL)
			nmismatch += _nmismatch_at_Pshift(unq_head, S,
					end1 - pdict_W - unq_head->nelt,
					fixedP, fixedS);
		if (unq_tail != NULL)
			nmismatch += _nmismatch_at_Pshift(unq_tail, S,
					end1, fixedP, fixedS);
		if (nmismatch > max_mm) {
			if (_MIndex_get_match_reporting_mode() == 0)
				continue;
			// We need to shrink the buffer we are walking on!
			// This is safe because shrinking an IntBuf object
			// should never trigger reallocation.
			_IntBuf_delete_at(ends_buf1, i--);
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

static void match_TBdna_tail(SEXP unq2dup,
		CachedXStringSet *cached_tail,
		const RoSeq *S,
		int max_mm, int fixedP, int fixedS, int is_count_only)
{
	IntBuf *matching_keys;
	int n1, i, j, *dup, k1, k2, nmatches;
	SEXP dups;
	RoSeq Ptail;

#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] ENTERING match_TBdna_tail()\n");
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
				Ptail = _get_CachedXStringSet_elt_asRoSeq(
						cached_tail, k2);
				nmatches = match_dup_tail(k1, k2, &Ptail,
						S, max_mm, fixedP, fixedS,
						is_count_only);
				if (nmatches != 0)
					_IntBuf_insert_at(matching_keys,
							  matching_keys->nelt,
							  k2);
			}
		}
		Ptail = _get_CachedXStringSet_elt_asRoSeq(cached_tail, k1);
		nmatches = match_unq_tail(k1, &Ptail,
				S, max_mm, fixedP, fixedS, is_count_only);
		if (nmatches == 0) {
			_IntBuf_delete_at(matching_keys, i--);
			n1--;
		}
	}
#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] LEAVING match_TBdna_tail()\n");
#endif
	return;
}

static void match_TBdna_headtail(SEXP unq2dup, int pdict_W,
		CachedXStringSet *cached_head, CachedXStringSet *cached_tail,
		const RoSeq *S,
		int max_mm, int fixedP, int fixedS, int is_count_only)
{
	IntBuf *matching_keys;
	int n1, i, j, *dup, k1, k2, nmatches;
	SEXP dups;
	RoSeq Phead, Ptail;
	const RoSeq *pPhead, *pPtail;

#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] ENTERING match_TBdna_headtail()\n");
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
				pPhead = pPtail = NULL;
				if (cached_head != NULL) {
				    Phead = _get_CachedXStringSet_elt_asRoSeq(
						cached_head, k2);
				    pPhead = &Phead;
				}
				if (cached_tail != NULL) {
				    Ptail = _get_CachedXStringSet_elt_asRoSeq(
						cached_tail, k2);
				    pPtail = &Ptail;
				}
				nmatches = match_dup_headtail(k1, k2,
						pdict_W, pPhead, pPtail,
						S, max_mm, fixedP, fixedS,
						is_count_only);
				if (nmatches != 0)
					_IntBuf_insert_at(matching_keys,
							  matching_keys->nelt,
							  k2);
			}
		}
		pPhead = pPtail = NULL;
		if (cached_head != NULL) {
			Phead = _get_CachedXStringSet_elt_asRoSeq(
					cached_head, k1);
			pPhead = &Phead;
		}
		if (cached_tail != NULL) {
			Ptail = _get_CachedXStringSet_elt_asRoSeq(
					cached_tail, k1);
			pPtail = &Ptail;
		}
		nmatches = match_unq_headtail(k1,
				pdict_W, pPhead, pPtail,
				S, max_mm, fixedP, fixedS, is_count_only);
		if (nmatches == 0) {
			_IntBuf_delete_at(matching_keys, i--);
			n1--;
		}
	}
#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] LEAVING match_TBdna_headtail()\n");
#endif
	return;
}

static void match_TBdna(SEXP pdict_data, SEXP unq2dup, int pdict_W,
		CachedXStringSet *cached_head, CachedXStringSet *cached_tail,
		const RoSeq *S,
		int max_mm, int fixedP, int fixedS, int is_count_only)
{
#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] ENTERING match_TBdna()\n");
#endif
	_match_pdict_ACtree(pdict_data, S, fixedS);
	if (cached_head == NULL && cached_tail == NULL)
		_MIndex_report_matches_for_dups(unq2dup);
	else if (cached_head == NULL)
		match_TBdna_tail(unq2dup, cached_tail,
			S, max_mm, fixedP, fixedS, is_count_only);
	else
		match_TBdna_headtail(unq2dup, pdict_W,
			cached_head, cached_tail,
			S, max_mm, fixedP, fixedS, is_count_only);
#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] LEAVING match_TBdna()\n");
#endif
	return;
}


/****************************************************************************
 * .Call entry point: XString_match_pdict_TBdna
 *                    XStringViews_match_pdict_TBdna
 *
 * Arguments:
 *   'pdict_data': a 2-elt list containing pdict@actree@nodes@xp and
 *                 pdict@actree@base_codes (in this order)
 *   'pdict_length': length(pdict)
 *   'pdict_width': width(pdict)
 *   'pdict_unq2dup': pdict@dups@unq2dup
 *   'pdict_head': pdict@head (XStringSet or NULL)
 *   'pdict_tail': pdict@tail (XStringSet or NULL)
 *   'subject': subject
 *   'max_mismatch': max.mismatch (max nb of mismatches in the tail)
 *   'fixed': logical vector of length 2
 *   'count_only': TRUE or FALSE
 *   'envir': environment to be populated with the matches
 *
 * Return an R object that will be assigned to the 'ends' slot of the
 * MIndex object returned by matchPDict(). Refer to the description
 * of this slot in the matchPDict.R file for the details.
 *
 ****************************************************************************/

SEXP XString_match_pdict_TBdna(SEXP pdict_data,
		SEXP pdict_length, SEXP pdict_width, SEXP pdict_unq2dup,
		SEXP pdict_head, SEXP pdict_tail,
		SEXP subject,
		SEXP max_mismatch, SEXP fixed, SEXP count_only,
		SEXP envir)
{
	CachedXStringSet cached_head, *pcached_head, cached_tail, *pcached_tail;
	RoSeq S;
	int pdict_L, pdict_W, no_head, no_tail,
	    fixedP, fixedS, is_count_only;

#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] ENTERING XString_match_pdict_TBdna()\n");
#endif
	pdict_L = INTEGER(pdict_length)[0];
	pdict_W = INTEGER(pdict_width)[0];
	pcached_head = pcached_tail = NULL;
	no_head = pdict_head == R_NilValue;
	if (!no_head) {
		error("matchPDict() doesn't support PDict objects with "
		      "a head yet, sorry!");
		cached_head = _new_CachedXStringSet(pdict_head);
		pcached_head = &cached_head;
	}
	no_tail = pdict_tail == R_NilValue;
	if (!no_tail) {
		cached_tail = _new_CachedXStringSet(pdict_tail);
		pcached_tail = &cached_tail;
	}
	S = _get_XString_asRoSeq(subject);
	fixedP = LOGICAL(fixed)[0];
	fixedS = LOGICAL(fixed)[1];
	is_count_only = LOGICAL(count_only)[0];

	_MIndex_init_match_reporting(is_count_only, !no_head || !no_tail, pdict_L);
	match_TBdna(pdict_data, pdict_unq2dup, pdict_W,
		pcached_head, pcached_tail,
		&S, INTEGER(max_mismatch)[0], fixedP, fixedS, is_count_only);
#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] LEAVING XString_match_pdict_TBdna()\n");
#endif
	return _MIndex_reported_matches_asSEXP(envir);
}

SEXP XStringViews_match_pdict_TBdna(SEXP pdict_data,
		SEXP pdict_length, SEXP pdict_width, SEXP pdict_unq2dup,
		SEXP pdict_head, SEXP pdict_tail,
		SEXP subject, SEXP views_start, SEXP views_width,
		SEXP max_mismatch, SEXP fixed,
		SEXP count_only, SEXP envir)
{
	CachedXStringSet cached_head, *pcached_head, cached_tail, *pcached_tail;
	RoSeq S, S_view;
	int pdict_L, pdict_W, no_head, no_tail,
	    fixedP, fixedS, is_count_only,
	    nviews, i, *view_start, *view_width, view_offset;
	IntBuf global_match_count;
	IntBBuf global_match_ends;

#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] ENTERING XStringViews_match_pdict_TBdna()\n");
#endif
	pdict_L = INTEGER(pdict_length)[0];
	pdict_W = INTEGER(pdict_width)[0];
	pcached_head = pcached_tail = NULL;
	no_head = pdict_head == R_NilValue;
	if (!no_head) {
		error("matchPDict() doesn't support PDict objects with "
		      "a head yet, sorry!");
		cached_head = _new_CachedXStringSet(pdict_head);
		pcached_head = &cached_head;
	}
	no_tail = pdict_tail == R_NilValue;
	if (!no_tail) {
		cached_tail = _new_CachedXStringSet(pdict_tail);
		pcached_tail = &cached_tail;
	}
	S = _get_XString_asRoSeq(subject);
	fixedP = LOGICAL(fixed)[0];
	fixedS = LOGICAL(fixed)[1];
	is_count_only = LOGICAL(count_only)[0];

	if (is_count_only)
		global_match_count = _new_IntBuf(pdict_L, pdict_L, 0);
	else
		global_match_ends = _new_IntBBuf(pdict_L, pdict_L);
	_MIndex_init_match_reporting(is_count_only, !no_head || !no_tail, pdict_L);
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
		match_TBdna(pdict_data, pdict_unq2dup, pdict_W,
			pcached_head, pcached_tail,
			&S_view, INTEGER(max_mismatch)[0], fixedP, fixedS,
			is_count_only);
		_MIndex_merge_matches(&global_match_count, &global_match_ends, view_offset);
	}

#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] LEAVING XStringViews_match_pdict_TBdna()\n");
#endif
	if (is_count_only)
		return _IntBuf_asINTEGER(&global_match_count);
	if (envir == R_NilValue)
		return _IntBBuf_asLIST(&global_match_ends, 1);
	return _IntBBuf_toEnvir(&global_match_ends, envir, 1);
}

