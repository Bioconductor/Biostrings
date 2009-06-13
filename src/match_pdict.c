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

static int is_count_only(SEXP matches_as)
{
	const char *matches_as0;

	matches_as0 = CHAR(STRING_ELT(matches_as, 0));
	if (strcmp(matches_as0, "MATCHES_AS_WHICH") == 0)
		return 1;
	if (strcmp(matches_as0, "MATCHES_AS_COUNTS") == 0)
		return 1;
	return 0;
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
		const RoSeq *S, int max_mm,
		int count_only, Seq2MatchBuf *seq2match_buf)
{
	IntAE *counts_buf, *ends_buf1, *ends_buf2;
	int HTdeltashift, i, Tshift, nmismatch, end2;

	counts_buf = &(seq2match_buf->match_counts);
	ends_buf1 = seq2match_buf->match_ends.elts + k1;
	ends_buf2 = seq2match_buf->match_ends.elts + k2;
	HTdeltashift = tb_width;
	if (dup_head != NULL)
		HTdeltashift += dup_head->nelt;
	for (i = 0; i < ends_buf1->nelt; i++) {
		Tshift = ends_buf1->elts[i];
		nmismatch = nmismatch_in_headtail(dup_head, dup_tail,
				S, Tshift - HTdeltashift, Tshift, max_mm);
		if (nmismatch > max_mm)
			continue;
		counts_buf->elts[k2]++;
		if (count_only)
			continue;
		end2 = Tshift;
		if (dup_tail != NULL)
			end2 += dup_tail->nelt;
		IntAE_insert_at(ends_buf2, ends_buf2->nelt, end2);
	}
	return counts_buf->elts[k2];
}

static int match_unq_headtail(int k1, int tb_width,
		const RoSeq *unq_head, const RoSeq *unq_tail,
		const RoSeq *S, int max_mm,
		int count_only, Seq2MatchBuf *seq2match_buf)
{
	IntAE *counts_buf, *ends_buf1;
	int HTdeltashift, i, Tshift, nmismatch;

	counts_buf = &(seq2match_buf->match_counts);
	ends_buf1 = seq2match_buf->match_ends.elts + k1;
	HTdeltashift = tb_width;
	if (unq_head != NULL)
		HTdeltashift += unq_head->nelt;
	for (i = 0; i < ends_buf1->nelt; i++) {
		Tshift = ends_buf1->elts[i];
		nmismatch = nmismatch_in_headtail(unq_head, unq_tail,
				S, Tshift - HTdeltashift, Tshift, max_mm);
		if (nmismatch > max_mm) {
			if (count_only)
				continue;
			// We need to shrink the buffer we are walking on!
			// This is safe because shrinking an IntAE object
			// should never trigger reallocation.
			IntAE_delete_at(ends_buf1, i--);
			continue;
		}
		counts_buf->elts[k1]++;
		if (count_only)
			continue;
		if (unq_tail != NULL)
			ends_buf1->elts[i] += unq_tail->nelt;
	}
	return counts_buf->elts[k1];
}

/* If 'cached_head' and 'cached_tail' are NULL then match_headtail() just
   propagates the matches to the duplicates */
static void match_headtail(SEXP low2high, int tb_width,
		CachedXStringSet *cached_head, CachedXStringSet *cached_tail,
		const RoSeq *S, int max_mm,
		int count_only, Seq2MatchBuf *seq2match_buf)
{
	IntAE *matching_keys;
	int n1, i, j, *dup, k1, k2, nmatches;
	SEXP dups;
	RoSeq Phead, Ptail;
	const RoSeq *H, *T;

#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] ENTERING match_headtail()\n");
#endif
	matching_keys = &(seq2match_buf->matching_keys);
	n1 = matching_keys->nelt;
	// The number of elements in matching_keys can increase or decrease
	// during the for loop!
	for (i = 0; i < n1; i++) {
		k1 = matching_keys->elts[i];
		dups = VECTOR_ELT(low2high, k1);
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
						S, max_mm,
						count_only, seq2match_buf);
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
				S, max_mm,
				count_only, seq2match_buf);
		if (nmatches == 0) {
			IntAE_delete_at(matching_keys, i--);
			n1--;
		}
	}
#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] LEAVING match_headtail()\n");
#endif
	return;
}

static void match_pdict(SEXP pptb,
		CachedXStringSet *cached_head, CachedXStringSet *cached_tail,
		const RoSeq *S,
		SEXP max_mismatch, SEXP fixed,
		int count_only, Seq2MatchBuf *seq2match_buf)
{
	int tb_width, max_mm, fixedP, fixedS;
	SEXP low2high;
	const char *type;

#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] ENTERING match_pdict()\n");
#endif
	tb_width = _get_PreprocessedTB_width(pptb);
	low2high = _get_PreprocessedTB_low2high(pptb);
	type = get_classname(pptb);
	max_mm = INTEGER(max_mismatch)[0];
	fixedP = LOGICAL(fixed)[0];
	fixedS = LOGICAL(fixed)[1];
	if (strcmp(type, "Twobit") == 0)
		_match_Twobit(pptb, S, fixedS, seq2match_buf);
	else if (strcmp(type, "ACtree") == 0)
		_match_ACtree(pptb, S, fixedS, seq2match_buf);
	else if (strcmp(type, "ACtree2") == 0)
		_match_ACtree2(pptb, S, fixedS, seq2match_buf);
	else
		error("%s: unsupported Trusted Band type in 'pdict'", type);
	_select_nmismatch_at_Pshift_fun(fixedP, fixedS);
	/* Call match_headtail() even if 'cached_head' and 'cached_tail' are
         * NULL in order to propagate the matches to the duplicates */
	match_headtail(low2high, tb_width,
		cached_head, cached_tail,
		S, max_mm, count_only, seq2match_buf);
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
 *   matches_as: "DEVNULL", "MATCHES_AS_WHICH", "MATCHES_AS_COUNTS"
 *               or "MATCHES_AS_ENDS";
 *   envir: NULL or environment to be populated with the matches.
 *
 ****************************************************************************/

SEXP XString_match_pdict(SEXP pptb, SEXP pdict_head, SEXP pdict_tail,
		SEXP subject,
		SEXP max_mismatch, SEXP fixed,
		SEXP matches_as, SEXP envir)
{
	int tb_length, count_only;
	CachedXStringSet *cached_head, *cached_tail;
	RoSeq S;
	Seq2MatchBuf seq2match_buf;

#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] ENTERING XString_match_pdict()\n");
#endif
	tb_length = _get_PreprocessedTB_length(pptb);
	cached_head = get_CachedXStringSet_ptr(pdict_head);
	cached_tail = get_CachedXStringSet_ptr(pdict_tail);
	S = _get_XString_asRoSeq(subject);
	count_only = is_count_only(matches_as);

	seq2match_buf = _new_Seq2MatchBuf(matches_as, tb_length);
	match_pdict(pptb, cached_head, cached_tail,
		&S, max_mismatch, fixed,
		count_only, &seq2match_buf);
#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] LEAVING XString_match_pdict()\n");
#endif
	return _Seq2MatchBuf_as_SEXP(&seq2match_buf, envir);
}

SEXP XStringViews_match_pdict(SEXP pptb, SEXP pdict_head, SEXP pdict_tail,
		SEXP subject, SEXP views_start, SEXP views_width,
		SEXP max_mismatch, SEXP fixed,
		SEXP matches_as, SEXP envir)
{
	int tb_length, count_only;
	CachedXStringSet *cached_head, *cached_tail;
	RoSeq S, S_view;
	int nviews, i, *view_start, *view_width, view_offset;
	Seq2MatchBuf seq2match_buf, global_seq2match_buf;

#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] ENTERING XStringViews_match_pdict()\n");
#endif
	tb_length = _get_PreprocessedTB_length(pptb);
	cached_head = get_CachedXStringSet_ptr(pdict_head);
	cached_tail = get_CachedXStringSet_ptr(pdict_tail);
	S = _get_XString_asRoSeq(subject);
	count_only = is_count_only(matches_as);

	seq2match_buf = _new_Seq2MatchBuf(matches_as, tb_length);
	global_seq2match_buf = _new_Seq2MatchBuf(matches_as, tb_length);
	nviews = LENGTH(views_start);
	for (i = 0, view_start = INTEGER(views_start), view_width = INTEGER(views_width);
	     i < nviews;
	     i++, view_start++, view_width++)
	{
		view_offset = *view_start - 1;
		if (view_offset < 0 || view_offset + *view_width > S.nelt)
			error("'subject' has \"out of limits\" views");
		S_view.elts = S.elts + view_offset;
		S_view.nelt = *view_width;
		match_pdict(pptb, cached_head, cached_tail,
			&S_view, max_mismatch, fixed,
			count_only, &seq2match_buf);
		_Seq2MatchBuf_append_and_flush(&global_seq2match_buf,
			&seq2match_buf, view_offset);
	}

#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] LEAVING XStringViews_match_pdict()\n");
#endif
	return _Seq2MatchBuf_as_SEXP(&global_seq2match_buf, envir);
}


/****************************************************************************
 * .Call entry point: XStringSet_vmatch_pdict
 ****************************************************************************/

static SEXP vwhich_pdict(SEXP pptb,
		CachedXStringSet *cached_head, CachedXStringSet * cached_tail,
		SEXP subject,
		SEXP max_mismatch, SEXP fixed,
		Seq2MatchBuf *seq2match_buf)
{
	int tb_length, S_length, j;
	CachedXStringSet S;
	SEXP ans;
	RoSeq S_elt;

	tb_length = _get_PreprocessedTB_length(pptb);
	S = _new_CachedXStringSet(subject);
	S_length = _get_XStringSet_length(subject);
	PROTECT(ans = NEW_LIST(S_length));
	for (j = 0; j < S_length; j++) {
		S_elt = _get_CachedXStringSet_elt_asRoSeq(&S, j);
		match_pdict(pptb, cached_head, cached_tail,
			&S_elt, max_mismatch, fixed, 1, seq2match_buf);
		SET_ELEMENT(ans, j, _Seq2MatchBuf_which_asINTEGER(seq2match_buf));
		_flush_Seq2MatchBuf(seq2match_buf);
	}
	UNPROTECT(1);
	return ans;
}

static SEXP vcount_pdict_notcollapsed(SEXP pptb,
		CachedXStringSet *cached_head, CachedXStringSet * cached_tail,
		SEXP subject,
		SEXP max_mismatch, SEXP fixed,
		Seq2MatchBuf *seq2match_buf)
{
	int tb_length, S_length, j, *current_col;
	CachedXStringSet S;
	SEXP ans;
	RoSeq S_elt;
	IntAE *counts_buf;

	tb_length = _get_PreprocessedTB_length(pptb);
	S = _new_CachedXStringSet(subject);
	S_length = _get_XStringSet_length(subject);
	PROTECT(ans = allocMatrix(INTSXP, tb_length, S_length));
	for (j = 0, current_col = INTEGER(ans);
	     j < S_length;
	     j++, current_col += tb_length)
	{
		S_elt = _get_CachedXStringSet_elt_asRoSeq(&S, j);
		match_pdict(pptb, cached_head, cached_tail,
			&S_elt, max_mismatch, fixed, 1, seq2match_buf);
		counts_buf = &(seq2match_buf->match_counts);
		/* counts_buf->nelt is tb_length */
		memcpy(current_col, counts_buf->elts, sizeof(int) * counts_buf->nelt);
		_flush_Seq2MatchBuf(seq2match_buf);
	}
	UNPROTECT(1);
	return ans;
}

static SEXP vcount_pdict_collapsed(SEXP pptb,
		CachedXStringSet *cached_head, CachedXStringSet * cached_tail,
		SEXP subject,
		SEXP max_mismatch, SEXP fixed,
		int collapse, SEXP weight,
		Seq2MatchBuf *seq2match_buf)
{
	int tb_length, S_length, ans_length, i, j;
	CachedXStringSet S;
	SEXP ans;
	RoSeq S_elt;
	IntAE *counts_buf;

	tb_length = _get_PreprocessedTB_length(pptb);
	S = _new_CachedXStringSet(subject);
	S_length = _get_XStringSet_length(subject);
	switch (collapse) {
	    case 1: ans_length = tb_length; break;
	    case 2: ans_length = S_length; break;
	    default: error("'collapse' must be FALSE, 1 or 2");
	}
	if (IS_INTEGER(weight)) {
		PROTECT(ans = NEW_INTEGER(ans_length));
		memset(INTEGER(ans), 0, ans_length * sizeof(int));
	} else {
		PROTECT(ans = NEW_NUMERIC(ans_length));
		for (i = 0; i < ans_length; i++)
			REAL(ans)[i] = 0.00;
	}
	for (j = 0; j < S_length; j++)
	{
		S_elt = _get_CachedXStringSet_elt_asRoSeq(&S, j);
		match_pdict(pptb, cached_head, cached_tail,
			&S_elt, max_mismatch, fixed, 1, seq2match_buf);
		counts_buf = &(seq2match_buf->match_counts);
		/* counts_buf->nelt is tb_length */
		for (i = 0; i < tb_length; i++)
			if (collapse == 1) {
				if (IS_INTEGER(weight))
					INTEGER(ans)[i] += counts_buf->elts[i] * INTEGER(weight)[j];
				else
					REAL(ans)[i] += counts_buf->elts[i] * REAL(weight)[j];
			} else {
				if (IS_INTEGER(weight))
					INTEGER(ans)[j] += counts_buf->elts[i] * INTEGER(weight)[i];
				else
					REAL(ans)[j] += counts_buf->elts[i] * REAL(weight)[i];
			}
		_flush_Seq2MatchBuf(seq2match_buf);
	}
	UNPROTECT(1);
	return ans;
}

SEXP XStringSet_vmatch_pdict(SEXP pptb, SEXP pdict_head, SEXP pdict_tail,
		SEXP subject,
		SEXP max_mismatch, SEXP fixed,
		SEXP collapse, SEXP weight,
		SEXP matches_as, SEXP envir)
{
	int tb_length, count_only, collapse0;
	CachedXStringSet *cached_head, *cached_tail;
	Seq2MatchBuf seq2match_buf;
	SEXP ans;

#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] ENTERING XStringSet_vmatch_pdict()\n");
#endif
	tb_length = _get_PreprocessedTB_length(pptb);
	cached_head = get_CachedXStringSet_ptr(pdict_head);
	cached_tail = get_CachedXStringSet_ptr(pdict_tail);
	count_only = is_count_only(matches_as);

	seq2match_buf = _new_Seq2MatchBuf(matches_as, tb_length);
	switch (seq2match_buf.matches_as) {
	    case MATCHES_AS_NULL:
		error("XStringSet_vmatch_pdict() does not support "
		      "'matches_as=\"%s\"' yet, sorry", seq2match_buf.matches_as);
	    break;
	    case MATCHES_AS_WHICH:
		PROTECT(ans = vwhich_pdict(pptb,
				cached_head, cached_tail,
				subject,
				max_mismatch, fixed,
				&seq2match_buf));
	    break;
	    case MATCHES_AS_COUNTS:
		collapse0 = INTEGER(collapse)[0];
		if (collapse0 == 0)
			PROTECT(ans = vcount_pdict_notcollapsed(pptb,
					cached_head, cached_tail,
					subject,
					max_mismatch, fixed,
					&seq2match_buf));
		else
			PROTECT(ans = vcount_pdict_collapsed(pptb,
					cached_head, cached_tail,
					subject,
					max_mismatch, fixed, collapse0, weight,
					&seq2match_buf));
	    break;
	    case MATCHES_AS_ENDS:
		error("vmatchPDict() is not supported yet, sorry");
	    break;
	}
#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] LEAVING XStringSet_vmatch_pdict()\n");
#endif
	UNPROTECT(1);
	return ans;
}

