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



/****************************************************************************
 * Helper functions for all the *match_pdict() .Call entry points.
 */

// FIXME: Pass 'const HeadTail *' instead of 'pdict_head'
// and 'pdict_tail'. Also change MatchPDictBuf struct so the
// head_widths and tail_widths members are IntAE buffers.
static MatchPDictBuf new_MatchPDictBuf_from_PDict3Parts(SEXP matches_as,
		SEXP pptb, SEXP pdict_head, SEXP pdict_tail)
{
	int tb_length, tb_width;
	const int *head_widths, *tail_widths;

	tb_length = _get_PreprocessedTB_length(pptb);
	tb_width = _get_PreprocessedTB_width(pptb);
	if (pdict_head == R_NilValue)
		head_widths = NULL;
	else
		head_widths = INTEGER(_get_XStringSet_width(pdict_head));
	if (pdict_tail == R_NilValue)
		tail_widths = NULL;
	else
		tail_widths = INTEGER(_get_XStringSet_width(pdict_tail));
	return _new_MatchPDictBuf(matches_as, tb_length, tb_width,
				head_widths, tail_widths);
}

static void match_pdict(SEXP pptb, HeadTail *headtail, const cachedCharSeq *S,
		SEXP max_mismatch, SEXP min_mismatch, SEXP fixed,
		MatchPDictBuf *matchpdict_buf)
{
	int max_nmis, min_nmis, fixedP, fixedS;
	SEXP low2high;
	const char *type;
	TBMatchBuf *tb_matches;

	max_nmis = INTEGER(max_mismatch)[0];
	min_nmis = INTEGER(min_mismatch)[0];
	fixedP = LOGICAL(fixed)[0];
	fixedS = LOGICAL(fixed)[1];
	_select_nmismatch_at_Pshift_fun(fixedP, fixedS);
	type = get_classname(pptb);
/*
	if (strcmp(type, "ACtree2") == 0) {
		_match_pdictACtree2(pptb, headtail, S,
			max_nmis, min_nmis, fixedP, fixedS,
			matchpdict_buf);
		return;
	}
*/
#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] ENTERING match_pdict()\n");
#endif
	low2high = _get_PreprocessedTB_low2high(pptb);
	tb_matches = &(matchpdict_buf->tb_matches);

	if (strcmp(type, "Twobit") == 0)
		_match_Twobit(pptb, S, fixedS, tb_matches);
	else if (strcmp(type, "ACtree") == 0)
		_match_ACtree(pptb, S, fixedS, tb_matches);
	else if (strcmp(type, "ACtree2") == 0)
		_match_tbACtree2(pptb, S, fixedS, tb_matches);
	else
		error("%s: unsupported Trusted Band type in 'pdict'", type);
	/* Call _match_pdict_all_flanks() even if 'headtail' is empty
	 * (i.e. headtail->max_HTwidth == 0) because we need to propagate
	 * the matches to the duplicates anyway */
	_match_pdict_all_flanks(low2high, headtail,
		S, max_nmis, min_nmis, matchpdict_buf);
#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] LEAVING match_pdict()\n");
#endif
	return;
}



/****************************************************************************
 * Helper functions for the vcount_*_XStringSet all() functions.
 */

static SEXP init_vcount_collapsed_ans(int np, int ns,
		int collapse, SEXP weight)
{
	int ans_length, i;
	SEXP ans;

	switch (collapse) {
	    case 1: ans_length = np; break;
	    case 2: ans_length = ns; break;
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
	UNPROTECT(1);
	return ans;
}

static void update_vcount_collapsed_ans(SEXP ans,
		int match_count, int i, int j,
		int collapse, SEXP weight)
{
	int tmp;

	/* If 'collapse' is 1, then collapse the matrix of match counts
	   horizontally by summing all the (weighted) columns together. */
	if (collapse != 1) {
		/* Otherwise, collapse vertically. */
		tmp = i;
		i = j;
		j = tmp;
	}
	if (IS_INTEGER(weight))
		INTEGER(ans)[i] += match_count * INTEGER(weight)[j];
	else
		REAL(ans)[i] += match_count * REAL(weight)[j];
	return;
}



/****************************************************************************
 * .Call entry points: match_PDict3Parts_XString()
 *                     match_XStringSet_XString()
 *
 * Arguments:
 *   o match_PDict3Parts_XString() only:
 *     - pptb: a PreprocessedTB object;
 *     - pdict_head: head(pdict) (XStringSet or NULL);
 *     - pdict_tail: tail(pdict) (XStringSet or NULL);
 *   o match_XStringSet_XString() only:
 *     - pattern: non-preprocessed pattern dict (XStringSet);
 *   o common arguments:
 *     - subject: reference sequence (XString);
 *     - max_mismatch: max.mismatch (max nb of mismatches *outside* the TB);
 *     - min_mismatch: min.mismatch (min nb of mismatches *outside* the TB);
 *     - fixed: logical vector of length 2;
 *     - matches_as: "MATCHES_AS_NULL", "MATCHES_AS_WHICH",
 *         "MATCHES_AS_COUNTS" or "MATCHES_AS_ENDS";
 *     - envir: NULL or environment to be populated with the matches.
 */

/* --- .Call ENTRY POINT --- */
SEXP match_PDict3Parts_XString(SEXP pptb, SEXP pdict_head, SEXP pdict_tail,
		SEXP subject,
		SEXP max_mismatch, SEXP min_mismatch, SEXP fixed,
		SEXP matches_as, SEXP envir)
{
	HeadTail headtail;
	cachedCharSeq S;
	MatchPDictBuf matchpdict_buf;

	headtail = _new_HeadTail(pdict_head, pdict_tail, pptb,
				max_mismatch, fixed, 1);
	S = cache_XRaw(subject);
	matchpdict_buf = new_MatchPDictBuf_from_PDict3Parts(matches_as,
				pptb, pdict_head, pdict_tail);
	match_pdict(pptb, &headtail,
		&S, max_mismatch, min_mismatch, fixed,
		&matchpdict_buf);
	return _MatchBuf_as_SEXP(matchpdict_buf.ms_code,
			&(matchpdict_buf.matches), envir);
}

/* --- .Call ENTRY POINT --- */
SEXP match_XStringSet_XString(SEXP pattern,
		SEXP subject,
		SEXP max_mismatch, SEXP min_mismatch, SEXP fixed,
		SEXP matches_as, SEXP envir)
{
/*
	const char *ms_mode;
	int ms_code;

	cachedXStringSet P;
	int P_length, i;
	cachedCharSeq P_elt, S;
	MatchBuf ans_buf;

	ms_mode = CHAR(STRING_ELT(matches_as, 0));
	ms_code = _get_match_storing_code(ms_mode);
	switch (ms_code) {
	    case MATCHES_AS_WHICH: ms_mode = "MATCHES_AS_COUNTS"; break;
	    case MATCHES_AS_ENDS: ms_mode = "MATCHES_AS_RANGES"; break;
	}
	P = _cache_XStringSet(pattern);
	P_length = _get_cachedXStringSet_length(&P);
	S = cache_XRaw(subject);
	_init_match_reporting(ms_mode);
	for (i = 0; i < P_length; i++) {
		P_elt = _get_cachedXStringSet_elt(&P, i);
		_match_pattern(&P_elt, &S, NULL,
			max_mismatch, min_mismatch, NULL, fixed);
	}
*/
	error("vmatchPDict() is not supported yet, sorry");
	return R_NilValue;
}


/****************************************************************************
 * .Call entry points: match_PDict3Parts_XStringViews()
 *                     match_XStringSet_XStringViews()
 *
 * Arguments:
 *   o match_PDict3Parts_XStringViews() only:
 *     - pptb: a PreprocessedTB object;
 *     - pdict_head: head(pdict) (XStringSet or NULL);
 *     - pdict_tail: tail(pdict) (XStringSet or NULL);
 *   o match_XStringSet_XStringViews() only:
 *     - pattern: non-preprocessed pattern dict (XStringSet);
 *   o common arguments:
 *     - subject: reference sequence (XString);
 *     - views_start, views_width: views defined on the subject;
 *     - max_mismatch: max.mismatch (max nb of mismatches *outside* the TB);
 *     - min_mismatch: min.mismatch (min nb of mismatches *outside* the TB);
 *     - fixed: logical vector of length 2;
 *     - matches_as: "MATCHES_AS_NULL", "MATCHES_AS_WHICH",
 *         "MATCHES_AS_COUNTS" or "MATCHES_AS_ENDS";
 *     - envir: NULL or environment to be populated with the matches.
 */

/* --- .Call ENTRY POINT --- */
SEXP match_PDict3Parts_XStringViews(SEXP pptb, SEXP pdict_head, SEXP pdict_tail,
		SEXP subject, SEXP views_start, SEXP views_width,
		SEXP max_mismatch, SEXP min_mismatch, SEXP fixed,
		SEXP matches_as, SEXP envir)
{
	HeadTail headtail;
	int tb_length;
	cachedCharSeq S, S_view;
	int nviews, i, *view_start, *view_width, view_offset;
	MatchPDictBuf matchpdict_buf;
	MatchBuf global_matchpdict_buf;

	tb_length = _get_PreprocessedTB_length(pptb);
	headtail = _new_HeadTail(pdict_head, pdict_tail, pptb,
				max_mismatch, fixed, 1);
	S = cache_XRaw(subject);
	matchpdict_buf = new_MatchPDictBuf_from_PDict3Parts(matches_as,
				pptb, pdict_head, pdict_tail);
	global_matchpdict_buf = _new_MatchBuf(matchpdict_buf.ms_code,
				tb_length);
	nviews = LENGTH(views_start);
	for (i = 0,
	     view_start = INTEGER(views_start),
	     view_width = INTEGER(views_width);
	     i < nviews;
	     i++, view_start++, view_width++)
	{
		view_offset = *view_start - 1;
		if (view_offset < 0 || view_offset + *view_width > S.length)
			error("'subject' has \"out of limits\" views");
		S_view.seq = S.seq + view_offset;
		S_view.length = *view_width;
		match_pdict(pptb, &headtail, &S_view,
			    max_mismatch, min_mismatch, fixed,
			    &matchpdict_buf);
		_MatchPDictBuf_append_and_flush(&global_matchpdict_buf,
			&matchpdict_buf, view_offset);
	}
	return _MatchBuf_as_SEXP(matchpdict_buf.ms_code,
				&global_matchpdict_buf, envir);
}

/* --- .Call ENTRY POINT --- */
SEXP match_XStringSet_XStringViews(SEXP pattern,
		SEXP subject, SEXP views_start, SEXP views_width,
		SEXP max_mismatch, SEXP min_mismatch, SEXP fixed,
		SEXP matches_as, SEXP envir)
{
	error("match_XStringSet_XStringViews(): IMPLEMENT ME!");
	return R_NilValue;
}


/****************************************************************************
 * .Call entry points: vmatch_PDict3Parts_XStringSet()
 *                     vmatch_XStringSet_XStringSet()
 *
 * Arguments:
 *   o vmatch_PDict3Parts_XStringSet() only:
 *     - pptb: a PreprocessedTB object;
 *     - pdict_head: head(pdict) (XStringSet or NULL);
 *     - pdict_tail: tail(pdict) (XStringSet or NULL);
 *   o vmatch_XStringSet_XStringSet() only:
 *     - pattern: non-preprocessed pattern dict (XStringSet);
 *   o common arguments:
 *     - subject: reference sequences (XStringSet);
 *     - max_mismatch: max.mismatch (max nb of mismatches *outside* the TB);
 *     - min_mismatch: min.mismatch (min nb of mismatches *outside* the TB);
 *     - fixed: logical vector of length 2;
 *     - collapse, weight: controls how the matrix of match counts must be
 *         collapsed before it's returned to the user (ignored if
 *         'matches_as' != MATCHES_AS_COUNTS);
 *     - matches_as: "MATCHES_AS_NULL", "MATCHES_AS_WHICH",
 *         "MATCHES_AS_COUNTS" or "MATCHES_AS_ENDS";
 *     - envir: NULL or environment to be populated with the matches.
 */

static SEXP vwhich_PDict3Parts_XStringSet(SEXP pptb, HeadTail *headtail,
		SEXP subject,
		SEXP max_mismatch, SEXP min_mismatch, SEXP fixed,
		MatchPDictBuf *matchpdict_buf)
{
	int S_length, j;
	cachedXStringSet S;
	SEXP ans, ans_elt;
	cachedCharSeq S_elt;

	S = _cache_XStringSet(subject);
	S_length = _get_cachedXStringSet_length(&S);
	PROTECT(ans = NEW_LIST(S_length));
	for (j = 0; j < S_length; j++) {
		S_elt = _get_cachedXStringSet_elt(&S, j);
		match_pdict(pptb, headtail, &S_elt,
			    max_mismatch, min_mismatch, fixed,
			    matchpdict_buf);
		PROTECT(ans_elt = _MatchBuf_which_asINTEGER(
					&(matchpdict_buf->matches)));
		SET_ELEMENT(ans, j, ans_elt);
		UNPROTECT(1);
		_MatchPDictBuf_flush(matchpdict_buf);
	}
	UNPROTECT(1);
	return ans;
}

static SEXP vwhich_XStringSet_XStringSet(SEXP pattern,
		SEXP subject,
		SEXP max_mismatch, SEXP min_mismatch, SEXP fixed)
{
	cachedXStringSet P, S;
	int P_length, S_length, i, j;
	cachedCharSeq P_elt, S_elt;
	IntAEAE ans_buf;

	P = _cache_XStringSet(pattern);
	P_length = _get_cachedXStringSet_length(&P);
	S = _cache_XStringSet(subject);
	S_length = _get_cachedXStringSet_length(&S);
	ans_buf = new_IntAEAE(S_length, S_length);
	for (j = 0; j < S_length; j++)
		ans_buf.elts[j].nelt = 0;
	_init_match_reporting("MATCHES_AS_COUNTS");
	for (i = 0; i < P_length; i++) {
		P_elt = _get_cachedXStringSet_elt(&P, i);
		for (j = 0; j < S_length; j++) {
			S_elt = _get_cachedXStringSet_elt(&S, j);
			_match_pattern(&P_elt, &S_elt, NULL,
				max_mismatch, min_mismatch, NULL, fixed);
			if (_get_match_count() != 0)
				IntAE_insert_at(ans_buf.elts + j,
					ans_buf.elts[j].nelt, i + 1);
			_drop_reported_matches();
		}
	}
	return IntAEAE_asLIST(&ans_buf, 1);
}

static SEXP vcount_PDict3Parts_XStringSet(SEXP pptb, HeadTail *headtail,
		SEXP subject,
		SEXP max_mismatch, SEXP min_mismatch, SEXP fixed,
		SEXP collapse, SEXP weight,
		MatchPDictBuf *matchpdict_buf)
{
	int tb_length, S_length, collapse0, i, j, match_count, *ans_col;
	cachedXStringSet S;
	SEXP ans;
	cachedCharSeq S_elt;
	const IntAE *count_buf;

	tb_length = _get_PreprocessedTB_length(pptb);
	S = _cache_XStringSet(subject);
	S_length = _get_cachedXStringSet_length(&S);
	collapse0 = INTEGER(collapse)[0];
	if (collapse0 == 0) {
		PROTECT(ans = allocMatrix(INTSXP, tb_length, S_length));
		ans_col = INTEGER(ans);
	} else {
		PROTECT(ans = init_vcount_collapsed_ans(tb_length, S_length,
					collapse0, weight));
	}
	for (j = 0; j < S_length; j++) {
		S_elt = _get_cachedXStringSet_elt(&S, j);
		match_pdict(pptb, headtail, &S_elt,
			max_mismatch, min_mismatch, fixed,
			matchpdict_buf);
		count_buf = &(matchpdict_buf->matches.match_counts);
		/* 'count_buf->nelt' is 'tb_length' */
		if (collapse0 == 0) {
			memcpy(ans_col, count_buf->elts,
				sizeof(int) * count_buf->nelt);
			ans_col += tb_length;
		} else {
			for (i = 0; i < tb_length; i++) {
				match_count = count_buf->elts[i];
				update_vcount_collapsed_ans(ans,
					match_count, i, j,
					collapse0, weight);
			}
		}
		_MatchPDictBuf_flush(matchpdict_buf);
	}
	UNPROTECT(1);
	return ans;
}

static SEXP vcount_XStringSet_XStringSet(SEXP pattern,
		SEXP subject,
		SEXP max_mismatch, SEXP min_mismatch, SEXP fixed,
		SEXP collapse, SEXP weight)
{
	cachedXStringSet P, S;
	int P_length, S_length, collapse0, i, j, match_count, *ans_elt;
	SEXP ans;
	cachedCharSeq P_elt, S_elt;

	P = _cache_XStringSet(pattern);
	P_length = _get_cachedXStringSet_length(&P);
	S = _cache_XStringSet(subject);
	S_length = _get_cachedXStringSet_length(&S);
	collapse0 = INTEGER(collapse)[0];
	if (collapse0 == 0)
		PROTECT(ans = allocMatrix(INTSXP, P_length, S_length));
	else
		PROTECT(ans = init_vcount_collapsed_ans(P_length, S_length,
					collapse0, weight));
	_init_match_reporting("MATCHES_AS_COUNTS");
	for (i = 0; i < P_length; i++) {
		P_elt = _get_cachedXStringSet_elt(&P, i);
		if (collapse0 == 0)
			ans_elt = INTEGER(ans) + i;
		for (j = 0; j < S_length; j++) {
			S_elt = _get_cachedXStringSet_elt(&S, j);
			_match_pattern(&P_elt, &S_elt, NULL,
				max_mismatch, min_mismatch, NULL, fixed);
			match_count = _get_match_count();
			if (collapse0 == 0) {
				*ans_elt = match_count;
				ans_elt += P_length;
			} else {
				update_vcount_collapsed_ans(ans,
					match_count, i, j,
					collapse0, weight);
			}
			_drop_reported_matches();
		}
	}
	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP vmatch_PDict3Parts_XStringSet(SEXP pptb, SEXP pdict_head, SEXP pdict_tail,
		SEXP subject,
		SEXP max_mismatch, SEXP min_mismatch, SEXP fixed,
		SEXP collapse, SEXP weight,
		SEXP matches_as, SEXP envir)
{
	HeadTail headtail;
	MatchPDictBuf matchpdict_buf;

	headtail = _new_HeadTail(pdict_head, pdict_tail, pptb,
				max_mismatch, fixed, 1);
	matchpdict_buf = new_MatchPDictBuf_from_PDict3Parts(matches_as,
				pptb, pdict_head, pdict_tail);
	switch (matchpdict_buf.ms_code) {
	    case MATCHES_AS_NULL:
		error("vmatch_PDict3Parts_XStringSet() does not support "
		      "'matches_as=\"%s\"' yet, sorry",
		      CHAR(STRING_ELT(matches_as, 0)));
	    case MATCHES_AS_WHICH:
		return vwhich_PDict3Parts_XStringSet(pptb, &headtail,
				subject,
				max_mismatch, min_mismatch, fixed,
				&matchpdict_buf);
	    case MATCHES_AS_COUNTS:
		return vcount_PDict3Parts_XStringSet(pptb, &headtail,
				subject,
				max_mismatch, min_mismatch, fixed,
				collapse, weight,
				&matchpdict_buf);
	}
	error("vmatchPDict() is not supported yet, sorry");
	return R_NilValue;
}

/* --- .Call ENTRY POINT --- */
SEXP vmatch_XStringSet_XStringSet(SEXP pattern,
		SEXP subject,
		SEXP max_mismatch, SEXP min_mismatch, SEXP fixed,
		SEXP collapse, SEXP weight,
		SEXP matches_as, SEXP envir)
{
	const char *ms_mode;
	int ms_code;

	ms_mode = CHAR(STRING_ELT(matches_as, 0));
	ms_code = _get_match_storing_code(ms_mode);
	switch (ms_code) {
	    case MATCHES_AS_NULL:
		error("vmatch_XStringSet_XStringSet() does not support "
		      "'matches_as=\"%s\"' yet, sorry", ms_mode);
	    case MATCHES_AS_WHICH:
		return vwhich_XStringSet_XStringSet(pattern,
				subject,
				max_mismatch, min_mismatch, fixed);
	    case MATCHES_AS_COUNTS:
		return vcount_XStringSet_XStringSet(pattern,
				subject,
				max_mismatch, min_mismatch, fixed,
				collapse, weight);
	}
	error("vmatchPDict() is not supported yet, sorry");
	return R_NilValue;
}

