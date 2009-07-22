/****************************************************************************
 *                    EXACT AND INEXACT PATTERN MATCHING                    *
 *		             Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"
#include "IRanges_interface.h"

static int debug = 0;

SEXP debug_match_pattern()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'match_pattern.c'\n",
                 debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'match_pattern.c'\n");
#endif
	return R_NilValue;
}


/****************************************************************************
 * A memcmp-based implementation of the "naive" method for exact matching.
 *
 * Here is how the "naive" (aka "memcmp", aka "blunt") method for finding
 * exact matches is described in Dan Gusfield book "Algorithms on strings,
 * trees, and sequences" (slightly modified):
 *   The naive method aligns the left end of P (the pattern) with the left
 *   end of S (the subject) and then compares the characters of P and S left
 *   to right until either two unequal characters are found or until P is
 *   exhausted, in which case an occurrence of P is reported. In either case,
 *   P is then shifted one place to the right, and the comparisons are
 *   restarted from the left end of P. This process repeats until the right
 *   end of P shifts past the right end of S.
 *
 * Why implement this inefficient "naive" method?
 * - For QC: we can validate other more sophisticated matching algo by
 *   comparing their results to those obtains with the "naive" method.
 * - To use as a reference when comparing performance.
 */

static void match_naive_exact(const cachedCharSeq *P, const cachedCharSeq *S)
{
	const char *p, *s;
	int plen, slen, start, n2;

	if (P->length <= 0)
		error("empty pattern");
	p = P->seq;
	plen = P->length;
	s = S->seq;
	slen = S->length;
	for (start = 1, n2 = plen; n2 <= slen; start++, n2++, s++) {
		if (memcmp(p, s, plen) == 0)
			_report_match(start, P->length);
	}
	return;
}


/****************************************************************************
 * An implementation of the "naive" method for inexact matching.
 */

static void match_naive_inexact(const cachedCharSeq *P, const cachedCharSeq *S,
		int max_mm, int fixedP, int fixedS)
{
	int Pshift, // position of pattern left-most char relative to the subject
	    n2, // 1 + position of pattern right-most char relative to the subject
	    min_Pshift, max_n2;

	if (P->length <= 0)
		error("empty pattern");
	_select_nmismatch_at_Pshift_fun(fixedP, fixedS);
	min_Pshift = P->length <= max_mm ? 1 - P->length : -max_mm;
	max_n2 = S->length - min_Pshift;
	for (Pshift = min_Pshift, n2 = min_Pshift + P->length; n2 <= max_n2; Pshift++, n2++)
		if (_selected_nmismatch_at_Pshift_fun(P, S, Pshift, max_mm) <= max_mm)
			_report_match(Pshift + 1, P->length);
	return;
}


/****************************************************************************
 * match_pattern()
 */

static void match_pattern(const cachedCharSeq *P, const cachedCharSeq *S, SEXP algorithm,
		SEXP max_mismatch, SEXP with_indels, SEXP fixed)
{
	const char *algo;
	int max_mm, fixedP, fixedS;

	max_mm = INTEGER(max_mismatch)[0];
	if (P->length > max_mm + S->length)
		return;
	algo = CHAR(STRING_ELT(algorithm, 0));
	fixedP = LOGICAL(fixed)[0];
	fixedS = LOGICAL(fixed)[1];
	if (P->length <= max_mm || strcmp(algo, "naive-inexact") == 0)
		match_naive_inexact(P, S, max_mm, fixedP, fixedS);
	else if (strcmp(algo, "naive-exact") == 0)
		match_naive_exact(P, S);
	else if (strcmp(algo, "boyer-moore") == 0)
		_match_pattern_boyermoore(P, S);
	else if (strcmp(algo, "shift-or") == 0)
		_match_pattern_shiftor(P, S, max_mm, fixedP, fixedS);
	else if (strcmp(algo, "indels") == 0)
		_match_pattern_indels(P, S, max_mm, fixedP, fixedS);
	else
		error("\"%s\": unknown algorithm", algo);
	return;
}


/****************************************************************************
 * --- .Call ENTRY POINTS ---
 *
 * Arguments:
 *   pattern, subject: XString objects;
 *   algorithm: single string;
 *   max_mismatch: (single integer) the max number of mismatching letters;
 *   with_indels: single logical;
 *   fixed: logical vector of length 2;
 *   count_only: single logical.
 *
 * If with_indels is FALSE: all matches have the length of the pattern.
 * Otherwise, matches are of variable length (>= length(pattern) - max_mismatch
 * and <= length(pattern) + max_mismatch).
 */

SEXP XString_match_pattern(SEXP pattern, SEXP subject, SEXP algorithm,
		SEXP max_mismatch, SEXP with_indels, SEXP fixed, SEXP count_only)
{
	cachedCharSeq P, S;
	int is_count_only;

	P = cache_XRaw(pattern);
	S = cache_XRaw(subject);
	is_count_only = LOGICAL(count_only)[0];

	_init_match_reporting(is_count_only ? mkString("COUNTONLY") : mkString("ASIRANGES"));
	match_pattern(&P, &S, algorithm, max_mismatch, with_indels, fixed);
	return _reported_matches_asSEXP();
}

/*
 * Arguments are the same as for XString_match_pattern() except for:
 *   subject: XString object;
 *   views_start, views_width: 2 integer vectors describing views on 'subject'.
 */
SEXP XStringViews_match_pattern(SEXP pattern,
		SEXP subject, SEXP views_start, SEXP views_width,
		SEXP algorithm,
		SEXP max_mismatch, SEXP with_indels, SEXP fixed, SEXP count_only)
{
	cachedCharSeq P, S, S_view;
	int is_count_only, nviews, i, *view_start, *view_width, view_offset;

	P = cache_XRaw(pattern);
	S = cache_XRaw(subject);
	is_count_only = LOGICAL(count_only)[0];

	_init_match_reporting(is_count_only ? mkString("COUNTONLY") : mkString("ASIRANGES"));
	nviews = LENGTH(views_start);
	for (i = 0, view_start = INTEGER(views_start), view_width = INTEGER(views_width);
	     i < nviews;
	     i++, view_start++, view_width++)
	{
		view_offset = *view_start - 1;
		if (view_offset < 0 || view_offset + *view_width > S.length)
			error("'subject' has \"out of limits\" views");
		S_view.seq = S.seq + view_offset;
		S_view.length = *view_width;
		_shift_match_on_reporting(view_offset);
		match_pattern(&P, &S_view, algorithm, max_mismatch, with_indels, fixed);
	}
	return _reported_matches_asSEXP();
}

/*
 * Arguments are the same as for XString_match_pattern() except for:
 *   subject: XStringSet object.
 */
SEXP XStringSet_vmatch_pattern(SEXP pattern, SEXP subject, SEXP algorithm,
		SEXP max_mismatch, SEXP with_indels, SEXP fixed, SEXP count_only)
{
	cachedCharSeq P, S_elt;
	cachedXStringSet S;
	int is_count_only, S_length, i;
	SEXP ans, ans_elt;

	P = cache_XRaw(pattern);
	S = _cache_XStringSet(subject);
	is_count_only = LOGICAL(count_only)[0];

	_init_match_reporting(is_count_only ? mkString("COUNTONLY") : mkString("ASIRANGES"));
	S_length = _get_XStringSet_length(subject);
	if (is_count_only)
		PROTECT(ans = NEW_INTEGER(S_length));
	else
		PROTECT(ans = NEW_LIST(S_length));
	for (i = 0; i < S_length; i++) {
		S_elt = _get_cachedXStringSet_elt(&S, i);
		match_pattern(&P, &S_elt, algorithm, max_mismatch, with_indels, fixed);
		PROTECT(ans_elt = _reported_matches_asSEXP());
		if (is_count_only)
			INTEGER(ans)[i] = INTEGER(ans_elt)[0];
		else
			SET_ELEMENT(ans, i, ans_elt);
		UNPROTECT(1);
		_drop_reported_matches();
	}
	UNPROTECT(1);
	return ans;
}

