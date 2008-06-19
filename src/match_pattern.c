/****************************************************************************
 *                    EXACT AND INEXACT PATTERN MATCHING                    *
 *		             Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
 *   exhausted, in which case an occurence of P is reported. In either case,
 *   P is then shifted one place to the right, and the comparisons are
 *   restarted from the left end of P. This process repeats until the right
 *   end of P shifts past the right end of S.
 *
 * Why implement this inefficient "naive" method?
 * - For QC: we can validate other more sophisticated matching algo by
 *   comparing their results to those obtains with the "naive" method.
 * - To use as a reference when comparing performance.
 */

static void match_naive_exact(const RoSeq *P, const RoSeq *S)
{
	const char *p, *s;
	int plen, slen, start, n2;

	if (P->nelt <= 0)
		error("empty pattern");
	p = P->elts;
	plen = P->nelt;
	s = S->elts;
	slen = S->nelt;
	for (start = 1, n2 = plen; n2 <= slen; start++, n2++, s++) {
		if (memcmp(p, s, plen) == 0)
			_report_match(start, 0);
	}
	return;
}


/****************************************************************************
 * An implementation of the "naive" method for inexact matching.
 */

static void match_naive_inexact(const RoSeq *P, const RoSeq *S,
		int max_mm, int fixedP, int fixedS)
{
	int n1, // position of pattern left-most char relative to the subject
	    n2, // 1 + position of pattern right-most char relative to the subject
	    min_n1, max_n2;

	if (P->nelt <= 0)
		error("empty pattern");
	min_n1 = P->nelt <= max_mm ? 1 - P->nelt : -max_mm;
	max_n2 = S->nelt - min_n1;
	for (n1 = min_n1, n2 = min_n1 + P->nelt; n2 <= max_n2; n1++, n2++)
		if (_is_matching_at_Pshift(P, S, n1, max_mm, fixedP, fixedS))
			_report_match(n1 + 1, 0);
	return;
}

static void match_pattern(const RoSeq *P, const RoSeq *S,
		SEXP algorithm,
		SEXP max_mismatch, SEXP fixed)
{
	const char *algo;
	int max_mm, fixedP, fixedS;

	max_mm = INTEGER(max_mismatch)[0];
	if (P->nelt > max_mm + S->nelt)
		return;
	algo = CHAR(STRING_ELT(algorithm, 0));
	fixedP = LOGICAL(fixed)[0];
	fixedS = LOGICAL(fixed)[1];
	if (P->nelt <= max_mm || strcmp(algo, "naive-inexact") == 0)
		match_naive_inexact(P, S, max_mm, fixedP, fixedS);
	else if (strcmp(algo, "naive-exact") == 0)
		match_naive_exact(P, S);
	else if (strcmp(algo, "boyer-moore") == 0)
		_match_pattern_boyermoore(P, S);
	else if (strcmp(algo, "shift-or") == 0)
		_match_pattern_shiftor(P, S, max_mm, fixedP, fixedS);
	else
		error("\"%s\": unknown algorithm", algo);
	return;
}


/****************************************************************************
 * --- .Call ENTRY POINTS ---
 *
 * Arguments:
 *   'pattern': pattern
 *   'subject': subject
 *   'algorithm': algorithm
 *   'max_mismatch': the max number of mismatching letters
 *   'fixed': logical vector of length 2
 *   'count_only': single logical
 *
 * All matches have the length of the pattern.
 */

SEXP XString_match_pattern(SEXP pattern, SEXP subject,
		SEXP algorithm,
		SEXP max_mismatch, SEXP fixed, SEXP count_only)
{
	RoSeq P, S;
	int is_count_only;

	P = _get_XString_asRoSeq(pattern);
	S = _get_XString_asRoSeq(subject);
	is_count_only = LOGICAL(count_only)[0];

	_init_match_reporting(is_count_only ? 1 : 2);
	match_pattern(&P, &S, algorithm, max_mismatch, fixed);
	return _reported_matches_asSEXP();
}

SEXP XStringViews_match_pattern(SEXP pattern,
		SEXP subject, SEXP views_start, SEXP views_width,
		SEXP algorithm,
		SEXP max_mismatch, SEXP fixed, SEXP count_only)
{
	RoSeq P, S, S_view;
	int is_count_only, nviews, i, *view_start, *view_width, view_offset;

	P = _get_XString_asRoSeq(pattern);
	S = _get_XString_asRoSeq(subject);
	is_count_only = LOGICAL(count_only)[0];

	_init_match_reporting(is_count_only ? 1 : 2);
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
		_set_match_shift(view_offset);
		match_pattern(&P, &S_view, algorithm, max_mismatch, fixed);
	}
	return _reported_matches_asSEXP();
}

SEXP XStringSet_vmatch_pattern(SEXP pattern, SEXP subject,
		SEXP algorithm,
		SEXP max_mismatch, SEXP fixed, SEXP count_only)
{
	RoSeq P, S_elt;
	CachedXStringSet S;
	int is_count_only, S_length, i;
	SEXP ans, ans_elt;

	P = _get_XString_asRoSeq(pattern);
	S = _new_CachedXStringSet(subject);
	is_count_only = LOGICAL(count_only)[0];

	_init_match_reporting(is_count_only ? 1 : 2);
	S_length = _get_XStringSet_length(subject);
	if (is_count_only)
		PROTECT(ans = NEW_INTEGER(S_length));
	else
		PROTECT(ans = NEW_LIST(S_length));
	for (i = 0; i < S_length; i++) {
		S_elt = _get_CachedXStringSet_elt_asRoSeq(&S, i);
		match_pattern(&P, &S_elt, algorithm, max_mismatch, fixed);
		PROTECT(ans_elt = _reported_matches_asSEXP());
		if (is_count_only)
			INTEGER(ans)[i] = INTEGER(ans_elt)[0];
		else
			SET_ELEMENT(ans, i, ans_elt);
		UNPROTECT(1);
		_drop_current_matches();
	}
	UNPROTECT(1);
	return ans;
}

