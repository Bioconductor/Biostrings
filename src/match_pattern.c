/****************************************************************************
 *                    EXACT AND INEXACT PATTERN MATCHING                    *
 *		             Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/****************************************************************************/
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

/* Return the number of matches */
static void match_naive_exact(RoSeq P, RoSeq S)
{
	int n1, n2;

	if (P.nelt <= 0)
		error("empty pattern");
	for (n1 = 0, n2 = P.nelt; n2 <= S.nelt; n1++, n2++, S.elts++) {
		if (memcmp(P.elts, S.elts, P.nelt) == 0)
			_Biostrings_report_match(n1, 0);
	}
	return;
}


/****************************************************************************
 * An implementation of the "naive" method for inexact matching.
 */

/* 
 * max_mm must be >= 0 (not safe otherwise)
 */
int _is_matching(RoSeq P, RoSeq S, int Pshift, int max_mm,
		int fixedP, int fixedS)
{
	int min_pm, mm, pm, i, chars_match;

	if (P.nelt <= max_mm)
		return 1;
	// 0 <= max_mm < P.nelt
	if (Pshift < 0) {
		max_mm += Pshift; // Pshift <= max_mm < P.nelt + Pshift < P.nelt
		if (max_mm < 0)
			return 0;
		// -P.nelt < Pshift < 0
		P.elts -= Pshift;
		P.nelt += Pshift; // 0 < P.nelt
	} else {
		S.elts += Pshift;
		S.nelt -= Pshift;
	}
	if (P.nelt > S.nelt) {
		max_mm -= P.nelt - S.nelt;
		if (max_mm < 0)
			return 0;
		P.nelt = S.nelt;
	}
	min_pm = P.nelt - max_mm;
	mm = pm = 0;
	// 0 = mm <= max_mm < P.nelt <= S.nelt
	// 0 = pm < min_pm <= P.nelt <= S.nelt
	// min_pm + max_mm = P.nelt
	for (i = 0; i < P.nelt; i++, P.elts++, S.elts++) {
		if (fixedP) {
			if (fixedS) {
				// *P.elts and *S.elts match iff they are equal
				chars_match = *P.elts == *S.elts;
			} else {
				// *P.elts and *S.elts match iff bits at 1
				// in *P.elts are are also at 1 in *S.elts
				chars_match = (*P.elts & ~(*S.elts)) == 0;
			}
		} else {
			if (fixedS) {
				// *P.elts and *S.elts match iff bits at 1
				// in *S.elts are also at 1 in *P.elts
				chars_match = (~(*P.elts) & *S.elts) == 0;
			} else {
				// *P.elts and *S.elts match iff they share
				// at least one bit at 1
				chars_match = *P.elts & *S.elts;
			}
		}
		if (chars_match) {
			if (++pm >= min_pm)
				return 1;
		} else {
			if (++mm > max_mm)
				return 0;
		}
	}
	error("Biostrings internal error in _is_matching(): "
	      "we should never be here");
	return -1;
}

static void match_naive_inexact(RoSeq P, RoSeq S,
		int max_mm, int fixedP, int fixedS)
{
	int n1, // position of pattern left-most char relative to the subject
	    n2, // 1 + position of pattern right-most char relative to the subject
	    min_n1, max_n2;

	if (P.nelt <= 0)
		error("empty pattern");
	min_n1 = P.nelt <= max_mm ? 1 - P.nelt : -max_mm;
	max_n2 = S.nelt - min_n1;
	for (n1 = min_n1, n2 = min_n1 + P.nelt; n2 <= max_n2; n1++, n2++)
		if (_is_matching(P, S, n1, max_mm, fixedP, fixedS))
			_Biostrings_report_match(n1, 0);
	return;
}

static void match_pattern(RoSeq P, RoSeq S, const char *algo,
		int max_mm, int fixedP, int fixedS)
{
	if (P.nelt > max_mm + S.nelt)
		return;
	if (P.nelt <= max_mm || strcmp(algo, "naive-inexact") == 0)
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
 *   'start': starting positions of the pattern relative to the subject
 *   'max_mismatch': the number of mismatches (integer vector of length 1)
 *   'fixed': logical vector of length 2
 *   'count_only': single logical
 *
 * is_matching() returns a logical vector of the same length as 'start'.
 * XString_match_pattern() returns an integer vector containing the relative
 * pos of the matches.
 * All matches have the length of the pattern.
 */

SEXP is_matching(SEXP pattern, SEXP subject, SEXP start,
		SEXP max_mismatch, SEXP fixed)
{
	RoSeq P, S;
	int start_len, max_mm, fixedP, fixedS, i, *start_elt, *ans_elt;
	SEXP ans;

	P = _get_XString_asRoSeq(pattern);
	S = _get_XString_asRoSeq(subject);
	start_len = LENGTH(start);
	max_mm = INTEGER(max_mismatch)[0];
	fixedP = LOGICAL(fixed)[0];
	fixedS = LOGICAL(fixed)[1];

	PROTECT(ans = NEW_LOGICAL(start_len));
	for (i = 0, start_elt = INTEGER(start), ans_elt = LOGICAL(ans);
             i < start_len;
             i++, start_elt++, ans_elt++) {
		if (*start_elt == NA_INTEGER) {
			*ans_elt = NA_LOGICAL;
			continue;
		}
 		*ans_elt = _is_matching(P, S, *start_elt - 1, max_mm,
				fixedP, fixedS);
	}
	UNPROTECT(1);
	return ans;
}

SEXP XString_match_pattern(SEXP pattern, SEXP subject, SEXP algorithm,
		SEXP max_mismatch, SEXP fixed,
		SEXP count_only)
{
	RoSeq P, S;
	const char *algo;
	int max_mm, fixedP, fixedS, is_count_only;

	P = _get_XString_asRoSeq(pattern);
	S = _get_XString_asRoSeq(subject);
	algo = CHAR(STRING_ELT(algorithm, 0));
	max_mm = INTEGER(max_mismatch)[0];
	fixedP = LOGICAL(fixed)[0];
	fixedS = LOGICAL(fixed)[1];
	is_count_only = LOGICAL(count_only)[0];

	_Biostrings_reset_viewsbuf(is_count_only ? 1 : 2);
	match_pattern(P, S, algo, max_mm, fixedP, fixedS);

	if (is_count_only)
		return _Biostrings_viewsbuf_count_asINTEGER();
	return _Biostrings_viewsbuf_start_asINTEGER();
}

SEXP XStringSet_match_pattern(SEXP pattern, SEXP subject, SEXP algorithm,
		SEXP max_mismatch, SEXP fixed, SEXP count_only)
{
        RoSeq P;
        CachedXStringSet S;
	const char *algo;
	int max_mm, fixedP, fixedS, is_count_only;

	P = _get_XString_asRoSeq(pattern);
	S = _new_CachedXStringSet(subject);
	algo = CHAR(STRING_ELT(algorithm, 0));
	max_mm = INTEGER(max_mismatch)[0];
	fixedP = LOGICAL(fixed)[0];
	fixedS = LOGICAL(fixed)[1];
	is_count_only = LOGICAL(count_only)[0];

	SEXP n_matches, counts;
	int length = _get_XStringSet_length(subject), i;

	if (is_count_only)
	    PROTECT(counts = NEW_INTEGER(length));
	else
	    PROTECT(counts = NEW_LIST(length));

	for (i = 0; i < length; ++i) {
	    _Biostrings_reset_viewsbuf(is_count_only ? 1 : 2);
	    match_pattern(P, _get_CachedXStringSet_elt_asRoSeq(&S, i), 
			  algo, max_mm, fixedP, fixedS);
	    if (is_count_only) {
		n_matches = _Biostrings_viewsbuf_count_asINTEGER();
		INTEGER(counts)[i] = INTEGER(n_matches)[0];
	    } else {
		PROTECT(n_matches = _Biostrings_viewsbuf_start_asINTEGER());
		SET_VECTOR_ELT(counts, i, n_matches);
		UNPROTECT(1);
	    }
	}

	UNPROTECT(1);
	return counts;
}

SEXP XStringViews_match_pattern(SEXP pattern,
		SEXP subject, SEXP views_start, SEXP views_width,
		SEXP algorithm, SEXP max_mismatch, SEXP fixed,
		SEXP count_only)
{
	RoSeq P, S, V;
	const char *algo;
	int max_mm, fixedP, fixedS, is_count_only,
	    nviews, i, *view_start, *view_width, view_offset;

	P = _get_XString_asRoSeq(pattern);
	S = _get_XString_asRoSeq(subject);
	algo = CHAR(STRING_ELT(algorithm, 0));
	max_mm = INTEGER(max_mismatch)[0];
	fixedP = LOGICAL(fixed)[0];
	fixedS = LOGICAL(fixed)[1];
	is_count_only = LOGICAL(count_only)[0];

	_Biostrings_reset_viewsbuf(is_count_only ? 1 : 2);
	nviews = LENGTH(views_start);
	for (i = 0, view_start = INTEGER(views_start), view_width = INTEGER(views_width);
	     i < nviews;
	     i++, view_start++, view_width++)
	{
		view_offset = *view_start - 1;
		if (view_offset < 0 || view_offset + *view_width > S.nelt)
			error("'subject' has out of limits views");
		V.elts = S.elts + view_offset;
		V.nelt = *view_width;
		_set_match_shift(view_offset);
		match_pattern(P, V, algo, max_mm, fixedP, fixedS);
	}

	if (is_count_only)
		return _Biostrings_viewsbuf_count_asINTEGER();
	return _Biostrings_viewsbuf_start_asINTEGER();
}

