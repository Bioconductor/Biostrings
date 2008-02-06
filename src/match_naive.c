/****************************************************************************
                 NAIVE METHODS FOR EXACT AND INEXACT MATCHING
		             Author: Herve Pages

 Here is how the "naive" (aka "memcmp", aka "blunt") method for finding exact
 matches is described in Dan Gusfield book "Algorithms on strings, trees, and
 sequences" (slightly modified):
     The naive method aligns the left end of P (the pattern) with the left
     end of S (the subject) and then compares the characters of P and S left
     to right until either two unequal characters are found or until P is
     exhausted, in which case an occurence of P is reported. In either case,
     P is then shifted one place to the right, and the comparisons are
     restarted from the left end of P. This process repeats until the right
     end of P shifts past the right end of S.
 
 Why implement this inefficient "naive" method?
   - For QC: we can validate other more sophisticated matching algo by
     comparing their results to those obtains with the "naive" method.
   - To use as a point of reference when comparing performance.

 ****************************************************************************/
#include "Biostrings.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/****************************************************************************/
static int debug = 0;

SEXP match_naive_debug()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'match_naive.c'\n", debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'match_naive.c'\n");
#endif
	return R_NilValue;
}


/****************************************************************************
 * A memcmp-based implementation of the "naive" method for exact matching
 * ======================================================================
 */

/* Return the number of matches */
static void naive_exact_search(const char *P, int nP, const char *S, int nS)
{
	int n1, n2;

	if (nP <= 0)
		error("empty pattern");
	for (n1 = 0, n2 = nP; n2 <= nS; n1++, n2++, S++) {
		if (memcmp(P, S, nP) == 0)
			_Biostrings_report_match(n1, 0);
	}
	return;
}


/****************************************************************************
 * An implementation of the "naive" method for inexact matching
 * ============================================================
 */

/* 
 * max_mm must be >= 0 (not safe otherwise)
 */
int _is_matching(const char *P, int nP, const char *S, int nS, int Pshift,
		int max_mm, int fixedP, int fixedS)
{
	int min_pm, mm, pm, i, chars_match;

	if (nP <= max_mm)
		return 1;
	// 0 <= max_mm < nP
	if (Pshift < 0) {
		max_mm += Pshift; // Pshift <= max_mm < nP + Pshift < nP
		if (max_mm < 0)
			return 0;
		// -nP < Pshift < 0
		P -= Pshift;
		nP += Pshift; // 0 < nP
	} else {
		S += Pshift;
		nS -= Pshift;
	}
	if (nP > nS) {
		max_mm -= nP - nS;
		if (max_mm < 0)
			return 0;
		nP = nS;
	}
	min_pm = nP - max_mm;
	mm = pm = 0;
	// 0 = mm <= max_mm < nP <= nS
	// 0 = pm < min_pm <= nP <= nS
	// min_pm + max_mm = nP
	for (i = 0; i < nP; i++, P++, S++) {
		if (fixedP) {
			if (fixedS) {
				// *P and *S match iff they are equal
				chars_match = *P == *S;
			} else {
				// *P and *S match iff bits at 1 in *P
				// are are also at 1 in *S
				chars_match = (*P & ~(*S)) == 0;
			}
		} else {
			if (fixedS) {
				// *P and *S match iff bits at 1 in *S
				// are also at 1 in *P
				chars_match = (~(*P) & *S) == 0;
			} else {
				// *P and *S match iff they share
				// at least one bit at 1
				chars_match = *P & *S;
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
	error("Biostrings internal error in _is_matching(): we should never be here");
	return -1;
}

static void naive_inexact_search(const char *P, int nP, const char *S, int nS,
		int max_mm, int fixedP, int fixedS)
{
	int n1, // position of pattern left-most char relative to the subject
	    n2, // 1 + position of pattern right-most char relative to the subject
	    max_n2;

	if (nP <= 0)
		error("empty pattern");
	max_n2 = nS + max_mm;
	for (n1 = -max_mm, n2 = nP - max_mm; n2 <= max_n2; n1++, n2++)
		if (_is_matching(P, nP, S, nS, n1, max_mm, fixedP, fixedS))
			_Biostrings_report_match(n1, 0);
	return;
}


/****************************************************************************
 * .Call entry points: "is_matching", "match_naive_exact" and "match_naive_inexact"
 *
 * Arguments:
 *   'pattern_BString': pattern
 *   'subject_BString': subject
 *   'start': starting positions of the pattern relative to the subject
 *   'max_mismatch': the number of mismatches (integer vector of length 1)
 *   'fixed': logical vector of length 2
 *   'count_only': single logical
 *
 * is_matching() return a logical vector of the same length as 'start'.
 * match_naive_exact() and match_naive_inexact() return an integer vector
 * containing the relative pos of the matches. All matches have the length of
 * the pattern.
 ****************************************************************************/

SEXP is_matching(SEXP pattern_BString, SEXP subject_BString, SEXP start,
		SEXP max_mismatch, SEXP fixed)
{
	const char *P, *S;
	int nP, nS, start_len, max_mm, fixedP, fixedS, i, *start_elt, *ans_elt;
	SEXP ans;

	P = getBString_charseq(pattern_BString, &nP);
	S = getBString_charseq(subject_BString, &nS);
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
		*ans_elt = _is_matching(P, nP, S, nS, *start_elt - 1, max_mm, fixedP, fixedS);
	}
	UNPROTECT(1);
	return ans;
}

SEXP match_naive_exact(SEXP pattern_BString, SEXP subject_BString,
		SEXP count_only)
{
	const char *P, *S;
	int nP, nS, is_count_only;
	SEXP ans;

	P = getBString_charseq(pattern_BString, &nP);
	S = getBString_charseq(subject_BString, &nS);
	is_count_only = LOGICAL(count_only)[0];

	_Biostrings_reset_viewsbuf(is_count_only ? 1 : 2);
	naive_exact_search(P, nP, S, nS);
	if (is_count_only)
		PROTECT(ans = _Biostrings_viewsbuf_count_asINTEGER());
	else
		PROTECT(ans = _Biostrings_viewsbuf_start_asINTEGER());
	UNPROTECT(1);
	return ans;
}

SEXP match_naive_inexact(SEXP pattern_BString, SEXP subject_BString,
		SEXP max_mismatch, SEXP fixed,
		SEXP count_only)
{
	const char *P, *S;
	int nP, nS, max_mm, fixedP, fixedS, is_count_only;
	SEXP ans;

	P = getBString_charseq(pattern_BString, &nP);
	S = getBString_charseq(subject_BString, &nS);
	max_mm = INTEGER(max_mismatch)[0];
	fixedP = LOGICAL(fixed)[0];
	fixedS = LOGICAL(fixed)[1];
	is_count_only = LOGICAL(count_only)[0];

	_Biostrings_reset_viewsbuf(is_count_only ? 1 : 2);
	naive_inexact_search(P, nP, S, nS, max_mm, fixedP, fixedS);
	if (is_count_only)
		PROTECT(ans = _Biostrings_viewsbuf_count_asINTEGER());
	else
		PROTECT(ans = _Biostrings_viewsbuf_start_asINTEGER());
	UNPROTECT(1);
	return ans;
}

