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

/* Return the number of matches */
static void naive_inexact_search(const char *P, int nP, const char *S, int nS,
		int mm_max, int fixedP, int fixedS)
{
	int n1, /* position of pattern left-most char relative to the subject */
	    n2, /* 1 + position of pattern right-most char relative to the subject */
	    n2_max, 
	    mm, /* current number of mismatches */
	    i, j;

	if (nP <= 0)
		error("empty pattern");
	n2_max = nS + mm_max;
	for (n1 = -mm_max, n2 = nP - mm_max; n2 <= n2_max; n1++, n2++) {
		mm = 0;
		for (i = n1, j = 0; j < nP; i++, j++) {
			if (i < 0 || nS <= i) {
				mm++;
			} else {
				if (fixedS) {
					if (fixedP) {
						/* S[i] and P[j] match iff they
						   are equal */
						if (S[i] != P[j]) mm++;
					} else {
						/* S[i] and P[j] match iff bits at 1
						   in S[i] are a subset of bits at 1
						   in P[j] */
						if (S[i] & (~P[j])) mm++;
					}
				} else {
					if (fixedP) {
						/* S[i] and P[j] match iff bits at 1
						   in P[j] are a subset of bits at 1
						   in S[i] */
						if ((~S[i]) & P[j]) mm++;
					} else {
						/* S[i] and P[j] match iff they have
						   at least one bit at 1 in common */
						if ((S[i] & P[j]) == 0) mm++;
					}
				}
			}
			if (mm > mm_max)
				break;
		}
		/*
		Rprintf("nS=%d nP=%d mm_max=%d n1=%d n2=%d mm=%d i=%d j=%d\n",
			nS, nP, mm_max, n1, n2, mm, i, j);
		*/
		if (j < nP)
			continue;
		_Biostrings_report_match(n1, 0);
	}
	return;
}


/****************************************************************************
 * .Call entry points: "match_naive_exact" and "match_naive_inexact"
 *
 * Arguments:
 *   'pattern_BString': pattern
 *   'subject_BString': subject
 *   'mismatch': the number of mismatches (integer vector of length 1)
 *   'fixed': logical vector of length 2
 *   'count_only': single logical
 *
 * The 2 functions return an integer vector containing the relative pos of
 * the matches. All matches have the length of the pattern.
 ****************************************************************************/

SEXP match_naive_exact(SEXP pattern_BString, SEXP subject_BString,
		SEXP count_only)
{
	const char *P, *S;
	int nP, nS, is_count_only;
	SEXP ans;

	P = get_BString_seq(pattern_BString, &nP);
	S = get_BString_seq(subject_BString, &nS);
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
		SEXP mismatch, SEXP fixed,
		SEXP count_only)
{
	const char *P, *S;
	int nP, nS, mm_max, fixedP, fixedS, is_count_only;
	SEXP ans;

	P = get_BString_seq(pattern_BString, &nP);
	S = get_BString_seq(subject_BString, &nS);
	mm_max = INTEGER(mismatch)[0];
	fixedP = LOGICAL(fixed)[0];
	fixedS = LOGICAL(fixed)[1];
	is_count_only = LOGICAL(count_only)[0];

	_Biostrings_reset_viewsbuf(is_count_only ? 1 : 2);
	naive_inexact_search(P, nP, S, nS, mm_max, fixedP, fixedS);
	if (is_count_only)
		PROTECT(ans = _Biostrings_viewsbuf_count_asINTEGER());
	else
		PROTECT(ans = _Biostrings_viewsbuf_start_asINTEGER());
	UNPROTECT(1);
	return ans;
}
