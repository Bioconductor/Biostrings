/****************************************************************************
                  NAIVE METHODS FOR EXACT AND FUZZY MATCHING
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
static int naive_exact_search(const char *P, int nP, const char *S, int nS,
		int is_count_only)
{
	int count = 0, n1, n2;

	for (n1 = 0, n2 = nP; n2 <= nS; n1++, n2++) {
		if (memcmp(P, S + n1, nP) != 0)
			continue;
		if (!is_count_only)
			_Biostrings_report_match(n1, 0);
		count++;
	}
	return count;
}


/****************************************************************************
 * An implementation of the "naive" method for fuzzy matching
 * ==========================================================
 */

/* Return the number of matches */
static int naive_fuzzy_search(const char *P, int nP, const char *S, int nS,
		int mm_max, int fixedP, int fixedS,
		int is_count_only)
{
	int count = 0,
	    n1, /* position of pattern left-most char relative to the subject */
	    n2, /* 1 + position of pattern right-most char relative to the subject */
	    n2_max, 
	    mm, /* current number of mismatches */
	    i, j;

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
		if (!is_count_only)
			_Biostrings_report_match(n1, 0);
		count++;
	}
	return count;
}


/****************************************************************************
 * .Call entry points: "match_naive_exact" and "match_naive_fuzzy"
 *
 * Arguments:
 *   'p_xp': pattern@data@xp
 *   'p_offset': pattern@offset
 *   'p_length': pattern@length
 *   's_xp': subject@data@xp
 *   's_offset': subject@offset
 *   's_length': subject@length
 *   'mismatch': the number of mismatches (integer vector of length 1)
 *   'fixed': logical vector of length 2
 *   'count_only': single logical
 *
 * The 2 functions return an integer vector containing the relative pos of
 * the matches. All matches have the length of the pattern.
 ****************************************************************************/

SEXP match_naive_exact(SEXP p_xp, SEXP p_offset, SEXP p_length,
		SEXP s_xp, SEXP s_offset, SEXP s_length,
		SEXP count_only)
{
	int pat_offset, pat_length, subj_offset, subj_length,
	    is_count_only, count;
	const Rbyte *pat, *subj;
	SEXP ans;

	pat_offset = INTEGER(p_offset)[0];
	pat_length = INTEGER(p_length)[0];
	pat = RAW(R_ExternalPtrTag(p_xp)) + pat_offset;
	subj_offset = INTEGER(s_offset)[0];
	subj_length = INTEGER(s_length)[0];
	subj = RAW(R_ExternalPtrTag(s_xp)) + subj_offset;
	is_count_only = LOGICAL(count_only)[0];

	if (!is_count_only)
		_Biostrings_reset_views_buffer();
	count = naive_exact_search(
			(char *) pat, pat_length,
			(char *) subj, subj_length,
			is_count_only);

	if (!is_count_only) {
		PROTECT(ans = NEW_INTEGER(count));
		memcpy(INTEGER(ans), _Biostrings_get_views_start(),
					sizeof(int) * count);
	} else {
		PROTECT(ans = NEW_INTEGER(1));
		INTEGER(ans)[0] = count;
	}
	UNPROTECT(1);
	return ans;
}

SEXP match_naive_fuzzy(SEXP p_xp, SEXP p_offset, SEXP p_length,
		SEXP s_xp, SEXP s_offset, SEXP s_length,
		SEXP mismatch, SEXP fixed,
		SEXP count_only)
{
	int pat_offset, pat_length, subj_offset, subj_length,
	    mm_max, fixedP, fixedS, is_count_only, count;
	const Rbyte *pat, *subj;
	SEXP ans;

	pat_offset = INTEGER(p_offset)[0];
	pat_length = INTEGER(p_length)[0];
	pat = RAW(R_ExternalPtrTag(p_xp)) + pat_offset;
	subj_offset = INTEGER(s_offset)[0];
	subj_length = INTEGER(s_length)[0];
	subj = RAW(R_ExternalPtrTag(s_xp)) + subj_offset;
	mm_max = INTEGER(mismatch)[0];
	fixedP = LOGICAL(fixed)[0];
	fixedS = LOGICAL(fixed)[1];
	is_count_only = LOGICAL(count_only)[0];

	if (!is_count_only)
		_Biostrings_reset_views_buffer();
	count = naive_fuzzy_search(
			(char *) pat, pat_length,
			(char *) subj, subj_length,
			mm_max, fixedP, fixedS,
			is_count_only);

	if (!is_count_only) {
		PROTECT(ans = NEW_INTEGER(count));
		memcpy(INTEGER(ans), _Biostrings_get_views_start(),
					sizeof(int) * count);
	} else {
		PROTECT(ans = NEW_INTEGER(1));
		INTEGER(ans)[0] = count;
	}
	UNPROTECT(1);
	return ans;
}
