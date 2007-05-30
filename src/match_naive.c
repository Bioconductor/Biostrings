/****************************************************************************
                      A NAIVE METHOD FOR EXACT MATCHING
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
 
 Why do we need a "naive" algo?
   - For QC: we can validate other more sophisticated matching algo by
     comparing their results to those obtains with the "naive" algo.
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
 * A memcmp-based implementation of the "naive" method
 * ===================================================
 */

/* Returns the number of matches */
static int naive_search(const char *P, int nP, const char *S, int nS, int is_count_only)
{
	int count = 0, n1, n2;

	n1 = 0;
	n2 = n1 + nP;
	for (n1 = 0, n2 = nP; n2 <= nS; n1++, n2++) {
		if (memcmp(P, S + n1, nP) != 0)
			continue;
		if (!is_count_only)
			Biostrings_reportMatch(n1);
		count++;
	}
	return count;
}

SEXP match_naive(SEXP p_xp, SEXP p_offset, SEXP p_length,
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
		Biostrings_resetMatchPosBuffer();
	count = naive_search((char *) pat, pat_length, (char *) subj, subj_length, is_count_only);
	if (!is_count_only) {
		PROTECT(ans = allocVector(INTSXP, count));
		memcpy(INTEGER(ans), Biostrings_resetMatchPosBuffer(),
					sizeof(int) * count);
	} else {
		PROTECT(ans = allocVector(INTSXP, 1));
		INTEGER(ans)[0] = count;
	}
	UNPROTECT(1);
	return ans;
}

