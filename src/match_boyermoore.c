/****************************************************************************
                      A BOYER-MOORE-LIKE MATCHING ALGO
		            Author: Herve Pages
 ****************************************************************************/
#include "Biostrings.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/****************************************************************************/
static int debug = 0;

SEXP match_boyermoore_debug()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'match_boyermoore.c'\n", debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'match_boyermoore.c'\n");
#endif
	return R_NilValue;
}


/****************************************************************************
 * A lookup table mapping any possible char to the position of its last
 * occurrence in the pattern P.
 * Positions are counted from the right (most right character in P being
 * at position 0). Chars not in P are mapped to position -1.
 */
static int rightPos_table[256];

static void prepare_rightPos_table(char *P, int nP)
{
	int u, j;
	unsigned char uc;
	char c;

	for (u = 0; u < 256; u++) {
		uc = (unsigned char) u;
		rightPos_table[uc] = -1;
		c = (char) uc;
		for (j = nP-1; j >= 0; j--) {
			if (P[j] == c) {
				rightPos_table[uc] = nP - 1 - j;
				break;
			}
		}
		/*
		printf("rightPos_table[%d] = %d", u, rightPos_table[uc]);
		if (rightPos_table[uc] != -1)
			printf(" (%c)", c);
		printf("\n");
		*/
	}
}


/****************************************************************************
 * The Good Right Shifts
 * =====================
 *
 * Definition
 * ----------
 *
 * For any j1 and j2 such that 0 <= j1 < j2 <= nP (nP being the length of the
 * pattern P), we note P(j1, j2) the subpattern of P starting at j1 and
 * ending at j2-1. Then GRS(j1, j2), the "Good Right Shift" for P(j1, j2), is
 * defined by:
 *   Imagine that, for a given alignement of P with the subject S, all
 *   letters in subpattern P(j1, j2) are matching S (we say that the current
 *   matching window is (j1, j2)). Then GRS(j1, j2) is the smallest (non-zero)
 *   amount of letters that we can shift P to the right without introducing an
 *   immediate mismatch.
 *
 * How to compute the GRS function
 * --------------------------------
 * 
 * The GRS function depends only on the pattern and therefore can be
 * preprocessed.
 * To compute GRS(j1, j2) we need to solve a matching problem again but within
 * P itself i.e. we must find the smallest (non-zero) amount of letters that we
 * need to shift P(j1, j2) to the _left_ until it matches P again.
 * IMPORTANT: This match must be a "full" match (when we are still within the
 * limits of P) and a "partial" match (when we are out of limits i.e. we've
 * moved to far to the left). For a "partial" match, all letters that are still
 * within P limits must match. When we have moved so far that all letters are
 * out of limits, then we have a "partial" match too.
 * This ensure that GRS(j1, j2) is always <= j2.
 * 
 * Example: P = acbaba
 *              012345
 * 
 *   GRS(0,1) = 1
 *   GRS(1,2) = 2
 *   GRS(2,3) = 3
 *   GRS(3.4) = GRS(2,4) = GRS(1,4) = GRS(0,4) = 3
 *   ...
 *   GRS(4,5) = GRS(5,6) = GRS(4,6) = 2
 *   GRS(0,6) = 5 (GRS max)
 * 
 * A nice property of the GRS function is that:
 *     if j1 <= j1' < j2' <= j2, then GRS(j1', j2') <= GRS(j1, j2)
 * Therefore the biggest value for GRS is reached for j1 = 0 and j2 = nP.
 * To summarize:
 *     1 <= GRS(j1, j2) <= min(j2, GRS(0, nP))
 * Another property of the GRS function is that if P1 and P2 are 2 patterns
 * such that P1 is a prefix of P2 then GRS1 and GRS2 (their respective GRS
 * functions) are equal on the set of points (j1, j2) that are valid for GRS1.
 * This last property is used by the prepare_GRS_table() C function below.
 * 
 * No need to preprocess the GRS function
 * --------------------------------------
 *
 * However this preprocessing would be too expensive since it requires the
 * evaluation of nP * (nP + 1) / 2 values, and each evaluation itself is a
 * string matching sub-problem with a cost of its own. So the trick is to
 * delay those evaluations until they are needed, and the fact is that, in
 * practise, very few of them are actually needed compared to the total
 * number of possible GRS(j1, j2) values.
 */

/* GRS_P contains the pattern associated with the current GRS_table. */
static char *GRS_P = NULL;
static int GRS_nP = 0;
/* GRS_table is a 2-dim array with nrow = ncol >= GRS_nP */
static int *GRS_table = NULL;
static int GRS_table_ncol = 0; /* GRS_table_ncol >= GRS_nP */

/* The layout of GRS_table is (only the values marked with an "x" are will
 * be potentially used):
 *
 *           1 2 3 4 5 6 j2
 *         0 x x x x - -
 *         1 - x x x - - 
 *         2 - - x x - -    GRS_nP = 4 <= GRS_table_ncol = 6
 *         3 - - - x - -
 *         4 - - - - - -
 *         5 - - - - - -
 *        j1
 *
 * The "x" region is defined by 0 <= j1 < j2 <= GRS_nP
 */

#define GRS(j1, j2)	(GRS_table[GRS_table_ncol * (j1) + (j2) - 1])

/*
 * The prepare_GRS_table() must ensure that the size of the GRS_P buffer
 * is always GRS_table_ncol */
static void prepare_GRS_table(char *P, int nP)
{
	int j1, j2, min_nP_GRS_nP;

	if (nP == 0) /* should never happen but safer anyway... */
		return;
	if (nP > 10000)
		error("pattern is too long");
	if (nP > GRS_table_ncol) {
		/* We need more memory */
		if (GRS_P != NULL)
			free(GRS_P);
		GRS_P = (char *) malloc(nP * sizeof(char));
		if (GRS_P == NULL)
			error("can't allocate memory for GRS_P");
		if (GRS_table != NULL)
			free(GRS_table);
		GRS_table = (int *) malloc(nP * nP * sizeof(int));
		if (GRS_table == NULL)
			error("can't allocate memory for GRS_table");
		GRS_table_ncol = nP;
		j2 = 1;
	} else {
		/* We have enough memory */
		if (nP < GRS_nP) min_nP_GRS_nP = nP; else min_nP_GRS_nP = GRS_nP;
		for (j2 = 0; j2 < min_nP_GRS_nP; j2++)
			if (P[j2] != GRS_P[j2])
				break;
		j2++;
	}
	memcpy(GRS_P, P, nP * sizeof(char));
	GRS_nP = nP;
	for ( ; j2 <= nP; j2++) {
		for (j1 = 0; j1 < j2; j1++) {
			GRS(j1, j2) = 0;
		}
	}
}

static int get_GRS(int j1, int j2)
{
	int grs, k1, k2, length;

	grs = GRS(j1, j2);
	if (grs != 0)
		return grs;
	for (grs = 1; grs < j2; grs++) {
		if (grs < j1) k1 = j1 - grs; else k1 = 0;
		k2 = j2 - grs;
		length = k2 - k1;
		if (memcmp(GRS_P + k1, GRS_P + j2 - length, length) == 0)
			break;
	}
	/* grs is j2 when the "for" loop is not interrupted by "break" */
	return GRS(j1, j2) = grs;
}


/****************************************************************************
 * The GRS-based exact matching algo
 * =================================
 * 
 * My poor understanding of the original Boyer-Moore algo as explained in Dan
 * Gusfield book "Algorithms on strings, trees, and sequences", is that it
 * will be inefficent in this situation:
 *   S = "Xbabababababab"
 *   P = "abababababab"
 * because after it has applied the strong good suffix rule and shifted P 2
 * letters to the right, it will start comparing P and S suffixes again which
 * is a waste of time.
 * The GRS_search algo below never compare twice the letters that are in the
 * current matching window (defined by (j1, j2) in P and (i1, i2) in S).
 */

/* Returns the number of matches */
static int GRS_search(char *P, int nP, char *S, int nS,
		int is_count_only, SEXP *p_matchpos, PROTECT_INDEX matchpos_pi)
{
	int count = 0, *index, n, i1, i2, j1, j2, shift0, grs, i, j;

	if (!is_count_only) {
		index = INTEGER(*p_matchpos);
	}
	prepare_rightPos_table(P, nP);
	prepare_GRS_table(P, nP);
	n = nP;
	i1 = i2 = n;
	j1 = j2 = nP;
	while (n <= nS) {
		if (j1 == j2) {
			/* No current partial match yet, we need to find one */
			shift0 = rightPos_table[(unsigned char) S[n-1]];
			if (shift0 == -1) {
				/* This is the best thing that can happen: the
				 * most-right letter of the pattern is aligned 
				 * with a letter in the subject that is not in
				 * the pattern.
				 */
				n += nP;
				i1 = i2 = n;
				j1 = j2 = nP;
				continue;
			}
			n += shift0; /* shift0 is always >= 0 and < nP */
			if (n > nS)
				break;
			i2 = n - shift0;
			j2 = nP - shift0;
			i1 = i2 - 1;
			j1 = j2 - 1;
			/* Now we have a partial match (of length 1) */
		}
		/* Let's try to extend the current partial match... */
		if (j1 > 0) {
			/* ... to the left */
			for (i = i1-1, j = j1-1; j >= 0; i--, j--)
				if (S[i] != P[j])
					break;
			i1 = i + 1;
			j1 = j + 1;
		}
		if (j2 < nP) {
			/* ... to the right */
			for ( ; j2 < nP; i2++, j2++)
				if (S[i2] != P[j2])
					break;
		}
		if (j1 == 0 && j2 == nP) {
			/* we have a full match! */
			if (!is_count_only) {
				if (count >= LENGTH(*p_matchpos)) {
					*p_matchpos = Biostrings_expandMatchIndex(
							*p_matchpos, i2, nS - i2);
					REPROTECT(*p_matchpos, matchpos_pi);
					index = INTEGER(*p_matchpos);
				}
				index[count] = i1;
			}
			count++;
		}
		grs = get_GRS(j1, j2); /* always >= 1 and <= min(j2, GRS(0, nP)) */
		n += grs;
		if (n > nS)
			break;
		if (grs <= j1) {
			j1 -= grs;
		} else {
			i1 += grs - j1;
			j1 = 0;
		}
	       	j2 -= grs;
	}
	return count;
}

SEXP match_boyermoore(SEXP p_xp, SEXP p_offset, SEXP p_length,
		SEXP s_xp, SEXP s_offset, SEXP s_length,
		SEXP count_only)
{
	int pat_offset, pat_length, subj_offset, subj_length,
	    is_count_only, count;
	char *pat, *subj;
	SEXP ans;
	int matchpos_length;
	SEXP matchpos = R_NilValue;
	PROTECT_INDEX matchpos_pi;

	pat_offset = INTEGER(p_offset)[0];
	pat_length = INTEGER(p_length)[0];
	pat = CHAR(R_ExternalPtrTag(p_xp)) + pat_offset;
	subj_offset = INTEGER(s_offset)[0];
	subj_length = INTEGER(s_length)[0];
	subj = CHAR(R_ExternalPtrTag(s_xp)) + subj_offset;
	is_count_only = LOGICAL(count_only)[0];

	if (!is_count_only) {
		matchpos_length = Biostrings_estimateExpectedMatchCount(
					pat_length, subj_length, 4);
		PROTECT_WITH_INDEX(matchpos, &matchpos_pi);
		matchpos = allocVector(INTSXP, matchpos_length);
		REPROTECT(matchpos, matchpos_pi);
	}
	count = GRS_search(pat, pat_length, subj, subj_length,
				is_count_only, &matchpos, matchpos_pi);
	if (is_count_only) {
		PROTECT(ans = allocVector(INTSXP, 1));
		INTEGER(ans)[0] = count;
		UNPROTECT(1);
	} else {
		PROTECT(ans = allocVector(INTSXP, count));
		memcpy(INTEGER(ans), INTEGER(matchpos), sizeof(int) * count);
		UNPROTECT(2);
	}
	return ans;
}

