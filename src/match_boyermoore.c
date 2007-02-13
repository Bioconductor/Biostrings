/****************************************************************************
                      A BOYER-MOORE-LIKE MATCHING ALGO
		            Author: Herve Pages
 ****************************************************************************/
#include "Biostrings.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* I've turned off the GRS feature. Doesn't seem to give very good
 * results :-( To turned it on, set MAXNP_WITH_GRS to a positive integer
 * (typically 128, 256 or 512, doesn't have to be a power of 2).
 * The original idea was that the "GRS feature" could boost the boyermoore()
 * function when the pattern is small.
 */
#define MAXNP_WITH_GRS	0


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
 * P0buffer holds a copy of the current pattern P0.
 * We must always have nP0 <= P0buffer_length
 */
static char *P0buffer = NULL;
static int P0buffer_length = 0, nP0 = 0;

/* Status meaning:
 *     -1: init_P0buffer() has changed the value of P0buffer_length.
 *   >= 0: init_P0buffer() didn't change the value of P0buffer_length.
 *         The non-negative integer is the length of the longest common
 *         suffix between old and new P0 (<= min(nP,nP0)).
 */
static int P0buffer_init_status;

static void init_P0buffer(char *P, int nP)
{
	int min_nP_nP0, j;

	if (nP == 0) {/* should never happen but safer anyway... */
		P0buffer_init_status = 0;
		return;
	}
	if (nP > 10000)
		error("pattern is too long");
	if (nP > P0buffer_length) {
		/* We need more memory */
		if (P0buffer != NULL)
			free(P0buffer);
		P0buffer_length = 0;
		P0buffer = (char *) malloc(nP * sizeof(char));
		if (P0buffer == NULL)
			error("can't allocate memory for P0buffer");
		P0buffer_length = nP;
		P0buffer_init_status = -1;
	} else {
		/* We have enough memory */
		if (nP < nP0) min_nP_nP0 = nP; else min_nP_nP0 = nP0;
		for (j = 0; j < min_nP_nP0; j++)
			if (P[j] != P0buffer[j])
				break;
		P0buffer_init_status = j;
	}
	memcpy(P0buffer, P, nP * sizeof(char));
	nP0 = nP;
	return;
}


/****************************************************************************
 * The Strong Good Suffix shifts
 * =============================
 */

/* Contains the Strong Good Suffix shifts for current pattern P0. */
static int *SGSshift_table = NULL;

/* SGSshift_table is a (256, P0buffer_length) matrix.
 * The layout of SGSshift_table is (only the values marked with an "x" will
 * be potentially used):
 *
 *           0 1 2 3 4 5 j
 *         0 x x x x - -
 *         1 x x x x - - 
 *         2 x x x x - -    nP0 = 4 <= P0buffer_length = 6
 *         .............
 *       256 x x x x - -
 *         c
 *
 * The "x" region is defined by 0 <= j < nP0
 */

#define SGSshift(c, j)	(SGSshift_table[P0buffer_length * ((unsigned char) c) + j])

/* Dumb and slow.
 * TODO: Use something faster like Cole's preprocessing algo.
 */
static int get_SGSshift(char c, int j)
{
	int shift, k, k1, k2, length;

	shift = SGSshift(c, j);
	if (shift != 0)
		return shift;
	for (shift = 1; shift < nP0; shift++) {
		if (shift <= j) {
			k = j - shift;
			if (P0buffer[k] != c)
				continue;
			k1 = k + 1;
		} else {
			k1 = 0;
		}
		k2 = nP0 - shift;
		if (k1 == k2)
			break;
		length = k2 - k1;
		if (memcmp(P0buffer + k1, P0buffer + k1 + shift, length) == 0)
			break;
	}
	/* shift is nP0 when the "for" loop is not interrupted by "break" */
	/*Rprintf("SGSshift(c=%c, j=%d) = %d\n", c, j, shift);*/
	return SGSshift(c, j) = shift;
}

static void init_SGSshift_table()
{
	int u, j;
	char c;

	if (P0buffer_init_status == -1 && SGSshift_table != NULL) {
		free(SGSshift_table);
		SGSshift_table = NULL;
	}
	if (P0buffer_length != 0 && SGSshift_table == NULL) {
		SGSshift_table = (int *) malloc(256 * P0buffer_length *
						sizeof(int));
		if (SGSshift_table == NULL)
			error("can't allocate memory for SGSshift_table");
	}
	for (u = 0; u < 256; u++) {
		for (j = 0; j < nP0; j++) {
			c = (char) u;
			SGSshift(c, j) = 0;
		}
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
 * This last property is used by the init_GRS_table() C function below.
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

/* Contains the GRS values for current pattern P0. */
static int *GRS_table = NULL;

/* GRS_table is a 2-dim array with nrow = ncol = P0buffer_length.
 * The layout of GRS_table is (only the values marked with an "x" will
 * be potentially used):
 *
 *           1 2 3 4 5 6 j2
 *         0 x x x x - -
 *         1 - x x x - - 
 *         2 - - x x - -    nP0 = 4 <= P0buffer_length = 6
 *         3 - - - x - -
 *         4 - - - - - -
 *         5 - - - - - -
 *        j1
 *
 * The "x" region is defined by 0 <= j1 < j2 <= nP0
 */

#define GRS(j1, j2)	(GRS_table[P0buffer_length * (j1) + (j2) - 1])

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
		if (memcmp(P0buffer + k1, P0buffer + j2 - length, length) == 0)
			break;
	}
	/* grs is j2 when the "for" loop is not interrupted by "break" */
	return GRS(j1, j2) = grs;
}

static void init_GRS_table()
{
	int j1, j2 = 1;

	if (P0buffer_init_status == -1 && GRS_table != NULL) {
		free(GRS_table);
		GRS_table = NULL;
	}
	if (P0buffer_length != 0 && GRS_table == NULL) {
		GRS_table = (int *) malloc(P0buffer_length * P0buffer_length *
						sizeof(int));
		if (GRS_table == NULL)
			error("can't allocate memory for GRS_table");
	}
	if (P0buffer_init_status != -1)
		j2 = P0buffer_init_status + 1;
	for ( ; j2 <= nP0; j2++) {
		for (j1 = 0; j1 < j2; j1++) {
			GRS(j1, j2) = 0;
		}
	}
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
 * The original idea of the "GRS feature" introduced in the boyermoore()
 * function below was to never compare twice the letters that are in the
 * current matching window (defined by (j1, j2) in P and (i1, i2) in S).
 */

#define ADJUSTMW(i, j, shift) \
{ \
	if ((shift) <= (j)) \
		(j) -= (shift); \
	else { \
		(i) += (shift) - (j); \
		(j) = 0; \
	} \
}

/* Returns the number of matches */
static int boyermoore(char *P, int nP, char *S, int nS, int is_count_only)
{
	int count = 0, n, i1, i2, j1, j2, shift0, shift1, i, j;
	char Pmrc, c; /* Pmrc is P most right char */

	init_P0buffer(P, nP);
	init_SGSshift_table();
	if (nP <= MAXNP_WITH_GRS)
		init_GRS_table();
	Pmrc = P[nP-1];
	n = nP - 1;
	j2 = 0;
	while (n < nS) {
		if (j2 == 0) {
			/* No matching window yet, we need to find one */
			c = S[n];
			if (c != Pmrc) {
				shift0 = get_SGSshift(c, nP-1);
				n += shift0;
				continue;
			}
			i1 = n;
			i2 = i1 + 1;
			j2 = nP;
			j1 = j2 - 1;
			/* Now we have a matching window (1-letter suffix) */
		}
		/* We try to extend the current matching window... */
		if (j1 > 0) {
			/* ... to the left */
			for (i = i1-1, j = j1-1; j >= 0; i--, j--)
				if ((c = S[i]) != P[j])
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
		if (j2 == nP) { /* matching window is a suffix */
			if (j1 == 0) {
				/* we have a full match! */
				if (!is_count_only)
					Biostrings_reportMatch(i1);
				count++;
				shift0 = get_SGSshift(P[0], 0); /* = max(GRS) */
			} else {
				shift0 = get_SGSshift(c, j1 - 1);
			}
		} else {
			shift0 = get_GRS(j1, j2);
			c = S[n];
			if (c != Pmrc) {
				shift1 = get_SGSshift(c, nP-1);
				if (shift1 > shift0)
					shift0 = shift1;
			}
		}
		n += shift0;
		if (nP <= MAXNP_WITH_GRS) {
			ADJUSTMW(i1, j1, shift0)
			ADJUSTMW(i2, j2, shift0)
		} else {
			j2 = 0; /* forget the current matching window */
		}
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

	pat_offset = INTEGER(p_offset)[0];
	pat_length = INTEGER(p_length)[0];
	pat = CHAR(R_ExternalPtrTag(p_xp)) + pat_offset;
	subj_offset = INTEGER(s_offset)[0];
	subj_length = INTEGER(s_length)[0];
	subj = CHAR(R_ExternalPtrTag(s_xp)) + subj_offset;
	is_count_only = LOGICAL(count_only)[0];

	if (!is_count_only)
		Biostrings_resetMatchPosBuffer();
	count = boyermoore(pat, pat_length, subj, subj_length, is_count_only);
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

