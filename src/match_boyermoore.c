/****************************************************************************
                      A BOYER-MOORE-LIKE MATCHING ALGO
		            Author: Herve Pages
 ****************************************************************************/
#include "Biostrings.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define LJ(j1, j2)	(LJ_table[LJ_table_N * (j1) + (j2)])

/*
 * A look up table mapping any possible char to the position of its last
 * occurrence in the pattern P.
 * Positions are counted from the right (most right character in P being
 * at position 0). Chars not in P are mapped to position -1.
 */
static int rightPos_table[256];

static int *LJ_table; /* 2-dim array with nrow = ncol */
static int LJ_table_N; /* = nrow = ncol = nP+1 */


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


/****************************************************************************/

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


static int prepare_LJ_table(char *P, int nP)
{
	int j1, j2, k1, k2, LJ_jump, jump, length;

	if (nP > 5000)
		error("pattern is too long");
	LJ_table_N = nP + 1;
	LJ_table = (int *) R_alloc(LJ_table_N * LJ_table_N, sizeof(int));
	for (j1 = 0; j1 < nP; j1++) {
		for (j2 = j1+1; j2 <= nP; j2++) {
			LJ_jump = j2;
			for (jump = 1; jump < LJ_jump; jump++) {
				if (jump < j1) k1 = j1 - jump; else k1 = 0;
				k2 = j2 - jump;
				length = k2 - k1;
				if (memcmp(P + k1, P + j2 - length, length) == 0) {
					LJ_jump = jump;
					break;
				}
			}
			LJ(j1, j2) = LJ_jump;
			/* printf("LJ(%d,%d) = %d\n", j1, j2, LJ(j1, j2)); */
		}
	}
	return 0;
}


/* Returns the number of matches */
static int LJsearch(char *P, int nP, char *S, int nS,
		int is_count_only, SEXP *p_matchpos, PROTECT_INDEX matchpos_pi)
{
	int count = 0, *index, n, i1, i2, j1, j2, jump, i, j;

	
	if (!is_count_only) {
		index = INTEGER(*p_matchpos);
	}
	prepare_rightPos_table(P, nP);
	prepare_LJ_table(P, nP);
	n = nP;
	i1 = i2 = n;
	j1 = j2 = nP;
	while (n <= nS) {
		if (j1 == j2) {
			/* No current partial match yet, we need to find one */
			jump = rightPos_table[(unsigned char) S[n-1]];
			if (jump == -1) {
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
			n += jump; /* jump is always >= 0 and < nP */
			if (n > nS)
				break;
			i2 = n - jump;
			j2 = nP - jump;
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
			printf("MATCH AT %d,%d\n", i1, i2);
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
		jump = LJ(j1, j2); /* always >= 1 and <= min(j2, LJ(0, nP)) */
		n += jump;
		if (n > nS)
			break;
		if (jump <= j1) {
			j1 -= jump;
		} else {
			i1 += jump - j1;
			j1 = 0;
		}
	       	j2 -= jump;
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
	count = LJsearch(pat, pat_length, subj, subj_length,
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

/*
int main()
{
	char *S, *P;
	int nS, nP;
	
	S = "XabcaabYabcaabcaabZ";
	P = "abcaab";
	nS = strlen(S);
	nP = strlen(P);
	printf("S = %s\n", S);
	printf("P = %s\n", P);
	LJsearch(P, nP, S, nS);

	// A bad case for the real Boyer-Moore algo
	S = "Xbabababababab";
	P = "abababababab";
	nS = strlen(S);
	nP = strlen(P);
	printf("S = %s\n", S);
	printf("P = %s\n", P);
	LJsearch(P, nP, S, nS);
	return 0;
}
*/


