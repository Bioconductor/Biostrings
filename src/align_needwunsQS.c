#include "Biostrings.h"

#include <stdio.h>


/* Returns the score of the alignment */
static int needwunsQS(char **pal1, char **pal2, int *pnal,
		const char *S1, int nS1, const char *S2, int nS2,
		const int *mat, int mat_nrow, const int *lkup, int lkup_length,
		int gap_cost, char gap_code)
{
	int score;
	static char *al1 = "aaC", *al2 = "zzX";

	*pal1 = al1;
	*pal2 = al2;
	*pnal = 3;
	score = -99;
	return score;
}

/*
 * 's1_xp', 's1_offset', 's1_length': left BString object
 * 's2_xp', 's2_offset', 's2_length': right BString object
 * 'mat': scoring matrix (integer square matrix)
 * 'lkup': lookup table for translating BString bytes to scoring matrix
 *         indices (integer vector)
 * 'gap_cost': gap cost or penalty (integer vector of length 1)
 * 'gap_code': encoded value of the '-' letter (raw vector of length 1)
 * Return a list of 3 elements: 2 "externalptr" objects describing the
 * alignments + the score.
 * Note that the 2 BString objects to align should contain no gaps.
 */
SEXP align_needwunsQS(SEXP s1_xp, SEXP s1_offset, SEXP s1_length,
		SEXP s2_xp, SEXP s2_offset, SEXP s2_length,
		SEXP mat, SEXP mat_nrow, SEXP lkup,
		SEXP gap_cost, SEXP gap_code)
{
	int s1_off, s1_len, s2_off, s2_len, nrow,
	    al_length, score;
	const Rbyte *s1, *s2;
	char *al1, *al2;
	SEXP ans, names, tag, ans_elt;

        s1_off = INTEGER(s1_offset)[0];
        s1_len = INTEGER(s1_length)[0];
        s1 = RAW(R_ExternalPtrTag(s1_xp)) + s1_off;
        s2_off = INTEGER(s2_offset)[0];
        s2_len = INTEGER(s2_length)[0];
        s2 = RAW(R_ExternalPtrTag(s2_xp)) + s2_off;
	nrow = INTEGER(mat_nrow)[0];
	score = needwunsQS(&al1, &al2, &al_length,
		   (char *) s1, s1_len, (char *) s2, s2_len,
		   INTEGER(mat), nrow, INTEGER(lkup), LENGTH(lkup),
		   INTEGER(gap_cost)[0], (char) RAW(gap_code)[0]);

	PROTECT(ans = NEW_LIST(3));
	/* set the names */
	PROTECT(names = allocVector(STRSXP, 3));
	SET_STRING_ELT(names, 0, mkChar("al1"));
	SET_STRING_ELT(names, 1, mkChar("al2"));
	SET_STRING_ELT(names, 2, mkChar("score"));
	SET_NAMES(ans, names);
	UNPROTECT(1);
	/* set the "al1" element */
	PROTECT(ans_elt = R_MakeExternalPtr(NULL, R_NilValue, R_NilValue));
	PROTECT(tag = NEW_RAW(al_length));
	memcpy(RAW(tag), al1, al_length * sizeof(Rbyte));
	R_SetExternalPtrTag(ans_elt, tag);
	UNPROTECT(1);
	SET_ELEMENT(ans, 0, ans_elt);
	UNPROTECT(1);
	/* set the "al2" element */
	PROTECT(ans_elt = R_MakeExternalPtr(NULL, R_NilValue, R_NilValue));
	PROTECT(tag = NEW_RAW(al_length));
	memcpy(RAW(tag), al2, al_length * sizeof(Rbyte));
	R_SetExternalPtrTag(ans_elt, tag);
	UNPROTECT(1);
	SET_ELEMENT(ans, 1, ans_elt);
	UNPROTECT(1);
	/* set the "score" element */
	PROTECT(ans_elt = NEW_INTEGER(1));
	INTEGER(ans_elt)[0] = score;
	SET_ELEMENT(ans, 2, ans_elt);
	UNPROTECT(1);
	/* ans is ready */
	UNPROTECT(1);
	return ans;
}

