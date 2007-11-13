#include "Biostrings.h"

#include <stdio.h>


#define SCO(i1, i2) (sco[n2 * (i1) + (i2)])
#define TRA(i1, i2) (tra[n2 * (i1) + (i2)])

#define SET_LKUP_VAL_INT(lkup, length, key) \
{ \
	unsigned char lkup_key = (unsigned char) (key); \
	if (lkup_key >= (length) || (lkup_val = (lkup)[lkup_key]) == NA_INTEGER) { \
		error("key %d not in lookup table", (int) lkup_key); \
	} \
}

static int nal = 0;
static char *al1_buf, *al2_buf, *al1, *al2;

/* Returns the score of the alignment */
static int needwunsQS(const char *S1, int nS1, const char *S2, int nS2,
		const int *mat, int mat_nrow, const int *lkup, int lkup_length,
		int gap_cost, char gap_code)
{
	int *sco, sc, n1, n2, i1, i2, j1, j2, al_buf_size;
	char *tra, tr;

	n1 = nS1 + 1;
	n2 = nS2 + 1;
	sco = (int *) R_alloc((long) n1 * n2, sizeof(int));
	tra = (char *) R_alloc((long) n1 * n2, sizeof(char));
	for (i1 = 0; i1 <= nS1; i1++)
		SCO(i1, 0) = - i1 * gap_cost;
	for (i2 = 1; i2 <= nS2; i2++)
		SCO(0, i2) = - i2 * gap_cost;
	for (i1 = 1, j1 = 0; i1 <= nS1; i1++, j1++) {
		for (i2 = 1, j2 = 0; i2 <= nS2; i2++, j2++) {
			int lkup_val, m1, m2, scR, scD, scI;
			SET_LKUP_VAL_INT(lkup, lkup_length, S1[j1]);
			m1 = lkup_val;
			SET_LKUP_VAL_INT(lkup, lkup_length, S2[j2]);
			m2 = lkup_val;
			scR = SCO(j1, j2) + mat[mat_nrow * m1 + m2];
			scD = SCO(j1, i2) - gap_cost;
			scI = SCO(i1, j2) - gap_cost;
			if (scD >= scI) {
				sc = scD;
				tr = 'D'; /* Deletion (gap in aligned string 2) */
			} else {
				sc = scI;
				tr = 'I'; /* Insertion (gap in aligned string 1) */
			}
			if (scR >= sc) {
				sc = scR;
				tr = 'R'; /* Replacement (can be an exact match) */
			}
			SCO(i1, i2) = sc;
			TRA(i1, i2) = tr;
		}
	}
/*
#ifdef DEBUG_BIOSTRINGS
	Rprintf("sco:\n");
	for (i1 = 0; i1 <= nS1; i1++) {
		for (i2 = 0; i2 <= nS2; i2++) {
			Rprintf("%4d", SCO(i1, i2));
		}
		Rprintf("\n");
	}
	Rprintf("tra:\n");
	for (i1 = 1; i1 <= nS1; i1++) {
		for (i2 = 1; i2 <= nS2; i2++) {
			Rprintf(" %c", TRA(i1, i2));
		}
		Rprintf("\n");
	}
#endif
*/
	al_buf_size = nS1 + nS2;
	al1_buf = (char *) R_alloc((long) al_buf_size, sizeof(char));
	al2_buf = (char *) R_alloc((long) al_buf_size, sizeof(char));
	nal = 0;
	al1 = al1_buf + al_buf_size;
	al2 = al2_buf + al_buf_size;
	i1 = nS1; i2 = nS2;
	while (i1 >= 1 || i2 >= 1) {
		nal++;
		al1--;
		al2--;
		j1 = i1 - 1;
		j2 = i2 - 1;
		if (i2 == 0)
			tr = 'D';
		else if (i1 == 0)
			tr = 'I';
		else
			tr = TRA(i1, i2);
		switch (tr) {
		    case 'D':
			*al1 = S1[j1];
			*al2 = gap_code;
			i1--;
			break;
		    case 'I':
			*al1 = gap_code;
			*al2 = S2[j2];
			i2--;
			break;
		    case 'R':
			*al1 = S1[j1];
			*al2 = S2[j2];
			i1--;
			i2--;
			break;
		    default:
			error("unknown traceback code %d", (int) tr);
			break;
		}
	}
	return sc;
}

/*
 * 's1_xp', 's1_offset', 's1_length': left BString object
 * 's2_xp', 's2_offset', 's2_length': right BString object
 * 'mat': scoring matrix (integer square matrix)
 * 'lkup': lookup table for translating BString bytes to scoring matrix
 *         indices (integer vector)
 * 'gap_cost': gap cost or penalty (integer vector of length 1)
 * 'gap_code': encoded value of the '-' letter (raw vector of length 1)
 * Return a named list with 3 elements: 2 "externalptr" objects describing
 * the alignments + the score.
 * Note that the 2 BString objects to align should contain no gaps.
 */
SEXP align_needwunsQS(SEXP s1_xp, SEXP s1_offset, SEXP s1_length,
		SEXP s2_xp, SEXP s2_offset, SEXP s2_length,
		SEXP mat, SEXP mat_nrow, SEXP lkup,
		SEXP gap_cost, SEXP gap_code)
{
	int s1_off, s1_len, s2_off, s2_len, nrow, score;
	const Rbyte *s1, *s2;
	SEXP ans, ans_names, tag, ans_elt;

	s1_off = INTEGER(s1_offset)[0];
	s1_len = INTEGER(s1_length)[0];
	s1 = RAW(R_ExternalPtrTag(s1_xp)) + s1_off;
	s2_off = INTEGER(s2_offset)[0];
	s2_len = INTEGER(s2_length)[0];
	s2 = RAW(R_ExternalPtrTag(s2_xp)) + s2_off;
	nrow = INTEGER(mat_nrow)[0];
	score = needwunsQS((char *) s1, s1_len, (char *) s2, s2_len,
		   INTEGER(mat), nrow, INTEGER(lkup), LENGTH(lkup),
		   INTEGER(gap_cost)[0], (char) RAW(gap_code)[0]);

	PROTECT(ans = NEW_LIST(3));
	/* set the names */
	PROTECT(ans_names = NEW_CHARACTER(3));
	SET_STRING_ELT(ans_names, 0, mkChar("al1"));
	SET_STRING_ELT(ans_names, 1, mkChar("al2"));
	SET_STRING_ELT(ans_names, 2, mkChar("score"));
	SET_NAMES(ans, ans_names);
	UNPROTECT(1);
	/* set the "al1" element */
	PROTECT(ans_elt = R_MakeExternalPtr(NULL, R_NilValue, R_NilValue));
	PROTECT(tag = NEW_RAW(nal));
	memcpy((char *) RAW(tag), al1, nal * sizeof(char));
	R_SetExternalPtrTag(ans_elt, tag);
	UNPROTECT(1);
	SET_ELEMENT(ans, 0, ans_elt);
	UNPROTECT(1);
	/* set the "al2" element */
	PROTECT(ans_elt = R_MakeExternalPtr(NULL, R_NilValue, R_NilValue));
	PROTECT(tag = NEW_RAW(nal));
	memcpy((char *) RAW(tag), al2, nal * sizeof(char));
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

