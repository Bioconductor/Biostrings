#include "Biostrings.h"

#include <stdio.h>


/* Return the length of the Longest Common Prefix. */
static int _lcprefix(const char *S1, int nS1, const char *S2, int nS2)
{
	int n;

	n = 0;
	while (n < nS1 && n < nS2) {
		if (*(S1++) != *(S2++))
			break;
		n++;
	}
	return n;
}

/* Return the length of the Longest Common Suffix. */
static int _lcsuffix(const char *S1, int nS1, const char *S2, int nS2)
{
	int n;

	n = 0;
	S1 += nS1 - 1;
	S2 += nS2 - 1;
	while (n < nS1 && n < nS2) {
		if (*(S1--) != *(S2--))
			break;
		n++;
	}
	return n;
}

/*
 * 's1_xp', 's1_offset', 's1_length': left BString object
 * 's2_xp', 's2_offset', 's2_length': right BString object
 */
SEXP lcprefix(SEXP s1_xp, SEXP s1_offset, SEXP s1_length,
		SEXP s2_xp, SEXP s2_offset, SEXP s2_length)
{
	int s1_off, s1_len, s2_off, s2_len, n;
	const Rbyte *s1, *s2;
	SEXP ans;

        s1_off = INTEGER(s1_offset)[0];
        s1_len = INTEGER(s1_length)[0];
        s1 = RAW(R_ExternalPtrTag(s1_xp)) + s1_off;
        s2_off = INTEGER(s2_offset)[0];
        s2_len = INTEGER(s2_length)[0];
        s2 = RAW(R_ExternalPtrTag(s2_xp)) + s2_off;

	n = _lcprefix((char *) s1, s1_len, (char *) s2, s2_len);

	PROTECT(ans = allocVector(INTSXP, 1));
	INTEGER(ans)[0] = n;
	UNPROTECT(1);
	return ans;
}

/*
 * 's1_xp', 's1_offset', 's1_length': left BString object
 * 's2_xp', 's2_offset', 's2_length': right BString object
 */
SEXP lcsuffix(SEXP s1_xp, SEXP s1_offset, SEXP s1_length,
		SEXP s2_xp, SEXP s2_offset, SEXP s2_length)
{
	int s1_off, s1_len, s2_off, s2_len, n;
	const Rbyte *s1, *s2;
	SEXP ans;

        s1_off = INTEGER(s1_offset)[0];
        s1_len = INTEGER(s1_length)[0];
        s1 = RAW(R_ExternalPtrTag(s1_xp)) + s1_off;
        s2_off = INTEGER(s2_offset)[0];
        s2_len = INTEGER(s2_length)[0];
        s2 = RAW(R_ExternalPtrTag(s2_xp)) + s2_off;

	n = _lcsuffix((char *) s1, s1_len, (char *) s2, s2_len);

	PROTECT(ans = allocVector(INTSXP, 1));
	INTEGER(ans)[0] = n;
	UNPROTECT(1);
	return ans;
}

