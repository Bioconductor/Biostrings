/****************************************************************************
 *       Base Occurence Count algorithm for exact and fuzzy matching        *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/****************************************************************************/
static int debug = 0;

SEXP match_BOC_debug()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'match_BOC.c'\n", debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'match_BOC.c'\n");
#endif
	return R_NilValue;
}


/****************************************************************************
 * Initialisation of the 3 OCbuffers ("occurence count buffers")
 * =============================================================
 * nP is assumed to be >= 1 and <= min(254, nS)
 * The 3 buffers are assumed to be long enough for the preprocessing i.e. to
 * have a length >= nS - nP + 1
 */

static void BOC_preprocess(const char *S, int nS, int nP,
		char c1, char *buf1,
		char c2, char *buf2,
		char c3, char *buf3,
		char c4, double *means)
{
	int cc1, cc2, cc3, n1, n2, last_nonbase_pos,
	    total, i, partsum1, partsum2, partsum3;
	char c;

	/* Rprintf("nS=%d nP=%d c1=%d c2=%d c3=%d c4=%d\n", nS, nP, c1, c2, c3, c4); */
	cc1 = cc2 = cc3 = total = i = partsum1 = partsum2 = partsum3 = 0;
	last_nonbase_pos = -1;
	means[0] = means[1] = means[2] = 0.0;
	for (n1 = -nP + 1, n2 = 0; n2 < nS; n1++, n2++) {
		c = S[n2];
		if (c == c1) cc1++;
		else if (c == c2) cc2++;
		else if (c == c3) cc3++;
		else if (c != c4) {
			last_nonbase_pos = n2;
			cc1 = cc2 = cc3 = 0;
		}
		if (n1 < 0)
			continue;
		if (n1 <= last_nonbase_pos) {
			buf1[n1] = buf2[n1] = buf3[n1] = 255;
			continue;
		}
		if (n1 >= 1) {
			c = S[n1 - 1];
			if (c == c1) cc1--;
			else if (c == c2) cc2--;
			else if (c == c3) cc3--;
		}
		total++;
		partsum1 += buf1[n1] = cc1;
		partsum2 += buf2[n1] = cc2;
		partsum3 += buf3[n1] = cc3;
		if (i++ < 100000000)
			continue;
		i = 0;
		means[0] += partsum1;
		means[1] += partsum2;
		means[2] += partsum3;
		partsum1 = partsum2 = partsum3 = 0;
	}
	means[0] += partsum1;
	means[1] += partsum2;
	means[2] += partsum3;
	means[0] /= total;
	means[1] /= total;
	means[2] /= total;
	means[3] = nP - means[0] - means[1] - means[2];
	return;
}


/****************************************************************************
 * An implementation of the "BOC" algo for exact matching
 * ======================================================
 */

static int BOC_exact_search(const char *P, int nP, const char *S, int nS,
		char c1, char * buf1,
		char c2, char * buf2,
		char c3, char * buf3,
		char c4, double *means, int is_count_only)
{
	int count = 0, n1, n2, cc1, cc2, cc3;
	char c;

	cc1 = cc2 = cc3 = 0;
	for (n2 = 0; n2 < nP; n2++) {
		c = P[n2];
		if (c == c1) cc1++;
		else if (c == c2) cc2++;
		else if (c == c3) cc3++;
		else if (c != c4)
			error("'pattern' contains non-base DNA letters");
	}
	return count;
}


/****************************************************************************
 * .Call entry point: "match_BOC_preprocess"
 *
 * Arguments:
 *   's_xp': subject@data@xp
 *   's_offset': subject@offset
 *   's_length': subject@length
 *   'p_length': pattern_length
 *   'code1': base1_code
 *   'buf1_xp': base1_OCbuffer@xp
 *   'code2': base2_code
 *   'buf2_xp': base2_OCbuffer@xp
 *   'code3': base3_code
 *   'buf3_xp': base3_OCbuffer@xp
 *   'code4': base4_code
 ****************************************************************************/

SEXP match_BOC_preprocess(SEXP s_xp, SEXP s_offset, SEXP s_length,
		SEXP p_length,
		SEXP code1, SEXP buf1_xp,
		SEXP code2, SEXP buf2_xp,
		SEXP code3, SEXP buf3_xp,
		SEXP code4)
{
	int subj_offset, subj_length, pat_length, c1, c2, c3, c4;
	const Rbyte *subj;
	SEXP buf1, buf2, buf3, ans;

	subj_offset = INTEGER(s_offset)[0];
	subj_length = INTEGER(s_length)[0];
	subj = RAW(R_ExternalPtrTag(s_xp)) + subj_offset;
	pat_length = INTEGER(p_length)[0];
	c1 = INTEGER(code1)[0];
	buf1 = R_ExternalPtrTag(buf1_xp);
	c2 = INTEGER(code2)[0];
	buf2 = R_ExternalPtrTag(buf2_xp);
	c3 = INTEGER(code3)[0];
	buf3 = R_ExternalPtrTag(buf3_xp);
	c4 = INTEGER(code4)[0];
	
	PROTECT(ans = NEW_NUMERIC(4));
	BOC_preprocess((char *) subj, subj_length, pat_length,
			(char) c1, (char *) RAW(buf1),
			(char) c2, (char *) RAW(buf2),
			(char) c3, (char *) RAW(buf3),
			(char) c4, REAL(ans));
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * .Call entry point: "match_BOC_exact"
 * 
 * Arguments:
 *   'p_xp': pattern@data@xp
 *   'p_offset': pattern@offset
 *   'p_length': pattern@length
 *   's_xp': boc_subject@subject@data@xp
 *   's_offset': boc_subject@subject@offset
 *   's_length': boc_subject@subject@length
 *   'code1': boc_subject@base1_code
 *   'buf1_xp': boc_subject@base1_OCbuffer@xp
 *   'code2': boc_subject@base2_code
 *   'buf2_xp': boc_subject@base2_OCbuffer@xp
 *   'code3': boc_subject@base3_code
 *   'buf3_xp': boc_subject@base3_OCbuffer@xp
 *   'code4': boc_subject@base4_code
 *   'means': boc_subject@OCmeans
 *   'count_only': single logical
 * 
 ****************************************************************************/

SEXP match_BOC_exact(SEXP p_xp, SEXP p_offset, SEXP p_length,
		SEXP s_xp, SEXP s_offset, SEXP s_length,
		SEXP code1, SEXP buf1_xp,
		SEXP code2, SEXP buf2_xp,
		SEXP code3, SEXP buf3_xp,
		SEXP code4, SEXP means, SEXP count_only)
{
	int pat_offset, pat_length, subj_offset, subj_length,
	    c1, c2, c3, c4, is_count_only, count;
	const Rbyte *pat, *subj;
	SEXP buf1, buf2, buf3, ans;

	pat_offset = INTEGER(p_offset)[0];
	pat_length = INTEGER(p_length)[0];
	pat = RAW(R_ExternalPtrTag(p_xp)) + pat_offset;
	subj_offset = INTEGER(s_offset)[0];
	subj_length = INTEGER(s_length)[0];
	subj = RAW(R_ExternalPtrTag(s_xp)) + subj_offset;
	c1 = INTEGER(code1)[0];
	buf1 = R_ExternalPtrTag(buf1_xp);
	c2 = INTEGER(code2)[0];
	buf2 = R_ExternalPtrTag(buf2_xp);
	c3 = INTEGER(code3)[0];
	buf3 = R_ExternalPtrTag(buf3_xp);
	c4 = INTEGER(code4)[0];
	is_count_only = LOGICAL(count_only)[0];

	if (!is_count_only)
		_Biostrings_reset_views_buffer();
	count = BOC_exact_search(
			(char *) pat, pat_length,
			(char *) subj, subj_length,
			(char) c1, (char *) RAW(buf1),
			(char) c2, (char *) RAW(buf2),
			(char) c3, (char *) RAW(buf3),
			(char) c4, REAL(means), is_count_only);

	if (!is_count_only) {
		PROTECT(ans = allocVector(INTSXP, count));
		memcpy(INTEGER(ans), _Biostrings_get_views_start(),
					sizeof(int) * count);
        } else {
		PROTECT(ans = allocVector(INTSXP, 1));
		INTEGER(ans)[0] = count;
	}
	UNPROTECT(1);
	return ans;
}

