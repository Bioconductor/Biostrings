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
		char c1, char * buf1,
		char c2, char * buf2,
		char c3, char * buf3,
		char c4)
{
	int cc1, cc2, cc3, i1, i2, last_nonbase_pos;
	char c;

	/* Rprintf("nS=%d nP=%d c1=%d c2=%d c3=%d c4=%d\n", nS, nP, c1, c2, c3, c4); */
	cc1 = cc2 = cc3 = 0;
	last_nonbase_pos = -1;
	for (i1 = -nP + 1, i2 = 0; i2 < nS; i1++, i2++) {
		c = S[i2];
		if (c == c1) cc1++;
		else if (c == c2) cc2++;
		else if (c == c3) cc3++;
		else if (c != c4) {
			last_nonbase_pos = i2;
			cc1 = cc2 = cc3 = 0;
		}
		if (i1 < 0)
			continue;
		if (i1 <= last_nonbase_pos) {
			buf1[i1] = buf2[i1] = buf3[i1] = 255;
			continue;
		}
		if (i1 >= 1) {
			c = S[i1 - 1];
			if (c == c1) cc1--;
			else if (c == c2) cc2--;
			else if (c == c3) cc3--;
		}
		buf1[i1] = cc1;
		buf2[i1] = cc2;
		buf3[i1] = cc3;
	}
	return;
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
	SEXP buf1, buf2, buf3;

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
	
	BOC_preprocess((char *) subj, subj_length, pat_length,
			(char) c1, (char *) RAW(buf1),
			(char) c2, (char *) RAW(buf2),
			(char) c3, (char *) RAW(buf3),
			(char) c4);
	return R_NilValue;
}

