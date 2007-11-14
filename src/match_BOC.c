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
 *   - nP is assumed to be >= 1 and <= min(254, nS).
 *   - The 3 buffers are assumed to be long enough for the preprocessing i.e. to
 *     have a length >= nS - nP + 1.
 *   - The 4 tables are of length nP + 1.
 */

static void BOC_preprocess(const char *S, int nS, int nP,
		char c1, char *buf1,
		char c2, char *buf2,
		char c3, char *buf3,
		char c4,
		double *means,
		int *table1, int *table2, int *table3, int *table4)
{
	int c1_oc, c2_oc, c3_oc, n1, n2, last_nonbase_pos,
	    total, i, partsum1, partsum2, partsum3;
	char c;

	/* Rprintf("nS=%d nP=%d c1=%d c2=%d c3=%d c4=%d\n", nS, nP, c1, c2, c3, c4); */
	for (i = 0; i <= nP; i++)
		table1[i] = table2[i] = table3[i] = table4[i] = 0;
	c1_oc = c2_oc = c3_oc = total = i = partsum1 = partsum2 = partsum3 = 0;
	last_nonbase_pos = -1;
	means[0] = means[1] = means[2] = 0.0;
	for (n1 = -nP + 1, n2 = 0; n2 < nS; n1++, n2++) {
		c = S[n2];
		if (c == c1) c1_oc++;
		else if (c == c2) c2_oc++;
		else if (c == c3) c3_oc++;
		else if (c != c4) {
			last_nonbase_pos = n2;
			c1_oc = c2_oc = c3_oc = 0;
		}
		if (n1 < 0)
			continue;
		if (n1 <= last_nonbase_pos) {
			buf1[n1] = buf2[n1] = buf3[n1] = 255;
			continue;
		}
		if (n1 >= 1) {
			c = S[n1 - 1];
			if (c == c1) c1_oc--;
			else if (c == c2) c2_oc--;
			else if (c == c3) c3_oc--;
		}
		total++;
		partsum1 += buf1[n1] = c1_oc;
		partsum2 += buf2[n1] = c2_oc;
		partsum3 += buf3[n1] = c3_oc;
		table1[c1_oc]++;
		table2[c2_oc]++;
		table3[c3_oc]++;
		table4[nP - c1_oc - c2_oc - c3_oc]++;
		if (i++ < 5000000)
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

/* 'order' and 'x' must be int arrays of length 3 */
static void order3(int *order, const int *x)
{
	if (x[1] < x[0]) {
		if (x[2] < x[0]) {
			order[2] = 0;
			if (x[1] < x[2]) {
				order[0] = 1;
				order[1] = 2;
			} else {
				order[0] = 2;
				order[1] = 1;
			}
		} else {
			order[0] = 1;
			order[1] = 0;
			order[2] = 2;
		}
	} else {
		if (x[2] < x[1]) {
			order[2] = 1;
			if (x[0] < x[2]) {
				order[0] = 0;
				order[1] = 2;
			} else {
				order[0] = 2;
				order[1] = 0;
			}
		} else {
			order[0] = 0;
			order[1] = 1;
			order[2] = 2;
		}
	}
	return;
}

static int switch_oc(int i, int c1_oc, int c2_oc, int c3_oc)
{
	switch (i) {
		case 0: return c1_oc;
		case 1: return c2_oc;
		case 2: return c3_oc;
	}
}

static const char *switch_buf(int i, const char *buf1, const char *buf2, const char *buf3)
{
	switch (i) {
		case 0: return buf1;
		case 1: return buf2;
		case 2: return buf3;
	}
}

static void order_bases(char *p_ocX, char *p_ocY, char *p_ocZ,
		const char **p_bufX, const char **p_bufY, const char **p_bufZ,
		int c1_oc, int c2_oc, int c3_oc,
		const char *buf1, const char *buf2, const char *buf3,
		int *nmers)
{
	int order[3];

	order3(order, nmers);
	// To turn off reordering, uncomment the following 3 lines
	/*
	order[0] = 0;
	order[1] = 1;
	order[2] = 2;
	*/
#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] order_bases: order[0]=%d order[1]=%d order[2]=%d\n",
			order[0], order[1], order[2]);
#endif
	*p_ocX = (char) switch_oc(order[0], c1_oc, c2_oc, c3_oc);
	*p_ocY = (char) switch_oc(order[1], c1_oc, c2_oc, c3_oc);
	*p_ocZ = (char) switch_oc(order[2], c1_oc, c2_oc, c3_oc);
	*p_bufX = switch_buf(order[0], buf1, buf2, buf3);
	*p_bufY = switch_buf(order[1], buf1, buf2, buf3);
	*p_bufZ = switch_buf(order[2], buf1, buf2, buf3);
	return;
}

static int BOC_exact_search(const char *P, int nP, const char *S, int nS,
		char c1, const char *buf1,
		char c2, const char *buf2,
		char c3, const char *buf3,
		char c4,
		const double *means,
		const int *table1, const int *table2, const int *table3, const int *table4,
		int is_count_only)
{
	int count = 0, n1, n1max, n2, c1_oc, c2_oc, c3_oc, nmers[3];
	char c, ocX, ocY, ocZ;
	const char *bufX, *bufY, *bufZ;
#ifdef DEBUG_BIOSTRINGS
	int count_memcmp = 0;
#endif

	c1_oc = c2_oc = c3_oc = 0;
	for (n2 = 0; n2 < nP; n2++) {
		c = P[n2];
		if (c == c1) c1_oc++;
		else if (c == c2) c2_oc++;
		else if (c == c3) c3_oc++;
		else if (c != c4)
			error("'pattern' contains non-base DNA letters");
	}
#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] pattern: c1_oc=%d c2_oc=%d c3_oc=%d\n", c1_oc, c2_oc, c3_oc);
#endif
	nmers[0] = table1[c1_oc];
	nmers[1] = table2[c2_oc];
	nmers[2] = table3[c3_oc];
	order_bases(&ocX, &ocY, &ocZ, &bufX, &bufY, &bufZ,
			c1_oc, c2_oc, c3_oc, buf1, buf2, buf3, nmers);
	n1max = nS - nP;
	for (n1 = 0; n1 <= n1max; n1++, bufX++, bufY++, bufZ++) {
		if (ocX != *bufX)
			continue;
		if (ocY != *bufY)
			continue;
		if (ocZ != *bufZ)
			continue;
#ifdef DEBUG_BIOSTRINGS
		count_memcmp++;
#endif
		if (memcmp(P, S + n1, nP) != 0)
			continue;
		if (!is_count_only)
			_Biostrings_report_match(n1, 0);
		count++;
	}
#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] count_memcmp=%d\n", count_memcmp);
#endif
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
 *
 * Return an R list with the following elements:
 *   - means: atomic vector of 4 doubles
 *   - table1, table2, table3, table4: atomic vectors of (p_length + 1) integers
 *
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
	SEXP buf1, buf2, buf3, ans, ans_names, ans_elt;

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

	PROTECT(ans = NEW_LIST(5));
	/* set the names */
	PROTECT(ans_names = NEW_CHARACTER(5));
	SET_STRING_ELT(ans_names, 0, mkChar("means"));
	SET_STRING_ELT(ans_names, 1, mkChar("table1"));
	SET_STRING_ELT(ans_names, 2, mkChar("table2"));
	SET_STRING_ELT(ans_names, 3, mkChar("table3"));
	SET_STRING_ELT(ans_names, 4, mkChar("table4"));
	SET_NAMES(ans, ans_names);
	UNPROTECT(1);
	/* set the "means" element */
	PROTECT(ans_elt = NEW_NUMERIC(4));
	SET_ELEMENT(ans, 0, ans_elt);
	UNPROTECT(1);
	/* set the "table1" element */
	PROTECT(ans_elt = NEW_INTEGER(pat_length + 1));
	SET_ELEMENT(ans, 1, ans_elt);
	UNPROTECT(1);
	/* set the "table2" element */
	PROTECT(ans_elt = NEW_INTEGER(pat_length + 1));
	SET_ELEMENT(ans, 2, ans_elt);
	UNPROTECT(1);
	/* set the "table3" element */
	PROTECT(ans_elt = NEW_INTEGER(pat_length + 1));
	SET_ELEMENT(ans, 3, ans_elt);
	UNPROTECT(1);
	/* set the "table4" element */
	PROTECT(ans_elt = NEW_INTEGER(pat_length + 1));
	SET_ELEMENT(ans, 4, ans_elt);
	UNPROTECT(1);

	BOC_preprocess((char *) subj, subj_length, pat_length,
			(char) c1, (char *) RAW(buf1),
			(char) c2, (char *) RAW(buf2),
			(char) c3, (char *) RAW(buf3),
			(char) c4,
			REAL(VECTOR_ELT(ans, 0)),
			INTEGER(VECTOR_ELT(ans, 1)),
			INTEGER(VECTOR_ELT(ans, 2)),
			INTEGER(VECTOR_ELT(ans, 3)),
			INTEGER(VECTOR_ELT(ans, 4)));

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
 *   'stats': boc_subject@stats
 *   'count_only': single logical
 * 
 ****************************************************************************/

SEXP match_BOC_exact(SEXP p_xp, SEXP p_offset, SEXP p_length,
		SEXP s_xp, SEXP s_offset, SEXP s_length,
		SEXP code1, SEXP buf1_xp,
		SEXP code2, SEXP buf2_xp,
		SEXP code3, SEXP buf3_xp,
		SEXP code4, SEXP stats, SEXP count_only)
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
			(char) c4,
			REAL(VECTOR_ELT(stats, 0)),
			INTEGER(VECTOR_ELT(stats, 1)),
			INTEGER(VECTOR_ELT(stats, 2)),
			INTEGER(VECTOR_ELT(stats, 3)),
			INTEGER(VECTOR_ELT(stats, 4)),
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

