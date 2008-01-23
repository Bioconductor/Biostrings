/****************************************************************************
 *                                Version 2                                 *
 *                                 of the                                   *
 *      Base Occurence Count algorithm for exact and inexact matching       *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"
#include <S.h> /* for Salloc() */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/****************************************************************************/
static int debug = 0;

SEXP match_BOC2_debug()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'match_BOC2.c'\n", debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'match_BOC2.c'\n");
#endif
	return R_NilValue;
}

static int make_32bit_signature(int c1_oc, int c2_oc, int c3_oc, char pre4)
{
	int signature = 0;

	signature += c1_oc;
	signature <<= 8;
	signature += c2_oc;
	signature <<= 8;
	signature += c3_oc;
	signature <<= 8;
	signature += (unsigned char) pre4;
	return signature;
}

static char make_pre4(const char *s, char c1, char c2, char c3, char c4)
{
	char pre4, c, twobit_code;
	int i;

	for (i = 0; i < 4; i++, s++) {
		c = *s;
		if (c == c1) twobit_code = 0;
		else if (c == c2) twobit_code = 1;
		else if (c == c3) twobit_code = 2;
		else twobit_code = 3;
		pre4 <<= 2;
		pre4 += twobit_code;
	}
	return pre4;
}


/****************************************************************************
 * Initialization of the BOC2_SubjectString object
 * ===============================================
 * - nP is assumed to be >= 1 and <= min(254, nS).
 * - The buffer is assumed to be long enough for the preprocessing i.e. to
 *   have a length >= nS - nP + 1.
 * - The 4 tables are of length nP + 1.
 */

static void BOC2_preprocess(const char *S, int nS, int nP,
		char c1, char c2, char c3, char c4,
		int *buf, double *means,
		int *table1, int *table2, int *table3, int *table4)
{
	int c1_oc, c2_oc, c3_oc, n1, n2, last_nonbase_pos,
	    total, i, partsum1, partsum2, partsum3;
	char c, pre4;

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
			buf[n1] = make_32bit_signature(255, 255, 255, 0);
			continue;
		}
		if (n1 >= 1) {
			c = S[n1 - 1];
			if (c == c1) c1_oc--;
			else if (c == c2) c2_oc--;
			else if (c == c3) c3_oc--;
		}
		total++;
		pre4 =  make_pre4(S + n1, c1, c2, c3, c4);
		buf[n1] = make_32bit_signature(c1_oc, c2_oc, c3_oc, pre4);
		partsum1 += c1_oc;
		partsum2 += c2_oc;
		partsum3 += c3_oc;
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
 * An implementation of the "BOC2" algo for exact matching
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

static int split4_offsets(char codes[4], int *offsets[4], int noffsets[4], const char *P, int nP)
{
	int tmp_codes[4], *tmp_offsets[4], tmp_noffsets[4], order[4], i, offset, ii, j;
	char c;

	for (i = 0; i < 4; i++)
		noffsets[i] = 0;
	for (offset = 0; offset < nP; offset++) {
		c = P[offset];
		for (i = 0; i < 4; i++) {
			if (c != codes[i])
				continue;
			offsets[i][noffsets[i]++] = offset;
			goto continue0;
		}
		return 1;
		continue0: ;
	}
	order3(order, noffsets);
	for (i = 3; i >= 1; i--) {
		if (noffsets[3] >= noffsets[order[i-1]])
			break;
		order[i] = order[i-1];
	}
	order[i] = 3;
	// Would be neat to perform in-place reordering. Something like this:
	//   http://www.priorartdatabase.com/IPCOM/000104324/
	for (i = 0; i < 4; i++) {
		tmp_codes[i] = codes[i];
		tmp_offsets[i] = offsets[i];
		tmp_noffsets[i] = noffsets[i];
	}
	for (i = 0; i < 4; i++) {
		ii = order[i];
		codes[i] = tmp_codes[ii];
		offsets[i] = tmp_offsets[ii];
		noffsets[i] = tmp_noffsets[ii];
#ifdef DEBUG_BIOSTRINGS
		if (debug) {
			Rprintf("[DEBUG] split4_offsets: codes[%d]=%d\n", i, codes[i]);
			Rprintf("[DEBUG] split4_offsets: noffsets[%d]=%d\n", i, noffsets[i]);
			Rprintf("[DEBUG] split4_offsets: offsets[%d]=", i);
			for (j = 0; j < noffsets[i]; j++)
				Rprintf(" %d", offsets[i][j]);
			Rprintf("\n");
		}
#endif
	}
	return 0;
}

static void BOC2_exact_search(const char *P, int nP, const char *S, int nS,
		char c1, char c2, char c3, char c4,
		const int *buf, const double *means,
		const int *table1, const int *table2, const int *table3, const int *table4)
{
	int n1, n1max, n2, c1_oc, c2_oc, c3_oc, Psignature,
	    nPsuf4, *Psuf4_offsets[4], Psuf4_noffsets[4], i, j, *offsets, noffsets;
	char c, Ppre4, codes[4];
	const char *Psuf4, *Ssuf4;
#ifdef DEBUG_BIOSTRINGS
	int count_preapprovals = 0;
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
	Ppre4 = make_pre4(P, c1, c2, c3, c4);
	Psignature = make_32bit_signature(c1_oc, c2_oc, c3_oc, Ppre4);
#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] pattern: c1_oc=%d c2_oc=%d c3_oc=%d Ppre4=%d\n",
			c1_oc, c2_oc, c3_oc, Ppre4);
#endif
	Psuf4 = P + 4;
	nPsuf4 = nP - 4;
	codes[0] = c1;
	codes[1] = c2;
	codes[2] = c3;
	codes[3] = c4;
	for (i = 0; i < 4; i++)
		Psuf4_offsets[i] = Salloc((long) nP, int);
	split4_offsets(codes, Psuf4_offsets, Psuf4_noffsets, Psuf4, nPsuf4);
	n1max = nS - nP;
	for (n1 = 0, Ssuf4 = S + 4; n1 <= n1max; n1++, Ssuf4++, buf++) {
		if (Psignature != *buf)
			continue;
#ifdef DEBUG_BIOSTRINGS
		count_preapprovals++;
#endif
		if (memcmp(Psuf4, Ssuf4, nPsuf4) != 0)
			continue; // same as goto continue0;
/*
		// Uncomment the 2 lines above if you want to use the fancy
		// comparison method below.
		for (i = 0; i < 3; i++) {
			c = codes[i];
			offsets = Psuf4_offsets[i];
			noffsets = Psuf4_noffsets[i];
			for (j = 0; j < noffsets; j++)
				if (c != Ssuf4[offsets[j]])
					goto continue0;
		}
*/
		_Biostrings_report_match(n1, 0);
		continue0: ;
	}
#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] count_preapprovals=%d\n", count_preapprovals);
#endif
	return;
}


/****************************************************************************
 * .Call entry point: "match_BOC2_preprocess"
 *
 * Arguments:
 *   's_xp': subject@data@xp
 *   's_offset': subject@offset
 *   's_length': subject@length
 *   'p_length': pattern_length
 *   'code1': base1_code
 *   'code2': base2_code
 *   'code3': base3_code
 *   'code4': base4_code
 *   'buf_xp': buffer@xp
 *
 * Return an R list with the following elements:
 *   - means: atomic vector of 4 doubles
 *   - table1, table2, table3, table4: atomic vectors of (p_length + 1) integers
 *
 ****************************************************************************/

SEXP match_BOC2_preprocess(SEXP s_xp, SEXP s_offset, SEXP s_length,
		SEXP p_length,
		SEXP code1, SEXP code2, SEXP code3, SEXP code4,
		SEXP buf_xp)
{
	int subj_offset, subj_length, pat_length, c1, c2, c3, c4;
	const Rbyte *subj;
	SEXP buf, ans, ans_names, ans_elt;

	subj_offset = INTEGER(s_offset)[0];
	subj_length = INTEGER(s_length)[0];
	subj = RAW(R_ExternalPtrTag(s_xp)) + subj_offset;
	pat_length = INTEGER(p_length)[0];
	c1 = INTEGER(code1)[0];
	c2 = INTEGER(code2)[0];
	c3 = INTEGER(code3)[0];
	c4 = INTEGER(code4)[0];
	buf = R_ExternalPtrTag(buf_xp);

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

	BOC2_preprocess((char *) subj, subj_length, pat_length,
			(char) c1, (char) c2, (char) c3, (char) c4,
			INTEGER(buf),
			REAL(VECTOR_ELT(ans, 0)),
			INTEGER(VECTOR_ELT(ans, 1)),
			INTEGER(VECTOR_ELT(ans, 2)),
			INTEGER(VECTOR_ELT(ans, 3)),
			INTEGER(VECTOR_ELT(ans, 4)));

	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * .Call entry point: "match_BOC2_exact"
 * 
 * Arguments:
 *   'p_xp': pattern@data@xp
 *   'p_offset': pattern@offset
 *   'p_length': pattern@length
 *   's_xp': boc_subject@subject@data@xp
 *   's_offset': boc_subject@subject@offset
 *   's_length': boc_subject@subject@length
 *   'code1': boc_subject@base1_code
 *   'code2': boc_subject@base2_code
 *   'code3': boc_subject@base3_code
 *   'code4': boc_subject@base4_code
 *   'buf_xp': boc_subject@buffer@xp
 *   'stats': boc_subject@stats
 *   'count_only': single logical
 * 
 ****************************************************************************/

SEXP match_BOC2_exact(SEXP p_xp, SEXP p_offset, SEXP p_length,
		SEXP s_xp, SEXP s_offset, SEXP s_length,
		SEXP code1, SEXP code2, SEXP code3, SEXP code4,
		SEXP buf_xp,
		SEXP stats, SEXP count_only)
{
	int pat_offset, pat_length, subj_offset, subj_length,
	    c1, c2, c3, c4, is_count_only;
	const Rbyte *pat, *subj;
	SEXP buf, ans;

	pat_offset = INTEGER(p_offset)[0];
	pat_length = INTEGER(p_length)[0];
	pat = RAW(R_ExternalPtrTag(p_xp)) + pat_offset;
	subj_offset = INTEGER(s_offset)[0];
	subj_length = INTEGER(s_length)[0];
	subj = RAW(R_ExternalPtrTag(s_xp)) + subj_offset;
	c1 = INTEGER(code1)[0];
	c2 = INTEGER(code2)[0];
	c3 = INTEGER(code3)[0];
	c4 = INTEGER(code4)[0];
	buf = R_ExternalPtrTag(buf_xp);
	is_count_only = LOGICAL(count_only)[0];

	_Biostrings_reset_viewsbuf(is_count_only ? 1 : 2);
	BOC2_exact_search(
		(char *) pat, pat_length,
		(char *) subj, subj_length,
		(char) c1, (char) c2, (char) c3, (char) c4,
		INTEGER(buf),
		REAL(VECTOR_ELT(stats, 0)),
		INTEGER(VECTOR_ELT(stats, 1)),
		INTEGER(VECTOR_ELT(stats, 2)),
		INTEGER(VECTOR_ELT(stats, 3)),
		INTEGER(VECTOR_ELT(stats, 4)));
	if (is_count_only)
		PROTECT(ans = _Biostrings_viewsbuf_count_asINTEGER());
	else 
		PROTECT(ans = _Biostrings_viewsbuf_start_asINTEGER());
	UNPROTECT(1);
	return ans;
}

