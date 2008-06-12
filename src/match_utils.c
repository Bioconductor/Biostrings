/****************************************************************************
 *              Utility functions related to pattern matching               *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int debug = 0;

SEXP debug_match_utils()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'match_utils.c'\n",
                 debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'match_utils.c'\n");
#endif
	return R_NilValue;
}

/*
 * fixedP | fixedS | letters p and s match iff...
 * ----------------------------------------------------------
 * TRUE   | TRUE   | ...they are equal
 * TRUE   | FALSE  | ...bits at 1 in p are are also at 1 in s
 * FALSE  | TRUE   | ...bits at 1 in s are also at 1 in p
 * FALSE  | FALSE  | ...they share at least one bit at 1
 */
#define LETTERS_MATCH(p, s, fixedP, fixedS, OK) \
{ \
	if (fixedP) { \
		if (fixedS) \
			OK = (p) == (s); \
		else \
			OK = ((p) & ~(s)) == 0; \
	} else { \
		if (fixedS) \
			OK = (~(p) & (s)) == 0; \
		else \
			OK = (p) & (s); \
	} \
}


/****************************************************************************
 * nmismatch_at()
 */

int _nmismatch_at_Pshift(const RoSeq *P, const RoSeq *S, int Pshift,
		int fixedP, int fixedS)
{
	int nmismatch, i, j, OK;
	const char *p, *s;

	nmismatch = 0;
	for (i = 0, j = Pshift, p = P->elts, s = S->elts + Pshift;
	     i < P->nelt;
	     i++, j++, p++, s++)
	{
		if (j < 0 || S->nelt <= j) {
			nmismatch++;
			continue;
		}
		LETTERS_MATCH(*p, *s, fixedP, fixedS, OK);
		if (!OK)
			nmismatch++;
	}
	return nmismatch;
}

/*
 * --- .Call ENTRY POINT ---
 */
SEXP nmismatch_at(SEXP pattern, SEXP subject, SEXP starting, SEXP at, SEXP fixed)
{
	RoSeq P, S;
	int is_start, at_len, fixedP, fixedS, i, *at_elt, *ans_elt, Pshift;
	SEXP ans;

	P = _get_XString_asRoSeq(pattern);
	S = _get_XString_asRoSeq(subject);
	is_start = LOGICAL(starting)[0];
	at_len = LENGTH(at);
	fixedP = LOGICAL(fixed)[0];
	fixedS = LOGICAL(fixed)[1];

	PROTECT(ans = NEW_INTEGER(at_len));
	for (i = 0, at_elt = INTEGER(at), ans_elt = INTEGER(ans);
             i < at_len;
             i++, at_elt++, ans_elt++)
	{
		if (*at_elt == NA_INTEGER) {
			*ans_elt = NA_INTEGER;
			continue;
		}
		Pshift = *at_elt - 1;
		if (is_start)
			Pshift = *at_elt - 1;
		else
			Pshift = *at_elt - P.nelt;
 		*ans_elt = _nmismatch_at_Pshift(&P, &S, Pshift, fixedP, fixedS);
	}
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * is_matching()
 */

/* max_mm must be >= 0 (not safe otherwise) */
int _is_matching_at_Pshift(const RoSeq *P, const RoSeq *S, int Pshift,
		int max_mm, int fixedP, int fixedS)
{
	const char *p, *s;
	int plen, slen, min_pm, mm, pm, i, OK;

#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] _is_matching_at_Pshift(): "
			"P->nelt=%d S->nelt=%d Pshift=%d\n",
			P->nelt, S->nelt, Pshift);
#endif
	if (P->nelt <= max_mm)
		return 1;
	p = P->elts;
	plen = P->nelt;
	s = S->elts;
	slen = S->nelt;
	// 0 <= max_mm < plen
	if (Pshift < 0) {
		max_mm += Pshift;
		// Pshift <= max_mm < plen + Pshift < plen
		if (max_mm < 0)
			return 0;
		// -plen < Pshift < 0
		p -= Pshift;
		plen += Pshift; // 0 < plen
	} else {
		s += Pshift;
		slen -= Pshift;
	}
	if (plen > slen) {
		max_mm -= plen - slen;
		if (max_mm < 0)
			return 0;
		plen = slen;
	}
	min_pm = plen - max_mm;
	mm = pm = 0;
	// 0 = mm <= max_mm < plen <= slen
	// 0 = pm < min_pm <= plen <= slen
	// min_pm + max_mm = plen
	for (i = 0; i < plen; i++, p++, s++) {
		LETTERS_MATCH(*p, *s, fixedP, fixedS, OK);
		if (OK) {
			if (++pm >= min_pm)
				return 1;
		} else {
			if (++mm > max_mm)
				return 0;
		}
	}
	error("Biostrings internal error in _is_matching_at_Pshift(): "
	      "we should never be here");
	return -1;
}

/*
 * --- .Call ENTRY POINT ---
 * is_matching() arguments are assumed to be:
 *   pattern: pattern
 *   subject: subject
 *   algorithm: algorithm
 *   start: starting positions of the pattern relative to the subject
 *   max_mismatch: the max number of mismatching letters
 *   fixed: logical vector of length 2
 * Return a logical vector of the same length as 'start'.
 */
SEXP is_matching(SEXP pattern, SEXP subject, SEXP start,
		SEXP max_mismatch, SEXP fixed)
{
	RoSeq P, S;
	int start_len, max_mm, fixedP, fixedS, i, *start_elt, *ans_elt;
	SEXP ans;

	P = _get_XString_asRoSeq(pattern);
	S = _get_XString_asRoSeq(subject);
	start_len = LENGTH(start);
	max_mm = INTEGER(max_mismatch)[0];
	fixedP = LOGICAL(fixed)[0];
	fixedS = LOGICAL(fixed)[1];

	PROTECT(ans = NEW_LOGICAL(start_len));
	for (i = 0, start_elt = INTEGER(start), ans_elt = LOGICAL(ans);
             i < start_len;
             i++, start_elt++, ans_elt++)
	{
		if (*start_elt == NA_INTEGER) {
			*ans_elt = NA_LOGICAL;
			continue;
		}
 		*ans_elt = _is_matching_at_Pshift(&P, &S, *start_elt - 1,
				max_mm, fixedP, fixedS);
	}
	UNPROTECT(1);
	return ans;
}

