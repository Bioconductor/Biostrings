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


/****************************************************************************
 * nmismatch_at()
 *
 * The 4 static functions below stop counting mismatches if the number
 * exceeds 'max_mm'. The caller can disable this by passing 'P->nelt' to
 * the 'max_mm' arg.
 *
 * fixedP | fixedS | letters *p and *s match iff...
 * --------------------------------------------------------
 * TRUE   | TRUE   | ...they are equal
 * TRUE   | FALSE  | ...bits at 1 in *p are also at 1 in *s
 * FALSE  | TRUE   | ...bits at 1 in *s are also at 1 in *p
 * FALSE  | FALSE  | ...they share at least one bit at 1
 */

static int nmismatch_at_Pshift_fixedPfixedS(const RoSeq *P, const RoSeq *S,
		int Pshift, int max_mm)
{
	int nmismatch, i, j;
	const char *p, *s;

	if (P == NULL)
		return 0;
	nmismatch = 0;
	for (i = 0, j = Pshift, p = P->elts, s = S->elts + Pshift;
	     i < P->nelt;
	     i++, j++, p++, s++)
	{
		if (j >= 0 && j < S->nelt && *p == *s)
			continue;
		if (nmismatch++ >= max_mm)
			break;
	}
	return nmismatch;
}

static int nmismatch_at_Pshift_fixedPnonfixedS(const RoSeq *P, const RoSeq *S,
		int Pshift, int max_mm)
{
	int nmismatch, i, j;
	const char *p, *s;

	if (P == NULL)
		return 0;
	nmismatch = 0;
	for (i = 0, j = Pshift, p = P->elts, s = S->elts + Pshift;
	     i < P->nelt;
	     i++, j++, p++, s++)
	{
		if (j >= 0 && j < S->nelt && ((*p) & ~(*s)) == 0)
			continue;
		if (nmismatch++ >= max_mm)
			break;
	}
	return nmismatch;
}

static int nmismatch_at_Pshift_nonfixedPfixedS(const RoSeq *P, const RoSeq *S,
		int Pshift, int max_mm)
{
	int nmismatch, i, j;
	const char *p, *s;

	if (P == NULL)
		return 0;
	nmismatch = 0;
	for (i = 0, j = Pshift, p = P->elts, s = S->elts + Pshift;
	     i < P->nelt;
	     i++, j++, p++, s++)
	{
		if (j >= 0 && j < S->nelt && (~(*p) & (*s)) == 0)
			continue;
		if (nmismatch++ >= max_mm)
			break;
	}
	return nmismatch;
}

static int nmismatch_at_Pshift_nonfixedPnonfixedS(const RoSeq *P, const RoSeq *S,
		int Pshift, int max_mm)
{
	int nmismatch, i, j;
	const char *p, *s;

	if (P == NULL)
		return 0;
	nmismatch = 0;
	for (i = 0, j = Pshift, p = P->elts, s = S->elts + Pshift;
	     i < P->nelt;
	     i++, j++, p++, s++)
	{
		if (j >= 0 && j < S->nelt && ((*p) & (*s)))
			continue;
		if (nmismatch++ >= max_mm)
			break;
	}
	return nmismatch;
}

int (*_selected_nmismatch_at_Pshift_fun)(const RoSeq *P, const RoSeq *S,
		int Pshift, int max_mm);

void _select_nmismatch_at_Pshift_fun(int fixedP, int fixedS)
{
	if (fixedP) {
		if (fixedS)
			_selected_nmismatch_at_Pshift_fun = &nmismatch_at_Pshift_fixedPfixedS;
		else
			_selected_nmismatch_at_Pshift_fun = &nmismatch_at_Pshift_fixedPnonfixedS;
	} else {
		if (fixedS)
			_selected_nmismatch_at_Pshift_fun = &nmismatch_at_Pshift_nonfixedPfixedS;
		else
			_selected_nmismatch_at_Pshift_fun = &nmismatch_at_Pshift_nonfixedPnonfixedS;
	}
	return;
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
	_select_nmismatch_at_Pshift_fun(fixedP, fixedS);

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
 		*ans_elt = _selected_nmismatch_at_Pshift_fun(&P, &S, Pshift, P.nelt);
	}
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * is_matching()
 */

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
	_select_nmismatch_at_Pshift_fun(fixedP, fixedS);

	PROTECT(ans = NEW_LOGICAL(start_len));
	for (i = 0, start_elt = INTEGER(start), ans_elt = LOGICAL(ans);
             i < start_len;
             i++, start_elt++, ans_elt++)
	{
		if (*start_elt == NA_INTEGER) {
			*ans_elt = NA_LOGICAL;
			continue;
		}
 		*ans_elt = _selected_nmismatch_at_Pshift_fun(&P, &S,
				*start_elt - 1, max_mm) <= max_mm;
	}
	UNPROTECT(1);
	return ans;
}

