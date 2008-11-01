/****************************************************************************
 *            An edit distance implementation with early bailout            *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int debug = 0;

SEXP debug_nedit_at()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in file %s\n",
		debug ? "on" : "off", __FILE__);
#else
	Rprintf("Debug mode not available in file %s\n", __FILE__);
#endif
	return R_NilValue;
}

/*
 * TODO: (maybe) replace static alloc of buffers by dynamic alloc.
 */
#define MAX_NEDIT 20
#define NEDIT_BUF_LENGTH (2*MAX_NEDIT+1)

static int nedit_buf1[NEDIT_BUF_LENGTH], nedit_buf2[NEDIT_BUF_LENGTH];

#define SWAP_NEDIT_BUFS(buf1, buf2) {int *tmp; tmp = (buf1); (buf1) = (buf2); (buf2) = tmp;}
#define PROPAGATE_NEDIT(buf2, B, buf1, S, j, Pc, Bmax) \
{ \
	int nedit, B2, nedit2; \
	nedit = (buf1)[(B)] + ((j) < 0 || (j) >= (S)->nelt || (S)->elts[(j)] != (Pc)); \
	if ((B2 = (B) - 1) >= 0 && (nedit2 = (buf2)[B2] + 1) < nedit) \
		nedit = nedit2; \
	if ((B2 = (B) + 1) <= (Bmax) && (nedit2 = (buf1)[B2] + 1) < nedit) \
		nedit = nedit2; \
	(buf2)[(B)] = nedit; \
}

static void print_buf2(const int *buf2, int Bmin, int Bmax)
{
	int B;

	for (B = 0; B <= Bmax; B++) {
		if (B < Bmin)
			Rprintf("%3s", "");
		else
			Rprintf("%3d", buf2[B]);
	}
	Rprintf("\n");
	return;
}

int _nedit_for_Pfirst_shift(const RoSeq *P, const RoSeq *S,
		int Pfirst_shift, int max_nedit, int loose_Pfirst_shift)
{
	int *buf1, *buf2, Bmax, B, b, i, iplus1, j, buf2_min;
	char Pc;

	if (P == NULL || P->nelt == 0)
		return 0;
	if (max_nedit > MAX_NEDIT)
		error("'max.nedit' too big");
	buf1 = nedit_buf1;
	buf2 = nedit_buf2;
	Bmax = 2 * max_nedit;
	for (B = max_nedit, b = 0; B <= Bmax; B++, b++)
		buf2[B] = b;
#ifdef DEBUG_BIOSTRINGS
	if (debug) print_buf2(buf2, max_nedit, Bmax);
#endif
	for (i = 0, iplus1 = 1; i < max_nedit; i++, iplus1++) {
		if (i >= P->nelt)
			return buf2_min;
		Pc = P->elts[i];
		SWAP_NEDIT_BUFS(buf1, buf2);
		B = max_nedit - iplus1;
		buf2_min = buf2[B++] = iplus1;
		for (j = Pfirst_shift; B <= Bmax; B++, j++) {
			PROPAGATE_NEDIT(buf2, B, buf1, S, j, Pc, Bmax);
			if (buf2[B] < buf2_min)
				buf2_min = buf2[B];
		}
#ifdef DEBUG_BIOSTRINGS
		if (debug) print_buf2(buf2, max_nedit - iplus1, Bmax);
#endif
		if (buf2_min > max_nedit)
			return buf2_min;
	}
	for ( ; i < P->nelt; i++) {
		Pc = P->elts[i];
		SWAP_NEDIT_BUFS(buf1, buf2);
		buf2_min = max_nedit + 1;
		for (B = 0, j = Pfirst_shift + i - max_nedit; B <= Bmax; B++, j++) {
			PROPAGATE_NEDIT(buf2, B, buf1, S, j, Pc, Bmax);
			if (buf2[B] < buf2_min)
				buf2_min = buf2[B];
		}
#ifdef DEBUG_BIOSTRINGS
		if (debug) print_buf2(buf2, 0, Bmax);
#endif
		if (buf2_min > max_nedit)
			return buf2_min;
	}
	return buf2_min;
}

/*
 * --- .Call ENTRY POINT ---
 */
SEXP nedit_starting_at(SEXP pattern, SEXP subject, SEXP at, SEXP fixed)
{
	RoSeq P, S;
	int at_len, fixedP, fixedS, i, *at_elt, *ans_elt, Pfirst_shift;
	SEXP ans;

	P = _get_XString_asRoSeq(pattern);
	S = _get_XString_asRoSeq(subject);
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
		Pfirst_shift = *at_elt - 1;
 		*ans_elt = _nedit_for_Pfirst_shift(&P, &S, Pfirst_shift, P.nelt, 1);
	}
	UNPROTECT(1);
	return ans;
}

