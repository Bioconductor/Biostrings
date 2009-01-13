/****************************************************************************
 *              Utility functions related to pattern matching               *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"

static int debug = 0;

SEXP debug_match_pattern_at()
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


/****************************************************************************
 * An edit distance implementation with early bailout.
 */

/*
 * TODO: (maybe) replace static alloc of buffers by dynamic alloc.
 */
#define MAX_NEDIT 20
#define MAX_ROW_LENGTH (2*MAX_NEDIT+1)

static int row1_buf[MAX_ROW_LENGTH], row2_buf[MAX_ROW_LENGTH];

#define SWAP_NEDIT_BUFS(prev_row, curr_row) \
{ \
	int *tmp; \
	tmp = (prev_row); \
	(prev_row) = (curr_row); \
	(curr_row) = tmp; \
}

#define PROPAGATE_NEDIT(curr_row, B, prev_row, S, j, Pc, row_length) \
{ \
	int nedit, B2, nedit2; \
	nedit = (prev_row)[(B)] + ((j) < 0 || (j) >= (S)->nelt || (S)->elts[(j)] != (Pc)); \
	if ((B2 = (B) - 1) >= 0 && (nedit2 = (curr_row)[B2] + 1) < nedit) \
		nedit = nedit2; \
	if ((B2 = (B) + 1) < (row_length) && (nedit2 = (prev_row)[B2] + 1) < nedit) \
		nedit = nedit2; \
	(curr_row)[(B)] = nedit; \
}

#ifdef DEBUG_BIOSTRINGS
static void print_curr_row(const char* margin, const int *curr_row, int Bmin, int row_length)
{
	int B;

	Rprintf("[DEBUG]   %s: ", margin);
	for (B = 0; B < row_length; B++) {
		if (B < Bmin)
			Rprintf("%3s", "");
		else
			Rprintf("%3d", curr_row[B]);
	}
	Rprintf("\n");
	return;
}
#endif

/*
 * P left-offset (Ploffset) is the offset of P's first letter in S.
 * P right-offset (Proffset) is the offset of P's last letter in S.
 * The min width (min_width) is the length of the shortest substring S'
 * of S starting at Ploffset (or ending at Proffset) for which nedit(P, S')
 * is minimal.
 * TODO: Implement the 'loose_Ploffset' feature (allowing or not an indel
 * on the first letter of the local alignement).
 */
int _nedit_for_Ploffset(const RoSeq *P, const RoSeq *S, int Ploffset,
		int max_nedit, int loose_Ploffset, int *min_width)
{
	int max_nedit_plus1, *prev_row, *curr_row, row_length,
	    B, b, i, iplus1, jmin, j, min_nedit;
	char Pc;

#ifdef DEBUG_BIOSTRINGS
	if (debug) Rprintf("[DEBUG] _nedit_for_Ploffset():\n");
#endif
	if (P == NULL || P->nelt == 0)
		return 0;
	if (max_nedit == 0)
		error("Biostrings internal error in _nedit_for_Ploffset(): ",
		      "use _selected_nmismatch_at_Pshift_fun() when 'max_nedit' is 0");
	max_nedit_plus1 = max_nedit + 1;
	if (max_nedit > P->nelt)
		max_nedit = P->nelt;
	// from now max_nedit <= P->nelt
	if (max_nedit > MAX_NEDIT)
		error("'max.nedit' too big");
	prev_row = row1_buf;
	curr_row = row2_buf;
	row_length = 2 * max_nedit + 1;
	jmin = Ploffset;

	// STAGE 0:
	for (B = max_nedit, b = 0; B < row_length; B++, b++)
		curr_row[B] = b;
#ifdef DEBUG_BIOSTRINGS
	if (debug) print_curr_row("STAGE0", curr_row, max_nedit, row_length);
#endif

	// STAGE 1 (1st for() loop): no attempt is made to bailout during
	// this stage because the smallest value in curr_row is guaranteed
	// to be <= iplus1 < max_nedit.
	for (iplus1 = 1, i = 0; iplus1 < max_nedit; iplus1++, i++) {
		Pc = P->elts[i]; // i < iplus1 < max_nedit <= P->nelt
		SWAP_NEDIT_BUFS(prev_row, curr_row);
		B = max_nedit - iplus1;
		curr_row[B++] = iplus1;
		for (j = jmin; B < row_length; B++, j++)
			PROPAGATE_NEDIT(curr_row, B, prev_row, S, j, Pc, row_length);
#ifdef DEBUG_BIOSTRINGS
		if (debug) print_curr_row("STAGE1", curr_row, max_nedit - iplus1, row_length);
#endif
	}

	// STAGE 2: no attempt is made to bailout during this stage either.
	Pc = P->elts[i];
	SWAP_NEDIT_BUFS(prev_row, curr_row);
	B = 0;
	curr_row[B++] = min_nedit = iplus1;
	*min_width = 0;
	for (j = jmin; B < row_length; B++, j++) {
		PROPAGATE_NEDIT(curr_row, B, prev_row, S, j, Pc, row_length);
		if (curr_row[B] < min_nedit) {
			min_nedit = curr_row[B];
			*min_width = j - Ploffset + 1;
		}
	}
#ifdef DEBUG_BIOSTRINGS
	if (debug) print_curr_row("STAGE2", curr_row, 0, row_length);
#endif
	iplus1++;
	i++;

	// STAGE 3 (2nd for() loop): with attempt to bailout.
	for ( ; i < P->nelt; i++, iplus1++, jmin++) {
		Pc = P->elts[i];
		SWAP_NEDIT_BUFS(prev_row, curr_row);
		min_nedit = iplus1;
		*min_width = 0;
		for (B = 0, j = jmin; B < row_length; B++, j++) {
			PROPAGATE_NEDIT(curr_row, B, prev_row, S, j, Pc, row_length);
			if (curr_row[B] < min_nedit) {
				min_nedit = curr_row[B];
				*min_width = j - Ploffset + 1;
			}
		}
#ifdef DEBUG_BIOSTRINGS
		if (debug) print_curr_row("STAGE3", curr_row, 0, row_length);
#endif
		if (min_nedit >= max_nedit_plus1) // should never be min_nedit > max_nedit_plus1
			break; // bailout
	}
	return min_nedit;
}

int _nedit_for_Proffset(const RoSeq *P, const RoSeq *S, int Proffset,
		int max_nedit, int loose_Proffset, int *min_width)
{
	int max_nedit_plus1, *prev_row, *curr_row, row_length,
	    B, b, i, iplus1, jmin, j, min_nedit;
	char Pc;

	min_nedit = 0;
	error("_nedit_for_Proffset() is not ready yet, sorry!");
	return min_nedit;
}


/****************************************************************************
 * --- .Call ENTRY POINT ---
 * Arguments:
 *   pattern
 *   subject
 *   at: positions of P's first or last letter in S;
 *   at_type: 0 if 'at' contains the positions of P's first letter in S,
 *            1 if 'at' contains the positions of P's last letter in S;
 *   max_mismatch:
 *   ans_type: 0 for a logical vector indicating whether there is a match or
 *             not at the specified positions, 1 for an integer vector giving
 *             the number of mismatches at the specified positions;
 */
SEXP match_pattern_at(SEXP pattern, SEXP subject, SEXP at, SEXP at_type,
		SEXP max_mismatch, SEXP with_indels, SEXP fixed, SEXP ans_type)
{
	RoSeq P, S;
	int at_length, at_type0, max_mm, indels, fixedP, fixedS, ans_type0,
	    i, *at_elt, *ans_elt, offset, nmismatch, min_width;
	SEXP ans;

	P = _get_XString_asRoSeq(pattern);
	S = _get_XString_asRoSeq(subject);
	at_length = LENGTH(at);
	at_type0 = INTEGER(at_type)[0];
	max_mm = INTEGER(max_mismatch)[0];
	indels = LOGICAL(with_indels)[0] && max_mm != 0;
	fixedP = LOGICAL(fixed)[0];
	fixedS = LOGICAL(fixed)[1];
	if (indels && !(fixedP && fixedS))
		error("when 'with.indels' is TRUE, only 'fixed=TRUE' is supported for now");
	ans_type0 = INTEGER(ans_type)[0];
	if (ans_type0) {
		PROTECT(ans = NEW_INTEGER(at_length));
		ans_elt = INTEGER(ans);
	} else {
		PROTECT(ans = NEW_LOGICAL(at_length));
		ans_elt = LOGICAL(ans);
	}
	if (!indels)
		_select_nmismatch_at_Pshift_fun(fixedP, fixedS);

	for (i = 0, at_elt = INTEGER(at); i < at_length; i++, at_elt++, ans_elt++)
	{
		if (*at_elt == NA_INTEGER) {
			*ans_elt = ans_type0 ? NA_INTEGER : NA_LOGICAL;
			continue;
		}
		if (indels) {
			offset = *at_elt - 1;
			if (at_type0 == 0)
				nmismatch = _nedit_for_Ploffset(&P, &S, offset, max_mm, 1, &min_width);
			else
				nmismatch = _nedit_for_Proffset(&P, &S, offset, max_mm, 1, &min_width);
		} else {
			if (at_type0 == 0)
				offset = *at_elt - 1;
			else
				offset = *at_elt - P.nelt;
			nmismatch = _selected_nmismatch_at_Pshift_fun(&P, &S, offset, max_mm);
		}
		*ans_elt = ans_type0 ? nmismatch : (nmismatch <= max_mm);
	}
	UNPROTECT(1);
	return ans;
}


/*
 * Arguments are the same as for match_pattern_at() except for:
 *   subject: XStringSet object.
 */
SEXP vmatch_pattern_at(SEXP pattern, SEXP subject, SEXP at, SEXP at_type,
		SEXP max_mismatch, SEXP with_indels, SEXP fixed, SEXP ans_type)
{
	RoSeq P, S_elt;
	CachedXStringSet S;
	int at_length, at_type0, max_mm, indels, fixedP, fixedS, ans_type0,
	    S_length, i, j, *at_elt, *ans_elt, offset, nmismatch, min_width;
	SEXP ans;

	P = _get_XString_asRoSeq(pattern);
	S = _new_CachedXStringSet(subject);
	at_length = LENGTH(at);
	at_type0 = INTEGER(at_type)[0];
	max_mm = INTEGER(max_mismatch)[0];
	indels = LOGICAL(with_indels)[0] && max_mm != 0;
	fixedP = LOGICAL(fixed)[0];
	fixedS = LOGICAL(fixed)[1];
	if (indels && !(fixedP && fixedS))
		error("when 'with.indels' is TRUE, only 'fixed=TRUE' is supported for now");
	ans_type0 = INTEGER(ans_type)[0];
	S_length = _get_XStringSet_length(subject);
	if (ans_type0) {
		PROTECT(ans = allocMatrix(INTSXP, S_length, at_length));
	} else {
		PROTECT(ans = allocMatrix(LGLSXP, S_length, at_length));
	}
	if (!indels)
		_select_nmismatch_at_Pshift_fun(fixedP, fixedS);

	for (i = 0; i < S_length; i++) {
		S_elt = _get_CachedXStringSet_elt_asRoSeq(&S, i);
		if (ans_type0) {
			ans_elt = INTEGER(ans) + i;
		} else {
			ans_elt = LOGICAL(ans) + i;
		}
		for (j = 0, at_elt = INTEGER(at); j < at_length; j++, at_elt++)
		{
			if (*at_elt == NA_INTEGER) {
				*ans_elt = ans_type0 ? NA_INTEGER : NA_LOGICAL;
			} else {
				if (indels) {
					offset = *at_elt - 1;
					if (at_type0 == 0)
						nmismatch = _nedit_for_Ploffset(&P, &S_elt, offset, max_mm, 1, &min_width);
					else
						nmismatch = _nedit_for_Proffset(&P, &S_elt, offset, max_mm, 1, &min_width);
				} else {
					if (at_type0 == 0)
						offset = *at_elt - 1;
					else
						offset = *at_elt - P.nelt;
					nmismatch = _selected_nmismatch_at_Pshift_fun(&P, &S_elt, offset, max_mm);
				}
				*ans_elt = ans_type0 ? nmismatch : (nmismatch <= max_mm);
			}
			ans_elt = ans_elt + S_length;
		}
	}
	UNPROTECT(1);
	return ans;
}
