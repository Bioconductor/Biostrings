/****************************************************************************
 *                      Low-level matching functions                        *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"
#include "IRanges_interface.h"

static int debug = 0;

SEXP debug_lowlevel_matching()
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
 * nmismatch_at_Pshift_*()
 *
 * The 4 static functions below stop counting mismatches if the number
 * exceeds 'max_nmis'. The caller can disable this by passing 'P->length' to
 * the 'max_nmis' arg.
 *
 * fixedP | fixedS | letters *p and *s match iff...
 * --------------------------------------------------------
 * TRUE   | TRUE   | ...they are equal
 * TRUE   | FALSE  | ...bits at 1 in *p are also at 1 in *s
 * FALSE  | TRUE   | ...bits at 1 in *s are also at 1 in *p
 * FALSE  | FALSE  | ...they share at least one bit at 1
 */

static int nmismatch_at_Pshift_fixedPfixedS(const cachedCharSeq *P,
		const cachedCharSeq *S, int Pshift, int max_nmis)
{
	int nmis, i, j;
	const char *p, *s;

	nmis = 0;
	for (i = 0, j = Pshift, p = P->seq, s = S->seq + Pshift;
	     i < P->length;
	     i++, j++, p++, s++)
	{
		if (j >= 0 && j < S->length && *p == *s)
			continue;
		if (nmis++ >= max_nmis)
			break;
	}
	return nmis;
}

static int nmismatch_at_Pshift_fixedPnonfixedS(const cachedCharSeq *P,
		const cachedCharSeq *S, int Pshift, int max_nmis)
{
	int nmis, i, j;
	const char *p, *s;

	nmis = 0;
	for (i = 0, j = Pshift, p = P->seq, s = S->seq + Pshift;
	     i < P->length;
	     i++, j++, p++, s++)
	{
		if (j >= 0 && j < S->length && ((*p) & ~(*s)) == 0)
			continue;
		if (nmis++ >= max_nmis)
			break;
	}
	return nmis;
}

static int nmismatch_at_Pshift_nonfixedPfixedS(const cachedCharSeq *P,
		const cachedCharSeq *S, int Pshift, int max_nmis)
{
	int nmis, i, j;
	const char *p, *s;

	nmis = 0;
	for (i = 0, j = Pshift, p = P->seq, s = S->seq + Pshift;
	     i < P->length;
	     i++, j++, p++, s++)
	{
		if (j >= 0 && j < S->length && (~(*p) & (*s)) == 0)
			continue;
		if (nmis++ >= max_nmis)
			break;
	}
	return nmis;
}

static int nmismatch_at_Pshift_nonfixedPnonfixedS(const cachedCharSeq *P,
		const cachedCharSeq *S, int Pshift, int max_nmis)
{
	int nmis, i, j;
	const char *p, *s;

	nmis = 0;
	for (i = 0, j = Pshift, p = P->seq, s = S->seq + Pshift;
	     i < P->length;
	     i++, j++, p++, s++)
	{
		if (j >= 0 && j < S->length && ((*p) & (*s)))
			continue;
		if (nmis++ >= max_nmis)
			break;
	}
	return nmis;
}

int (*_selected_nmismatch_at_Pshift_fun)(const cachedCharSeq *P,
		const cachedCharSeq *S, int Pshift, int max_nmis);

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
#define MAX_NEDIT 100
#define MAX_ROW_LENGTH (2*MAX_NEDIT+1)

static int row1_buf[MAX_ROW_LENGTH], row2_buf[MAX_ROW_LENGTH];

#define SWAP_NEDIT_BUFS(prev_row, curr_row) \
{ \
	int *tmp; \
	tmp = (prev_row); \
	(prev_row) = (curr_row); \
	(curr_row) = tmp; \
}

#define PROPAGATE_NEDIT(curr_row, B, prev_row, S, Si, Pc, row_length) \
{ \
	int nedit, B2, nedit2; \
	nedit = (prev_row)[(B)] + ((Si) < 0 || (Si) >= (S)->length || (S)->seq[(Si)] != (Pc)); \
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
int _nedit_for_Ploffset(const cachedCharSeq *P, const cachedCharSeq *S,
		int Ploffset, int max_nedit, int loose_Ploffset, int *min_width)
{
	int max_nedit_plus1, *prev_row, *curr_row, row_length,
	    a, B, b, min_Si, min_nedit,
	    Pi, Si; // 0-based letter pos in P and S, respectively
	char Pc;

#ifdef DEBUG_BIOSTRINGS
	if (debug) Rprintf("[DEBUG] _nedit_for_Ploffset():\n");
#endif
	if (P->length == 0)
		return 0;
	if (max_nedit == 0)
		error("Biostrings internal error in _nedit_for_Ploffset(): ",
		      "use _selected_nmismatch_at_Pshift_fun() when 'max_nedit' is 0");
	max_nedit_plus1 = max_nedit + 1;
	if (max_nedit > P->length)
		max_nedit = P->length;
	// from now max_nedit <= P->length
	if (max_nedit > MAX_NEDIT)
		error("'max.nedit' too big");
	prev_row = row1_buf;
	curr_row = row2_buf;
	row_length = 2 * max_nedit + 1;
	min_Si = Ploffset;

	// STAGE 0:
	for (B = max_nedit, b = 0; B < row_length; B++, b++)
		curr_row[B] = b;
#ifdef DEBUG_BIOSTRINGS
	if (debug) print_curr_row("STAGE0", curr_row, max_nedit, row_length);
#endif

	// STAGE 1 (1st for() loop): no attempt is made to bailout during
	// this stage because the smallest value in curr_row is guaranteed
	// to be <= a < max_nedit.
	for (a = 1, Pi = 0; a < max_nedit; a++, Pi++) {
		Pc = P->seq[Pi]; // Pi < a < max_nedit <= P->length
		SWAP_NEDIT_BUFS(prev_row, curr_row);
		B = max_nedit - a;
		curr_row[B++] = a;
		for (Si = min_Si; B < row_length; B++, Si++)
			PROPAGATE_NEDIT(curr_row, B, prev_row, S, Si, Pc, row_length);
#ifdef DEBUG_BIOSTRINGS
		if (debug) print_curr_row("STAGE1", curr_row, max_nedit - a, row_length);
#endif
	}

	// STAGE 2: no attempt is made to bailout during this stage either.
	Pc = P->seq[Pi];
	SWAP_NEDIT_BUFS(prev_row, curr_row);
	B = 0;
	curr_row[B++] = min_nedit = a;
	*min_width = 0;
	for (Si = min_Si; B < row_length; B++, Si++) {
		PROPAGATE_NEDIT(curr_row, B, prev_row, S, Si, Pc, row_length);
		if (curr_row[B] < min_nedit) {
			min_nedit = curr_row[B];
			*min_width = Si - Ploffset + 1;
		}
	}
#ifdef DEBUG_BIOSTRINGS
	if (debug) print_curr_row("STAGE2", curr_row, 0, row_length);
#endif
	a++;
	Pi++;

	// STAGE 3 (2nd for() loop): with attempt to bailout.
	for ( ; Pi < P->length; Pi++, a++, min_Si++) {
		Pc = P->seq[Pi];
		SWAP_NEDIT_BUFS(prev_row, curr_row);
		min_nedit = a;
		*min_width = 0;
		for (B = 0, Si = min_Si; B < row_length; B++, Si++) {
			PROPAGATE_NEDIT(curr_row, B, prev_row, S, Si, Pc, row_length);
			if (curr_row[B] < min_nedit) {
				min_nedit = curr_row[B];
				*min_width = Si - Ploffset + 1;
			}
		}
#ifdef DEBUG_BIOSTRINGS
		if (debug) print_curr_row("STAGE3", curr_row, 0, row_length);
#endif
		// 'min_nedit > max_nedit_plus1' should actually never happen.
		if (min_nedit >= max_nedit_plus1)
			break; // bailout
	}
	return min_nedit;
}

int _nedit_for_Proffset(const cachedCharSeq *P, const cachedCharSeq *S,
		int Proffset, int max_nedit, int loose_Proffset, int *min_width)
{
	int max_nedit_plus1, *prev_row, *curr_row, row_length,
	    a, B, b, max_Si, min_nedit,
	    Pi, Si; // 0-based letter pos in P and S, respectively
	char Pc;

#ifdef DEBUG_BIOSTRINGS
	if (debug) Rprintf("[DEBUG] _nedit_for_Proffset():\n");
#endif
	if (P->length == 0)
		return 0;
	if (max_nedit == 0)
		error("Biostrings internal error in _nedit_for_Proffset(): ",
		      "use _selected_nmismatch_at_Pshift_fun() when 'max_nedit' is 0");
	max_nedit_plus1 = max_nedit + 1;
	if (max_nedit > P->length)
		max_nedit = P->length;
	// from now max_nedit <= P->length
	if (max_nedit > MAX_NEDIT)
		error("'max.nedit' too big");
	prev_row = row1_buf;
	curr_row = row2_buf;
	row_length = 2 * max_nedit + 1;
	max_Si = Proffset;
	min_nedit = 0;

	// STAGE 0:
	for (B = max_nedit, b = 0; B < row_length; B++, b++)
		curr_row[B] = b;
#ifdef DEBUG_BIOSTRINGS
	if (debug) print_curr_row("STAGE0", curr_row, max_nedit, row_length);
#endif

	// STAGE 1 (1st for() loop): no attempt is made to bailout during
	// this stage because the smallest value in curr_row is guaranteed
	// to be <= a < max_nedit.
	for (a = 1, Pi = P->length - 1; a < max_nedit; a++, Pi--) {
		Pc = P->seq[Pi];
		SWAP_NEDIT_BUFS(prev_row, curr_row);
		B = max_nedit - a;
		curr_row[B++] = a;
		for (Si = max_Si; B < row_length; B++, Si--)
			PROPAGATE_NEDIT(curr_row, B, prev_row, S, Si, Pc, row_length);
#ifdef DEBUG_BIOSTRINGS
		if (debug) print_curr_row("STAGE1", curr_row, max_nedit - a, row_length);
#endif
	}

	// STAGE 2: no attempt is made to bailout during this stage either.
	Pc = P->seq[Pi];
	SWAP_NEDIT_BUFS(prev_row, curr_row);
	B = 0;
	curr_row[B++] = min_nedit = a;
	*min_width = 0;
	for (Si = max_Si; B < row_length; B++, Si--) {
		PROPAGATE_NEDIT(curr_row, B, prev_row, S, Si, Pc, row_length);
		if (curr_row[B] < min_nedit) {
			min_nedit = curr_row[B];
			*min_width = Proffset - Si + 1;
		}
	}
#ifdef DEBUG_BIOSTRINGS
	if (debug) print_curr_row("STAGE2", curr_row, 0, row_length);
#endif
	a++;
	Pi--;

	// STAGE 3 (2nd for() loop): with attempt to bailout.
	for ( ; Pi >= 0; Pi--, a++, max_Si--) {
		Pc = P->seq[Pi];
		SWAP_NEDIT_BUFS(prev_row, curr_row);
		min_nedit = a;
		*min_width = 0;
		for (B = 0, Si = max_Si; B < row_length; B++, Si--) {
			PROPAGATE_NEDIT(curr_row, B, prev_row, S, Si, Pc, row_length);
			if (curr_row[B] < min_nedit) {
				min_nedit = curr_row[B];
				*min_width = Proffset - Si + 1;
			}
		}
#ifdef DEBUG_BIOSTRINGS
		if (debug) print_curr_row("STAGE3", curr_row, 0, row_length);
#endif
		// 'min_nedit > max_nedit_plus1' should actually never happen.
		if (min_nedit >= max_nedit_plus1)
			break; // bailout
	}
	return min_nedit;
}


/****************************************************************************
 * nedit_at()
 */

static int nedit_at(const cachedCharSeq *P, const cachedCharSeq *S,
		int at, int at_type0, int max_nmis, int with_indels0,
		int fixedP, int fixedS)
{
	int offset, nmis, min_width;

	if (!with_indels0 || max_nmis == 0) {
		if (at_type0 == 0)
			offset = at - 1;
		else
			offset = at - P->length;
		_select_nmismatch_at_Pshift_fun(fixedP, fixedS);
		return _selected_nmismatch_at_Pshift_fun(P, S, offset, max_nmis);
	}
	if (!(fixedP && fixedS))
		error("when 'with.indels' is TRUE, only 'fixed=TRUE' is supported for now");
	offset = at - 1;
	if (at_type0 == 0)
		nmis = _nedit_for_Ploffset(P, S, offset, max_nmis, 1, &min_width);
	else
		nmis = _nedit_for_Proffset(P, S, offset, max_nmis, 1, &min_width);
	return nmis;
}


/****************************************************************************
 * match_pattern_at()
 */

static void check_mismatch_lengths(int at_length,
		SEXP max_mismatch, SEXP min_mismatch, int ans_type0)
{
	if ((at_length == 0 && LENGTH(max_mismatch) > 1)
         || (at_length != 0 && LENGTH(max_mismatch) > at_length))
		warning("'max_mismatch' is longer than 'at' "
			"(remaining elements are ignored)");
	if ((at_length == 0 && LENGTH(min_mismatch) > 1)
         || (at_length != 0 && LENGTH(min_mismatch) > at_length))
		warning("'min_mismatch' is longer than 'at' "
			"(remaining elements are ignored)");
	if (at_length == 0)
		return;
	if (LENGTH(max_mismatch) == 0)
		error("'max_mismatch' must have at least 1 element");
	if (ans_type0 == 0)
		return;
	if (LENGTH(min_mismatch) == 0)
		error("'min_mismatch' must have at least 1 element");
	return;
}

static void match_pattern_at(const cachedCharSeq *P, const cachedCharSeq *S,
		SEXP at, int at_type0,
		SEXP max_mismatch, SEXP min_mismatch, int with_indels0,
		int fixedP, int fixedS, int ans_type0, int *ans_elt,
		int auto_reduce_pattern0)
{
	int at_length, i, k1, k2, *at_elt, max_nmis, min_nmis, nmis, is_matching;
	cachedCharSeq my_p = *P;

	at_length = LENGTH(at);
	if (ans_type0 >= 2)
		*ans_elt = NA_INTEGER;
	for (i = 1, k1 = k2 = 0, at_elt = INTEGER(at);
	     i <= at_length;
	     i++, k1++, k2++, at_elt++)
	{
		if (k1 >= LENGTH(max_mismatch))
			k1 = 0; /* recycle */
		if (k2 >= LENGTH(min_mismatch))
			k2 = 0; /* recycle */
		if (*at_elt == NA_INTEGER) {
			switch (ans_type0) {
				case 0: *(ans_elt++) = NA_INTEGER; break;
				case 1: *(ans_elt++) = NA_LOGICAL; break;
			}
			continue;
		}
		max_nmis = INTEGER(max_mismatch)[k1];
		if (max_nmis == NA_INTEGER)
			max_nmis = my_p.length;
		nmis = nedit_at(&my_p, S, *at_elt, at_type0,
				max_nmis, with_indels0, fixedP, fixedS);
		if (auto_reduce_pattern0 && i < at_length) {
			if (at_type0 == 0)
				my_p.seq++;
			my_p.length--;
		}
		if (ans_type0 == 0) {
			*(ans_elt++) = nmis;
			continue;
		}
		min_nmis = INTEGER(min_mismatch)[k2];
		if (min_nmis == NA_INTEGER)
			min_nmis = 0;
		is_matching = nmis <= max_nmis && nmis >= min_nmis;
		if (ans_type0 == 1) {
			*(ans_elt++) = is_matching;
			continue;
		}
		if (is_matching) {
			*ans_elt = ans_type0 == 2 ? i : *at_elt;
			break;
		}
	}
	return;
}



/****************************************************************************
 *                        --- .Call ENTRY POINTS ---                        *
 ****************************************************************************/

/*
 * XString_match_pattern_at() arguments:
 *   pattern: XString object of same base type as 'subject'.
 *   subject: XString object.
 *   at: the 1-based positions of 'pattern' with respect to 'subject'.
 *   at_type: how to interpret the positions in 'at'. If 0, then they are
 *            those of the first letter of 'pattern'. If 1, then they are
 *            those of its last letter.
 *   max_mismatch, min_mismatch: integer vectors of length >= 1 recycled to
 *            the length of 'at'. If the number of effective mismatches is
 *            <= 'max_mismatch', then it is reported accurately. Otherwise
 *            any number > 'max_mismatch' could be reported. This is to allow
 *            the matching functions used as backends to which
 *            XString_match_pattern_at() delegates to implement early bailout
 *            strategies.
 *   with_indels: TRUE or FALSE. If TRUE, then the "number of mismatches" at
 *            a given position means the smallest edit distance between the
 *            'pattern' and all the substrings in 'subject' that start (if
 *            'at_type' is 0) or end (if 'at_type' is 1) at this position.
 *   fixed: a logical of length 2.
 *   ans_type: a single integer specifying the type of answer to return:
 *       0: ans is an integer vector of the same length as 'at';
 *       1: ans is a logical vector of the same length as 'at';
 *       2: ans is the lowest *index* (1-based position) in 'at' for which
 *          a match occurred (or NA if no match occurred);
 *       3: ans is the first *value* in 'at' for which a match occurred
 *          (or NA if no match occurred).
 *   auto_reduce_pattern: a logical scalar passed to match_pattern_at
 *	     instructing it whether to reduce the pattern by 1 letter
 *	     for each successive (at, max_mismatch) "pair", from its
 *	     beginning with at_type 0 and from its end with at_type 1.
 *	     Currently, the at vector is always constant of length
 *	     nchar(pattern), and furthermore equal to 1 with at_type 0,
 *	     i.e. "starting.at=1", and equal to the subject length with
 *	     at_type 1, i.e. "ending.at=nchar(subject)".
 */
SEXP XString_match_pattern_at(SEXP pattern, SEXP subject, SEXP at, SEXP at_type,
		SEXP max_mismatch, SEXP min_mismatch, SEXP with_indels, SEXP fixed,
		SEXP ans_type, SEXP auto_reduce_pattern)
{
	cachedCharSeq P, S;
	int at_length, at_type0, with_indels0, fixedP, fixedS, ans_type0, *ans_elt;
	int auto_reduce_pattern0 = LOGICAL(auto_reduce_pattern)[0];
	SEXP ans;

	P = cache_XRaw(pattern);
	S = cache_XRaw(subject);
	at_length = LENGTH(at);
	at_type0 = INTEGER(at_type)[0];
	with_indels0 = LOGICAL(with_indels)[0];
	fixedP = LOGICAL(fixed)[0];
	fixedS = LOGICAL(fixed)[1];
	ans_type0 = INTEGER(ans_type)[0];
	check_mismatch_lengths(at_length,
			max_mismatch, min_mismatch, ans_type0);
	switch (ans_type0) {
		case 0:
			PROTECT(ans = NEW_INTEGER(at_length));
			ans_elt = INTEGER(ans);
			break;
		case 1:
			PROTECT(ans = NEW_LOGICAL(at_length));
			ans_elt = LOGICAL(ans);
			break;
		case 2: case 3:
			PROTECT(ans = NEW_INTEGER(1));
			ans_elt = INTEGER(ans);
			break;
		default: error("invalid 'ans_type' value (%d)", ans_type0);
	}
	match_pattern_at(&P, &S, at, at_type0,
			 max_mismatch, min_mismatch, with_indels0,
			 fixedP, fixedS, ans_type0, ans_elt,
			 auto_reduce_pattern0);
	UNPROTECT(1);
	return ans;
}


/*
 * XStringSet_vmatch_pattern_at() arguments are the same as for
 * XString_match_pattern_at() except for:
 *   subject: XStringSet object.
 */
SEXP XStringSet_vmatch_pattern_at(SEXP pattern, SEXP subject, SEXP at, SEXP at_type,
		SEXP max_mismatch, SEXP min_mismatch, SEXP with_indels, SEXP fixed,
		SEXP ans_type, SEXP auto_reduce_pattern)
{
	cachedCharSeq P, S_elt;
	cachedXStringSet S;
	int S_length, at_length, at_type0, with_indels0, fixedP, fixedS,
	    ans_type0, *ans_elt, ans_nrow, i;
	int auto_reduce_pattern0 = LOGICAL(auto_reduce_pattern)[0];
	SEXP ans;

	P = cache_XRaw(pattern);
	S = _cache_XStringSet(subject);
	S_length = _get_cachedXStringSet_length(&S);
	at_length = LENGTH(at);
	at_type0 = INTEGER(at_type)[0];
	with_indels0 = LOGICAL(with_indels)[0];
	fixedP = LOGICAL(fixed)[0];
	fixedS = LOGICAL(fixed)[1];
	ans_type0 = INTEGER(ans_type)[0];
	check_mismatch_lengths(at_length,
			max_mismatch, min_mismatch, ans_type0);
	switch (ans_type0) {
		case 0:
			PROTECT(ans = allocMatrix(INTSXP, at_length, S_length));
			ans_elt = INTEGER(ans);
			ans_nrow = at_length;
			break;
		case 1:
			PROTECT(ans = allocMatrix(LGLSXP, at_length, S_length));
			ans_elt = LOGICAL(ans);
			ans_nrow = at_length;
			break;
		case 2: case 3:
			PROTECT(ans = NEW_INTEGER(S_length));
			ans_elt = INTEGER(ans);
			ans_nrow = 1;
			break;
		default: error("invalid 'ans_type' value (%d)", ans_type0);
	}
	for (i = 0; i < S_length; i++, ans_elt += ans_nrow) {
		S_elt = _get_cachedXStringSet_elt(&S, i);
		match_pattern_at(&P, &S_elt, at, at_type0,
				 max_mismatch, min_mismatch, with_indels0,
				 fixedP, fixedS, ans_type0, ans_elt,
				 auto_reduce_pattern0);
	}
	UNPROTECT(1);
	return ans;
}


/*
 * XStringSet_dist_hamming() used by stringDist, method = "hamming".
 */
SEXP XStringSet_dist_hamming(SEXP x)
{
	cachedCharSeq x_i, x_j;
	cachedXStringSet X;
	int X_length, *ans_elt, i, j;
	unsigned long ans_length;
	SEXP ans;

	X = _cache_XStringSet(x);
	X_length = _get_cachedXStringSet_length(&X);
	if (X_length < 2)
		return NEW_INTEGER(0);

	x_i = _get_cachedXStringSet_elt(&X, 0);
	for (j = 1; j < X_length; j++) {
		x_j = _get_cachedXStringSet_elt(&X, j);
		if (x_i.length != x_j.length)
		      error("Hamming distance requires equal length strings");
	}

	ans_length = ((unsigned long) X_length) *
		     ((unsigned long) X_length - 1) / 2;
	if (ans_length > INT_MAX)
		error("result would be too big an object");
	PROTECT(ans = NEW_INTEGER((int) ans_length));
	ans_elt = INTEGER(ans);

	for (i = 0; i < (X_length - 1); i++) {
		x_i = _get_cachedXStringSet_elt(&X, i);
		for (j = (i+1); j < X_length; j++, ans_elt++) {
			x_j = _get_cachedXStringSet_elt(&X, j);
			*ans_elt = nedit_at(&x_i, &x_j, 1, 0,
					    0, 0, 1, 1);
		}
	}
	UNPROTECT(1);
	return ans;
}

