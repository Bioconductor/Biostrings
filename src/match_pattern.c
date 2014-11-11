/****************************************************************************
 *                    EXACT AND INEXACT PATTERN MATCHING                    *
 *		             Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"
#include "XVector_interface.h"
#include "IRanges_interface.h"

static int debug = 0;

SEXP debug_match_pattern()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'match_pattern.c'\n",
                 debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'match_pattern.c'\n");
#endif
	return R_NilValue;
}


/****************************************************************************
 * A memcmp-based implementation of the "naive" method for exact matching.
 *
 * Here is how the "naive" (aka "memcmp", aka "blunt") method for finding
 * exact matches is described in Dan Gusfield book "Algorithms on strings,
 * trees, and sequences" (slightly modified):
 *   The naive method aligns the left end of P (the pattern) with the left
 *   end of S (the subject) and then compares the characters of P and S left
 *   to right until either two unequal characters are found or until P is
 *   exhausted, in which case an occurrence of P is reported. In either case,
 *   P is then shifted one place to the right, and the comparisons are
 *   restarted from the left end of P. This process repeats until the right
 *   end of P shifts past the right end of S.
 *
 * Why implement this inefficient "naive" method?
 * - For QC: we can validate other more sophisticated matching algo by
 *   comparing their results to those obtains with the "naive" method.
 * - To use as a reference when comparing performance.
 */

static void match_naive_exact(const Chars_holder *P, const Chars_holder *S)
{
	const char *p, *s;
	int plen, slen, start, n2;

	if (P->length <= 0)
		error("empty pattern");
	p = P->seq;
	plen = P->length;
	s = S->seq;
	slen = S->length;
	for (start = 1, n2 = plen; n2 <= slen; start++, n2++, s++) {
		if (memcmp(p, s, plen) == 0)
			_report_match(start, P->length);
	}
	return;
}


/****************************************************************************
 * An implementation of the "naive" method for inexact matching.
 */

static void match_naive_inexact(const Chars_holder *P, const Chars_holder *S,
		int max_nmis, int min_nmis, int fixedP, int fixedS)
{
	int Pshift, // position of pattern left-most char relative to the subject
	    n2, // 1 + position of pattern right-most char relative to the subject
	    min_Pshift, max_n2, nmis;
	const BytewiseOpTable *bytewise_match_table;

	if (P->length <= 0)
		error("empty pattern");
	bytewise_match_table = _select_bytewise_match_table(fixedP, fixedS);
	min_Pshift = P->length <= max_nmis ? 1 - P->length : -max_nmis;
	max_n2 = S->length - min_Pshift;
	for (Pshift = min_Pshift, n2 = min_Pshift + P->length;
	     n2 <= max_n2;
	     Pshift++, n2++)
	{
		nmis = _nmismatch_at_Pshift(P, S, Pshift, max_nmis,
					    bytewise_match_table);
		if (nmis <= max_nmis && nmis >= min_nmis)
			_report_match(Pshift + 1, P->length);
	}
	return;
}


/****************************************************************************
 * _match_pattern_XString() and _match_pattern_XStringViews()
 */

void _match_pattern_XString(const Chars_holder *P, const Chars_holder *S,
		SEXP max_mismatch, SEXP min_mismatch,
		SEXP with_indels, SEXP fixed,
		const char *algo)
{
	int max_nmis, min_nmis, fixedP, fixedS;

	max_nmis = INTEGER(max_mismatch)[0];
	min_nmis = INTEGER(min_mismatch)[0];
	if (max_nmis < P->length - S->length
	 || min_nmis > P->length)
		return;
	fixedP = LOGICAL(fixed)[0];
	fixedS = LOGICAL(fixed)[1];
	if (P->length <= max_nmis || strcmp(algo, "naive-inexact") == 0)
		match_naive_inexact(P, S, max_nmis, min_nmis, fixedP, fixedS);
	else if (strcmp(algo, "naive-exact") == 0)
		match_naive_exact(P, S);
	else if (strcmp(algo, "boyer-moore") == 0)
		_match_pattern_boyermoore(P, S, -1, 0);
	else if (strcmp(algo, "shift-or") == 0)
		_match_pattern_shiftor(P, S, max_nmis, fixedP, fixedS);
	else if (strcmp(algo, "indels") == 0)
		_match_pattern_indels(P, S, max_nmis, fixedP, fixedS);
	else
		error("\"%s\": unknown algorithm", algo);
	return;
}

void _match_pattern_XStringViews(const Chars_holder *P,
		const Chars_holder *S, SEXP views_start, SEXP views_width,
		SEXP max_mismatch, SEXP min_mismatch,
		SEXP with_indels, SEXP fixed,
		const char *algo)
{
	Chars_holder S_view;
	int nviews, v, *view_start, *view_width, view_offset;

	nviews = LENGTH(views_start);
	for (v = 0,
	     view_start = INTEGER(views_start),
	     view_width = INTEGER(views_width);
	     v < nviews;
	     v++, view_start++, view_width++)
	{
		view_offset = *view_start - 1;
		if (view_offset < 0 || view_offset + *view_width > S->length)
			error("'subject' has \"out of limits\" views");
		S_view.seq = S->seq + view_offset;
		S_view.length = *view_width;
		_set_match_shift(view_offset);
		_match_pattern_XString(P, &S_view,
			max_mismatch, min_mismatch, with_indels, fixed,
			algo);
	}
	return;
}


/****************************************************************************
 * --- .Call ENTRY POINTS ---
 *
 * Arguments:
 *   pattern, subject: XString objects;
 *   max_mismatch: (single integer) the max number of mismatching letters;
 *   min_mismatch: (single integer) the min number of mismatching letters;
 *   with_indels: single logical;
 *   fixed: logical vector of length 2;
 *   algorithm: single string;
 *   count_only: single logical.
 *
 * If with_indels is FALSE: all matches have the length of the pattern.
 * Otherwise, matches are of variable length (>= length(pattern) - max_mismatch
 * and <= length(pattern) + max_mismatch).
 */

/* --- .Call ENTRY POINT --- */
SEXP XString_match_pattern(SEXP pattern, SEXP subject,
		SEXP max_mismatch, SEXP min_mismatch,
		SEXP with_indels, SEXP fixed,
		SEXP algorithm, SEXP count_only)
{
	Chars_holder P, S;
	const char *algo;
	int is_count_only;

	P = hold_XRaw(pattern);
	S = hold_XRaw(subject);
	algo = CHAR(STRING_ELT(algorithm, 0));
	is_count_only = LOGICAL(count_only)[0];
	_init_match_reporting(is_count_only ?
		"MATCHES_AS_COUNTS" : "MATCHES_AS_RANGES", 1);
	_match_pattern_XString(&P, &S,
		max_mismatch, min_mismatch, with_indels, fixed,
		algo);
	return _reported_matches_asSEXP();
}

/* --- .Call ENTRY POINT ---
 * Arguments are the same as for XString_match_pattern() except for:
 *   subject: XString object;
 *   views_start, views_width: 2 integer vectors describing views on 'subject'.
 */
SEXP XStringViews_match_pattern(SEXP pattern,
		SEXP subject, SEXP views_start, SEXP views_width,
		SEXP max_mismatch, SEXP min_mismatch,
		SEXP with_indels, SEXP fixed,
		SEXP algorithm, SEXP count_only)
{
	Chars_holder P, S;
	const char *algo;
	int is_count_only;

	P = hold_XRaw(pattern);
	S = hold_XRaw(subject);
	algo = CHAR(STRING_ELT(algorithm, 0));
	is_count_only = LOGICAL(count_only)[0];
	_init_match_reporting(is_count_only ?
		"MATCHES_AS_COUNTS" : "MATCHES_AS_RANGES", 1);
	_match_pattern_XStringViews(&P,
		&S, views_start, views_width,
		max_mismatch, min_mismatch, with_indels, fixed,
		algo);
	return _reported_matches_asSEXP();
}

/* --- .Call ENTRY POINT ---
 * Arguments are the same as for XString_match_pattern() except for:
 *   subject: XStringSet object.
 */
SEXP XStringSet_vmatch_pattern(SEXP pattern, SEXP subject,
		SEXP max_mismatch, SEXP min_mismatch,
		SEXP with_indels, SEXP fixed,
		SEXP algorithm, SEXP ms_mode)
{
	Chars_holder P, S_elt;
	XStringSet_holder S;
	int S_length, j;
	const char *algo;

	P = hold_XRaw(pattern);
	S = _hold_XStringSet(subject);
	S_length = _get_XStringSet_length(subject);
	algo = CHAR(STRING_ELT(algorithm, 0));
	_init_match_reporting(CHAR(STRING_ELT(ms_mode, 0)), S_length);
	for (j = 0; j < S_length; j++) {
		S_elt = _get_elt_from_XStringSet_holder(&S, j);
		_set_active_PSpair(j);
		_match_pattern_XString(&P, &S_elt,
			max_mismatch, min_mismatch, with_indels, fixed,
			algo);
	}
	return _MatchBuf_as_SEXP(_get_internal_match_buf(), R_NilValue);
}

