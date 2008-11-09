/****************************************************************************
               A SIMPLE MATCHING ALGO WITH SUPPORT FOR INDELS
	                     Author: Herve Pages
 ****************************************************************************/
#include "Biostrings.h"

static void test_match_pattern_indels(const char *p, const char *s, int max_mm, const char *expected_matches)
{
	RoSeq P, S;

	Rprintf("P=%s S=%s max_mm=%d expected_matches=%s\n", p, s, max_mm, expected_matches);
	P.elts = p;
	P.nelt = strlen(P.elts);
	S.elts = s;
	S.nelt = strlen(S.elts);
	_match_pattern_indels(&P, &S, max_mm, 1, 1);
	return;
}

static int debug = 0;

SEXP debug_match_pattern_indels()
{
#ifdef DEBUG_BIOSTRINGS
	const char *p = "ABCDE", *s = "BCDExAxBCDDxDABCxExxABDCExExAABCDEE";
	debug = !debug;
	Rprintf("Debug mode turned %s in file %s\n",
		debug ? "on" : "off", __FILE__);
	if (debug == 1) {
		_init_match_reporting(2);
		test_match_pattern_indels(p, s, 0, "30:34");
		test_match_pattern_indels(p, s, 1, "1:4, 14:18, 30:34");
		test_match_pattern_indels(p, s, 2, "1:4, 6:10, 14:18, 21:25, 30:34");
	}
#else
	Rprintf("Debug mode not available in file %s\n", __FILE__);
#endif
	return R_NilValue;
}

static ByteTrTable byte2offset;

void _match_pattern_indels(const RoSeq *P, const RoSeq *S,
		int max_mm, int fixedP, int fixedS)
{
	int i0, j0, max_mm1, nedit1, width1;
	char c0;
	RoSeq P1;

	if (P->nelt <= 0)
		error("empty pattern");
	_select_nmismatch_at_Pshift_fun(fixedP, fixedS);
	if (!fixedP || !fixedS)
		error("'fixed' must be TRUE when 'algorithm=\"indels\"' (for now)");
	// Before we can support fixedP=FALSE or fixedS=FALSE in
	// _match_pattern_indels(), we need to support them in
	// _init_byte2offset_with_RoSeq() and _nedit_for_Ploffset().
	_init_byte2offset_with_RoSeq(byte2offset, P, 0);
	j0 = 0;
	while (j0 < S->nelt) {
		while (1) {
			c0 = S->elts[j0];
			i0 = byte2offset[(unsigned char) c0];
			if (i0 != NA_INTEGER) break;
			j0++;
			if (j0 >= S->nelt) return;
		}
		P1.elts = P->elts + i0 + 1;
		P1.nelt = P->nelt - i0 - 1;
		max_mm1 = max_mm - i0;
/*
#ifdef DEBUG_BIOSTRINGS
		if (debug) {
			Rprintf("[DEBUG] _match_pattern_indels(): "
				"j0=%d c0=%c i0=%d max_mm1=%d\n", j0, c0, i0, max_mm1);
		}
#endif
*/
		if (max_mm1 >= 0) {
			if (max_mm1 == 0) {
				nedit1 = _selected_nmismatch_at_Pshift_fun(&P1, S, j0 + 1, max_mm1);
				width1 = P1.nelt;
			} else {
				nedit1 = _nedit_for_Ploffset(&P1, S, j0 + 1, max_mm1, 1, &width1);
			}
			if (nedit1 <= max_mm1) {
				_report_match(j0 + 1, 0);
#ifdef DEBUG_BIOSTRINGS
				if (debug) {
					int start, end, width, nedit0, width0;
					char mbuf[1001];

					width = width1 + 1;
					if (width >= sizeof(mbuf))
						error("sizeof(mbuf) too small");
					start = j0 + 1; // because j0 is 0-based
					end = start + width - 1;
					snprintf(mbuf, width + 1, "%s", S->elts + j0);
					nedit0 = _nedit_for_Ploffset(P, S, j0, P->nelt, 1, &width0);
					Rprintf("[DEBUG] _match_pattern_indels(): match at "
						"start=%d end=%d (%s) nedit0=%d\n", start, end, mbuf, nedit0);
				}
#endif
			}
		}
		j0++;
	}
	return;
}

