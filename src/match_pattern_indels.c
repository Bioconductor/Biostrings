/****************************************************************************
               A SIMPLE MATCHING ALGO WITH SUPPORT FOR INDELS
	                     Author: Herve Pages
 ****************************************************************************/
#include "Biostrings.h"

static void test_match_pattern_indels(const char *p, const char *s,
		int max_nmis, const char *expected_matches)
{
	cachedCharSeq P, S;

	Rprintf("P=%s S=%s max_nmis=%d expected_matches=%s\n", p, s,
		max_nmis, expected_matches);
	P.seq = p;
	P.length = strlen(P.seq);
	S.seq = s;
	S.length = strlen(S.seq);
	_match_pattern_indels(&P, &S, max_nmis, 1, 1);
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
		_init_match_reporting("MATCHES_AS_NULL", 1);
		test_match_pattern_indels(p, s, 0, "30:34");
		test_match_pattern_indels(p, s, 1, "1:4, 14:18, 30:34");
		test_match_pattern_indels(p, s, 2, "1:4, 8:10, 14:18, 21:23, 30:34");
	}
#else
	Rprintf("Debug mode not available in file %s\n", __FILE__);
#endif
	return R_NilValue;
}

#ifdef DEBUG_BIOSTRINGS
static void print_match(int start, int width,
		const cachedCharSeq *P, const cachedCharSeq *S,
		const BytewiseOpTable *bytewise_match_table)
{
	int end, j0, nedit0, width0;
	char mbuf[1001];

	if (width >= sizeof(mbuf))
		error("sizeof(mbuf) too small");
	j0 = start - 1;
	end = j0 + width;
	snprintf(mbuf, width + 1, "%s", S->seq + j0);
	nedit0 = _nedit_for_Ploffset(P, S, j0, P->length, 1, &width0,
				     bytewise_match_table);
	Rprintf("start=%d end=%d (%s) nedit0=%d\n", start, end, mbuf, nedit0);
	return;
}
#endif

/*
 * In order to avoid pollution by redundant matches, we report only the "best
 * local matches". A substring S' of S is a best local match iff:
 *   (a) nedit(P, S') <= max_nmis
 *   (b) for every substring S1 of S', nedit(P, S1) > nedit(P, S')
 *   (c) for every substring S2 of S that contains S', nedit(P, S2) >= nedit(P, S')
 * One nice property of best local matches is that their first and last letters
 * are guaranteed to match letters in P.
 * The report_provisory_match() function will store a provisory match and will
 * hold it until it is replaced by a better one or until it's guaranteed to be 
 * a best local match (then it's reported as a match).
 */
static int provisory_match_start, provisory_match_end, provisory_match_width,
	   provisory_match_nedit;

static void report_provisory_match(int start, int width, int nedit)
{
	int end;

	end = start + width - 1;
	if (provisory_match_nedit != -1) {
		// Given how we walk on S, 'start' is always guaranteed to be >
		// 'provisory_match_start'.
		if (end > provisory_match_end)
			_report_match(provisory_match_start, provisory_match_width);
		else if (nedit > provisory_match_nedit)
			return;
	}
	provisory_match_start = start;
	provisory_match_end = end;
	provisory_match_width = width;
	provisory_match_nedit = nedit;
	return;
}

static ByteTrTable byte2offset;

void _match_pattern_indels(const cachedCharSeq *P, const cachedCharSeq *S,
		int max_nmis, int fixedP, int fixedS)
{
	int i0, j0, max_nmis1, nedit1, width1;
	char c0;
	const BytewiseOpTable *bytewise_match_table;
	cachedCharSeq P1;

	if (P->length <= 0)
		error("empty pattern");
	bytewise_match_table = _select_bytewise_match_table(fixedP, fixedS);
	if (!fixedP || !fixedS)
		error("'fixed' must be TRUE when 'algorithm=\"indels\"' (for now)");
	// Before we can support fixedP=FALSE or fixedS=FALSE in
	// _match_pattern_indels(), we need to support them in
	// _init_byte2offset_with_cachedCharSeq().
	_init_byte2offset_with_cachedCharSeq(&byte2offset, P, 0);
	provisory_match_nedit = -1; // means no provisory match yet
	j0 = 0;
	while (j0 < S->length) {
		while (1) {
			c0 = S->seq[j0];
			i0 = byte2offset.byte2code[(unsigned char) c0];
			if (i0 != NA_INTEGER) break;
			j0++;
			if (j0 >= S->length) goto done;
		}
		P1.seq = P->seq + i0 + 1;
		P1.length = P->length - i0 - 1;
		max_nmis1 = max_nmis - i0;
/*
#ifdef DEBUG_BIOSTRINGS
		if (debug) {
			Rprintf("[DEBUG] _match_pattern_indels(): "
				"j0=%d c0=%c i0=%d max_nmis1=%d\n", j0, c0, i0, max_nmis1);
		}
#endif
*/
		if (max_nmis1 >= 0) {
			if (max_nmis1 == 0) {
				nedit1 = _nmismatch_at_Pshift(&P1, S, j0 + 1,
							max_nmis1,
							bytewise_match_table);
				width1 = P1.length;
			} else {
				nedit1 = _nedit_for_Ploffset(&P1, S, j0 + 1,
							max_nmis1, 1, &width1,
							bytewise_match_table);
			}
			if (nedit1 <= max_nmis1) {
#ifdef DEBUG_BIOSTRINGS
				if (debug) {
					Rprintf("[DEBUG] _match_pattern_indels(): "
						"provisory match found at ");
					print_match(j0 + 1, width1 + 1, P, S,
						    bytewise_match_table);
				}
#endif
				report_provisory_match(j0 + 1, width1 + 1, nedit1 + i0);
			}
		}
		j0++;
	}
	done:
	if (provisory_match_nedit != -1)
		_report_match(provisory_match_start, provisory_match_width);
	return;
}

