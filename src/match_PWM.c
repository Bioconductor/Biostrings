/****************************************************************************
                A SIMPLE POSITION WEIGHT MATRIX MATCHING ALGO
                             Author: Herve Pages
 ****************************************************************************/
#include "Biostrings.h"
#include "IRanges_interface.h"

/*
 * Table used for fast look up between A, C, G, T internal codes and the
 * corresponding 0-based row index (the row offset) in the PWM:
 *   A internal code     -> 0
 *   C internal code     -> 1
 *   G internal code     -> 2
 *   T internal code     -> 3
 *   other internal code -> NA_INTEGER
 */
static ByteTrTable byte2offset;

static double compute_pwm_score(const double *pwm, int pwm_ncol,
		const char *S, int nS, int pwm_shift)
{
	int i, rowoffset;
	double score;

	S += pwm_shift;
	nS -= pwm_shift;
	if (pwm_shift < 0 || nS < pwm_ncol)
		error("trying to compute the score from an invalid starting position");
	score = 0.00;
	for (i = 0; i < pwm_ncol; i++, pwm += 4, S++) {
		rowoffset = byte2offset[(unsigned char) *S];
		if (rowoffset == NA_INTEGER)
			continue;
		score += pwm[rowoffset];
	}
	return score;
}

static void _match_PWM_XString(const double *pwm, int pwm_ncol,
		const cachedCharSeq *S, int minscore)
{
	int n1, n2;
	double score;

	for (n1 = 0, n2 = pwm_ncol; n2 <= S->length; n1++, n2++) {
		score = compute_pwm_score(pwm, pwm_ncol, S->seq, S->length, n1);
		if (score >= minscore)
			_report_match(n1 + 1, pwm_ncol);
	}
	return;
}

/*
 * --- .Call ENTRY POINT ---
 * PWM_score_starting_at() arguments are assumed to be:
 *   pwm: matrix of doubles with row names A, C, G and T;
 *   subject: DNAString object containing the subject sequence;
 *   starting_at: an integer vector of arbitrary length (NAs accepted);
 *   base_codes: named integer vector of length 4 obtained with
 *       'xscodes(subject, baseOnly=TRUE)'.
 */
SEXP PWM_score_starting_at(SEXP pwm, SEXP subject, SEXP starting_at,
		SEXP base_codes)
{
	cachedCharSeq S;
	int pwm_ncol, ans_length, i, *start_elt;
	SEXP ans;
	double *ans_elt;

	if (INTEGER(GET_DIM(pwm))[0] != 4)
		error("'pwm' must have 4 rows");
	pwm_ncol = INTEGER(GET_DIM(pwm))[1];
	S = cache_XRaw(subject);
	_init_byte2offset_with_INTEGER(byte2offset, base_codes, 1);
	ans_length = LENGTH(starting_at);
	PROTECT(ans = NEW_NUMERIC(ans_length));
	for (i = 0, start_elt = INTEGER(starting_at), ans_elt = REAL(ans);
	     i < ans_length;
	     i++, start_elt++, ans_elt++)
	{
		if (*start_elt == NA_INTEGER) {
			*ans_elt = NA_REAL;
			continue;
		}
		*ans_elt = compute_pwm_score(REAL(pwm), pwm_ncol,
				S.seq, S.length, *start_elt - 1);
	}
	UNPROTECT(1);
	return ans;
}

/*
 * --- .Call ENTRY POINT ---
 * XString_match_PWM() arguments are assumed to be:
 *   pwm: matrix of doubles with row names A, C, G and T;
 *   subject: DNAString object containing the subject sequence;
 *   min_score: single double (not NA);
 *   count_only: single logical (not NA);
 *   base_codes: named integer vector of length 4 obtained with
 *       'xscodes(subject, baseOnly=TRUE)'.
 */
SEXP XString_match_PWM(SEXP pwm, SEXP subject,
		SEXP min_score, SEXP count_only, SEXP base_codes)
{
	cachedCharSeq S;
	int pwm_ncol, is_count_only;
	double minscore;

	if (INTEGER(GET_DIM(pwm))[0] != 4)
		error("'pwm' must have 4 rows");
	pwm_ncol = INTEGER(GET_DIM(pwm))[1];
	S = cache_XRaw(subject);
	minscore = REAL(min_score)[0];
	is_count_only = LOGICAL(count_only)[0];
	_init_byte2offset_with_INTEGER(byte2offset, base_codes, 1);
	_init_match_reporting(is_count_only ?
		"MATCHES_AS_COUNTS" : "MATCHES_AS_RANGES", 1);
	_match_PWM_XString(REAL(pwm), pwm_ncol, &S, minscore);
	return _reported_matches_asSEXP();
}

/*
 * --- .Call ENTRY POINT ---
 * XStringViews_match_PWM() arguments are assumed to be:
 *   pwm: matrix of doubles with row names A, C, G and T;
 *   subject: DNAString object containing the subject sequence;
 *   views_start, views_width: integer vectors describing views on 'subject';
 *   min_score: single double (not NA);
 *   count_only: single logical (not NA);
 *   base_codes: named integer vector of length 4 obtained with
 *       'xscodes(subject, baseOnly=TRUE)'.
 */
SEXP XStringViews_match_PWM(SEXP pwm,
		SEXP subject, SEXP views_start, SEXP views_width,
		SEXP min_score, SEXP count_only, SEXP base_codes)
{
	cachedCharSeq S, S_view;
	int pwm_ncol, is_count_only;
	int nviews, v, *start_p, *width_p, view_offset;
	double minscore;

	if (INTEGER(GET_DIM(pwm))[0] != 4)
		error("'pwm' must have 4 rows");
	pwm_ncol = INTEGER(GET_DIM(pwm))[1];
	S = cache_XRaw(subject);
	minscore = REAL(min_score)[0];
	is_count_only = LOGICAL(count_only)[0];
	_init_byte2offset_with_INTEGER(byte2offset, base_codes, 1);

	_init_match_reporting(is_count_only ?
		"MATCHES_AS_COUNTS" : "MATCHES_AS_RANGES", 1);
	nviews = LENGTH(views_start);
	for (v = 0,
	     start_p = INTEGER(views_start),
	     width_p = INTEGER(views_width);
	     v < nviews;
	     v++, start_p++, width_p++)
	{
		view_offset = *start_p - 1;
		if (view_offset < 0 || view_offset + *width_p > S.length)
			error("'subject' has \"out of limits\" views");
		S_view.seq = S.seq + view_offset;
		S_view.length = *width_p;
		_set_match_shift(view_offset);
		_match_PWM_XString(REAL(pwm), pwm_ncol, &S_view, minscore);
	}
	return _reported_matches_asSEXP();
}

