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

static double compute_pwm_score(const double *pwm, int pwm_ncol, const char *S, int nS, int pwm_shift)
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

/*
 * --- .Call ENTRY POINT ---
 * PWM_score_starting_at() arguments are assumed to be:
 *   pwm: the Position Weight Matrix (numeric matrix with row names A, C, G and T)
 *   subject: a DNAString object containing the subject sequence
 *   base_codes: named integer vector of length 4 obtained with
 *       xscodes(subject, baseOnly=TRUE)
 *   starting_at: an integer vector of arbitrary length (NAs accepted)
 */
SEXP PWM_score_starting_at(SEXP pwm, SEXP subject, SEXP base_codes, SEXP starting_at)
{
	cachedCharSeq S;
	int pwm_ncol, i, *start_elt;
	SEXP ans;
	double *ans_elt;

	if (INTEGER(GET_DIM(pwm))[0] != 4)
		error("'pwm' must have 4 rows");
	pwm_ncol = INTEGER(GET_DIM(pwm))[1];
	S = cache_XRaw(subject);
	_init_byte2offset_with_INTEGER(byte2offset, base_codes, 1);
	PROTECT(ans = NEW_NUMERIC(LENGTH(starting_at)));
	for (i = 0, start_elt = INTEGER(starting_at), ans_elt = REAL(ans);
	     i < LENGTH(starting_at);
	     i++, start_elt++, ans_elt++) {
		if (*start_elt == NA_INTEGER) {
			*ans_elt = NA_REAL;
			continue;
		}
		*ans_elt = compute_pwm_score(REAL(pwm), pwm_ncol, S.seq, S.length, *start_elt - 1);
	}
	UNPROTECT(1);
	return ans;
}

/*
 * --- .Call ENTRY POINT ---
 * XString_match_PWM() arguments are assumed to be:
 *   pwm: the Position Weight Matrix (numeric matrix with row names A, C, G and T)
 *   subject: a DNAString object containing the subject sequence
 *   base_codes: named integer vector of length 4 obtained with
 *       xscodes(subject, baseOnly=TRUE)
 *   min_score: a single double (not NA)
 *   count_only: a single logical (not NA)
 */
SEXP XString_match_PWM(SEXP pwm, SEXP subject, SEXP base_codes, SEXP min_score, SEXP count_only)
{
	cachedCharSeq S;
	int pwm_ncol, is_count_only, n1, n2;
	double minscore;

	if (INTEGER(GET_DIM(pwm))[0] != 4)
		error("'pwm' must have 4 rows");
	pwm_ncol = INTEGER(GET_DIM(pwm))[1];
	S = cache_XRaw(subject);
	_init_byte2offset_with_INTEGER(byte2offset, base_codes, 1);
	minscore = REAL(min_score)[0];
	is_count_only = LOGICAL(count_only)[0];

	_init_match_reporting(is_count_only ?
		"MATCHES_AS_COUNTS" : "MATCHES_AS_RANGES");
	for (n1 = 0, n2 = pwm_ncol; n2 <= S.length; n1++, n2++) {
		if (compute_pwm_score(REAL(pwm), pwm_ncol, S.seq, S.length, n1) >= minscore)
			_report_match(n1 + 1, pwm_ncol);
	}
	return _reported_matches_asSEXP();
}

/*
 * --- .Call ENTRY POINT ---
 * XStringViews_match_PWM() arguments are assumed to be:
 *   pwm: the Position Weight Matrix (numeric matrix with row names A, C, G and T)
 *   subject: an XStringViews object containing the subject sequence
 *   views_start, views_width: 2 integer vectors describing views on 'subject'.
 *   base_codes: named integer vector of length 4 obtained with
 *       xscodes(subject, baseOnly=TRUE)
 *   min_score: a single double (not NA)
 *   count_only: a single logical (not NA)
 */
SEXP XStringViews_match_PWM(SEXP pwm,
		SEXP subject, SEXP views_start, SEXP views_width,
		SEXP base_codes, SEXP min_score, SEXP count_only)
{
	cachedCharSeq S, S_view;
	int pwm_ncol, is_count_only, n1, n2;
	int nviews, i, *view_start, *view_width, view_offset;
	double minscore;

	if (INTEGER(GET_DIM(pwm))[0] != 4)
		error("'pwm' must have 4 rows");
	pwm_ncol = INTEGER(GET_DIM(pwm))[1];
	S = cache_XRaw(subject);
	_init_byte2offset_with_INTEGER(byte2offset, base_codes, 1);
	minscore = REAL(min_score)[0];
	is_count_only = LOGICAL(count_only)[0];

	_init_match_reporting(is_count_only ?
		"MATCHES_AS_COUNTS" : "MATCHES_AS_RANGES");
	nviews = LENGTH(views_start);
	for (i = 0, view_start = INTEGER(views_start), view_width = INTEGER(views_width);
	     i < nviews;
	     i++, view_start++, view_width++)
	{
		view_offset = *view_start - 1;
		if (view_offset < 0 || view_offset + *view_width > S.length)
			error("'subject' has \"out of limits\" views");
		S_view.seq = S.seq + view_offset;
		S_view.length = *view_width;
		_shift_match_on_reporting(view_offset);
		for (n1 = 0, n2 = pwm_ncol; n2 <= S_view.length; n1++, n2++) {
			if (compute_pwm_score(REAL(pwm), pwm_ncol, S_view.seq, S_view.length, n1) >= minscore)
				_report_match(n1 + 1, pwm_ncol);
		}
	}
	return _reported_matches_asSEXP();
}

