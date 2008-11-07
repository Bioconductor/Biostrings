/****************************************************************************
                         A SIMPLE PWM MATCHING ALGO
                             Author: Herve Pages
 ****************************************************************************/
#include "Biostrings.h"

/*
 * Table used for fast look up between A, C, G, T internal codes and the
 * corresponding 0-based row indice (the row offset) in the PWM:
 *   A internal code     -> 0
 *   C internal code     -> 1
 *   G internal code     -> 2
 *   T internal code     -> 3
 *   other internal code -> NA_INTEGER
 */
static ByteTrTable DNAcode2PWMrowoffset;

static void init_DNAcode2PWMrowoffset()
{
	int i;

	for (i = 0; i < BYTETRTABLE_LENGTH; i++)
		DNAcode2PWMrowoffset[i] = NA_INTEGER;
	DNAcode2PWMrowoffset[(unsigned char) _DNAencode('A')] = 0;
	DNAcode2PWMrowoffset[(unsigned char) _DNAencode('C')] = 1;
	DNAcode2PWMrowoffset[(unsigned char) _DNAencode('G')] = 2;
	DNAcode2PWMrowoffset[(unsigned char) _DNAencode('T')] = 3;
	return;
}

static int compute_score(const int *pwm, int pwm_ncol, const char *S, int nS, int pwm_shift)
{
	int score, i, rowoffset;

	S += pwm_shift;
	nS -= pwm_shift;
	if (pwm_shift < 0 || nS < pwm_ncol)
		error("trying to compute the score from an invalid starting position");
	score = 0;
	for (i = 0; i < pwm_ncol; i++, pwm += 4, S++) {
		rowoffset = DNAcode2PWMrowoffset[(unsigned char) *S];
		if (rowoffset == NA_INTEGER)
			continue;
		score += pwm[rowoffset];
	}
	return score;
}

/*
 * --- .Call ENTRY POINT ---
 * PWM_score() arguments are assumed to be:
 *   pwm: the Position Weight Matrix (integer matrix with row names A, C, G and T)
 *   subject: a DNAString object containing the subject sequence
 *   start: an integer vector of arbitrary length (NAs accepted)
 */
SEXP PWM_score(SEXP pwm, SEXP subject, SEXP start)
{
	RoSeq S;
	int pwm_ncol, i, *start_elt, *ans_elt;
	SEXP ans;

	if (INTEGER(GET_DIM(pwm))[0] != 4)
		error("'pwm' must have 4 rows");
	pwm_ncol = INTEGER(GET_DIM(pwm))[1];
	S = _get_XString_asRoSeq(subject);
	init_DNAcode2PWMrowoffset();
	PROTECT(ans = NEW_INTEGER(LENGTH(start)));
	for (i = 0, start_elt = INTEGER(start), ans_elt = INTEGER(ans);
	     i < LENGTH(start);
	     i++, start_elt++, ans_elt++) {
		if (*start_elt == NA_INTEGER) {
			*ans_elt = NA_INTEGER;
			continue;
		}
		*ans_elt = compute_score(INTEGER(pwm), pwm_ncol, S.elts, S.nelt, *start_elt - 1);
	}
	UNPROTECT(1);
	return ans;
}

/*
 * --- .Call ENTRY POINT ---
 * match_PWM() arguments are assumed to be:
 *   pwm: the Position Weight Matrix (integer matrix with row names A, C, G and T)
 *   subject: a DNAString object containing the subject sequence
 *   min_score: an integer vector of length 1 (not NA)
 *   count_only: a logical vector of length 1 (not NA)
 */
SEXP match_PWM(SEXP pwm, SEXP subject, SEXP min_score, SEXP count_only)
{
	RoSeq S;
	int pwm_ncol, minscore, is_count_only, n1, n2;

	if (INTEGER(GET_DIM(pwm))[0] != 4)
		error("'pwm' must have 4 rows");
	pwm_ncol = INTEGER(GET_DIM(pwm))[1];
	S = _get_XString_asRoSeq(subject);
	minscore = INTEGER(min_score)[0];
	is_count_only = LOGICAL(count_only)[0];
	init_DNAcode2PWMrowoffset();	
	_init_match_reporting(is_count_only ? COUNT_MRMODE : START_MRMODE);
	for (n1 = 0, n2 = pwm_ncol; n2 <= S.nelt; n1++, n2++) {
		if (compute_score(INTEGER(pwm), pwm_ncol, S.elts, S.nelt, n1) >= minscore) {
			// The second arg (end) is ignored in match reporting
			// modes COUNT_MRMODE and START_MRMODE
			_report_match(n1 + 1, -1);
		}
	}
	// The SEXP returned by reported_matches_asSEXP() is UNPROTECTED
	// but you don't have to PROTECT it here since you are returning it
	// right away.
	return _reported_matches_asSEXP();
}

