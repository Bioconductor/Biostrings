/****************************************************************************
 *                           The Twobit algorithm                           *
 *                   for constant width DNA dictionaries                    *
 *                                                                          *
 *                           Author: Herve Pages                            *
 *                                                                          *
 * Note: a constant width dictionary is a non-empty set of non-empty        *
 * words of the same length.                                                *
 ****************************************************************************/
#include "Biostrings.h"
#include "XVector_interface.h"
#include "IRanges_interface.h"

static int debug = 0;


SEXP debug_match_pdict_Twobit()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'match_pdict_Twobit.c'\n",
		debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'match_pdict_Twobit.c'\n");
#endif
	return R_NilValue;
}



/****************************************************************************
 *                                                                          *
 *                             A. PREPROCESSING                             *
 *                                                                          *
 ****************************************************************************/

static void init_twobit_sign2pos(SEXP twobit_sign2pos, int val0)
{
	int i;

	for (i = 0; i < LENGTH(twobit_sign2pos); i++)
		INTEGER(twobit_sign2pos)[i] = val0;
	return;
}

static int pp_pattern(SEXP twobit_sign2pos, TwobitEncodingBuffer *teb,
		const Chars_holder *pattern, int poffset)
{
	int twobit_sign, *pos0;

	//printf("poffset=%d: ", poffset);
	twobit_sign = _get_twobit_signature(teb, pattern);
	if (twobit_sign == NA_INTEGER)
		return -1;
	//printf("twobit_sign=%d\n", twobit_sign);
	pos0 = INTEGER(twobit_sign2pos) + twobit_sign;
	if (*pos0 == NA_INTEGER)
		*pos0 = poffset + 1;
	else
		_report_ppdup(poffset, *pos0);
	return 0;
}


/****************************************************************************
 * Turning our local data structures into an R list (SEXP)
 * -------------------------------------------------------
 */

/*
 * Twobit_asLIST() returns an R list with the following elements:
 *   - sign2pos: XInteger object;
 *   - high2low: an integer vector containing the mapping between duplicated and
 *         primary reads.
 */

static SEXP Twobit_asLIST(SEXP twobit_sign2pos)
{
	SEXP ans, ans_names, ans_elt;

	PROTECT(ans = NEW_LIST(2));

	/* set the names */
	PROTECT(ans_names = NEW_CHARACTER(2));
	SET_STRING_ELT(ans_names, 0, mkChar("sign2pos"));
	SET_STRING_ELT(ans_names, 1, mkChar("high2low"));
	SET_NAMES(ans, ans_names);
	UNPROTECT(1);

	/* set the "sign2pos" element */
	PROTECT(ans_elt = new_XInteger_from_tag("XInteger", twobit_sign2pos));
	SET_ELEMENT(ans, 0, ans_elt);
	UNPROTECT(1);

	/* set the "high2low" element */
	PROTECT(ans_elt = _get_ppdups_buf_asINTEGER());
	SET_ELEMENT(ans, 1, ans_elt);
	UNPROTECT(1);

	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * .Call entry point for preprocessing
 * -----------------------------------
 *
 * Arguments:
 *   tb:         the Trusted Band extracted from the original dictionary as a
 *               DNAStringSet object;
 *   pp_exclude: NULL or an integer vector of the same length as 'tb' where
 *               non-NA values indicate the elements to exclude from
 *               preprocessing;
 *   base_codes: the internal codes for A, C, G and T.
 *
 * See Twobit_asLIST() for a description of the returned SEXP.
 */

SEXP build_Twobit(SEXP tb, SEXP pp_exclude, SEXP base_codes)
{
	int tb_length, tb_width, poffset, twobit_len;
	XStringSet_holder tb_holder;
	Chars_holder pattern;
	TwobitEncodingBuffer teb;
	SEXP ans, twobit_sign2pos;

	tb_length = _get_XStringSet_length(tb);
	_init_ppdups_buf(tb_length);
	tb_width = -1;
	tb_holder = _hold_XStringSet(tb);
	for (poffset = 0; poffset < tb_length; poffset++) {
		/* Skip duplicated patterns */
		if (pp_exclude != R_NilValue
		 && INTEGER(pp_exclude)[poffset] != NA_INTEGER)
			continue;
		pattern = _get_elt_from_XStringSet_holder(&tb_holder, poffset);
		if (pattern.length == 0)
			error("empty trusted region for pattern %d",
			      poffset + 1);
		if (tb_width == -1) {
			tb_width = pattern.length;
			if (tb_width > 14)
				error("the width of the Trusted Band must "
				      "be <= 14 when 'type=\"Twobit\"'");
			teb = _new_TwobitEncodingBuffer(base_codes, tb_width, 0);
			twobit_len = 1 << (tb_width * 2); // 4^tb_width
			PROTECT(twobit_sign2pos = NEW_INTEGER(twobit_len));
			init_twobit_sign2pos(twobit_sign2pos, NA_INTEGER);
		} else if (pattern.length != tb_width) {
			error("all the trusted regions must have "
			      "the same length");
		}
		if (pp_pattern(twobit_sign2pos, &teb, &pattern, poffset) != 0) {
			UNPROTECT(1);
			error("non-base DNA letter found in Trusted Band "
			      "for pattern %d", poffset + 1);
		}
			
	}
	PROTECT(ans = Twobit_asLIST(twobit_sign2pos));
	UNPROTECT(2);
	return ans;
}



/****************************************************************************
 *                                                                          *
 *                             B. MATCH FINDING                             *
 *                                                                          *
 ****************************************************************************/

void walk_subject(const int *twobit_sign2pos, TwobitEncodingBuffer *teb,
		const Chars_holder *S, TBMatchBuf *tb_matches)
{
	int n, twobit_sign, P_id;
	const char *s;

	_reset_twobit_signature(teb);
	for (n = 1, s = S->ptr; n <= S->length; n++, s++) {
		twobit_sign = _shift_twobit_signature(teb, *s);
		if (twobit_sign == NA_INTEGER)
			continue;
		P_id = twobit_sign2pos[twobit_sign];
		if (P_id == NA_INTEGER)
			continue;
		_TBMatchBuf_report_match(tb_matches, P_id - 1, n);
	}
	return;
}

void _match_Twobit(SEXP pptb, const Chars_holder *S, int fixedS,
		TBMatchBuf *tb_matches)
{
	int tb_width;
	const int *twobit_sign2pos;
	SEXP base_codes;
	TwobitEncodingBuffer teb;

#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] ENTERING _match_Twobit()\n");
#endif
	tb_width = _get_PreprocessedTB_width(pptb);
	twobit_sign2pos = INTEGER(_get_Twobit_sign2pos_tag(pptb));
	base_codes = _get_PreprocessedTB_base_codes(pptb);
	teb = _new_TwobitEncodingBuffer(base_codes, tb_width, 0);
	if (!fixedS)
		error("cannot treat IUPAC extended letters in the subject "
		      "as ambiguities when 'pdict' is a PDict object of "
		      "the \"Twobit\" type");
	walk_subject(twobit_sign2pos, &teb, S, tb_matches);
#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] LEAVING _match_Twobit()\n");
#endif
	return;
}

