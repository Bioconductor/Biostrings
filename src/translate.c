#include "Biostrings.h"
#include "IRanges_interface.h"
#include "XVector_interface.h"


#define TRANSLATE_ERROR	1
#define TRANSLATE_SOLVE	2
#define TRANSLATE_TO_X	3

static char errmsg_buf[200];

/*
 * Returns -1 if error, or the nb of trailing letters that were ignored
 * if successful (0, 1, or 2).
 */
static int fast_translate(const Chars_holder *dna, Chars_holder *aa,
			  char skip_code,
			  TwobitEncodingBuffer *teb, SEXP lkup)
{
	int phase, i, lkup_key;
	const char *c;
	char aa_letter;

	aa->length = phase = 0;
	_reset_twobit_signature(teb);
	for (i = 0, c = dna->seq; i < dna->length; i++, c++) {
		if (*c == skip_code)
			continue;
		lkup_key = _shift_twobit_signature(teb, *c);
		if (teb->lastin_twobit == NA_INTEGER) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "not a base at pos %d", i + 1);
			return -1;
		}
		if (phase < 2) {
			phase++;
			continue;
		}
		aa_letter = (char) INTEGER(lkup)[lkup_key];
		/* aa->seq is a const char * so we need to cast it to
		   char * before we can write to it */
		((char *) aa->seq)[aa->length++] = aa_letter;
		phase = 0;
	}
	return phase;
}

/*
 * Returns -1 if error, or the nb of trailing letters that were ignored
 * if successful (0, 1, or 2).
 */
static int translate(const Chars_holder *dna, Chars_holder *aa,
		     char skip_code,
		     int ncodes, ByteTrTable *byte2offset, SEXP lkup,
		     int if_non_ambig, int if_ambig)
{
	int phase, i, lkup_key, offset, is_fuzzy;
	const char *c;
	char aa_letter;

	aa->length = phase = 0;
	for (i = 0, c = dna->seq; i < dna->length; i++, c++) {
		if (*c == skip_code)
			continue;
		offset = byte2offset->byte2code[(unsigned char) *c];
		if (offset == NA_INTEGER) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "not a base at pos %d", i + 1);
			return -1;
		}
		if (phase == 0) {
			lkup_key = offset;
			is_fuzzy = 0;
			phase++;
			continue;
		}
		lkup_key *= ncodes;
		lkup_key += offset;
		if (offset >= 4)
			is_fuzzy = 1;
		if (phase < 2) {
			phase++;
			continue;
		}
		aa_letter = (char) INTEGER(lkup)[lkup_key];
		if (is_fuzzy) {
			/* codon is fuzzy */
			if (aa_letter != 'X') {
				/* non-ambiguous fuzzy codon */
				if (if_non_ambig == TRANSLATE_ERROR) {
					snprintf(errmsg_buf,
						 sizeof(errmsg_buf),
						 "non-ambiguous fuzzy codon "
						 "starting at pos %d", i - 1);
					return -1;
				}
				if (if_non_ambig == TRANSLATE_TO_X)
					aa_letter = 'X';
			} else {
				/* ambiguous fuzzy codon */
				if (if_ambig == TRANSLATE_ERROR) {
					snprintf(errmsg_buf,
						 sizeof(errmsg_buf),
						 "ambiguous fuzzy codon "
						 "starting at pos %d", i - 1);
					return -1;
				}
			}
		}
		/* aa->seq is a const char * so we need to cast it to
		   char * before we can write to it */
		((char *) aa->seq)[aa->length++] = aa_letter;
		phase = 0;
	}
	return phase;
}

/*
 * --- .Call ENTRY POINT ---
 * Return an AAStringSet object.
 */
SEXP DNAStringSet_translate(SEXP x, SEXP skip_code, SEXP dna_codes, SEXP lkup,
		SEXP if_non_ambig, SEXP if_ambig)
{
	char skip_code0;
	int ncodes, if_non_ambig0, if_ambig0, ans_length, i, errcode;
	TwobitEncodingBuffer teb;
	ByteTrTable byte2offset;
	const char *s1, *s2;
	XStringSet_holder X, Y;
	Chars_holder X_elt, Y_elt;
	SEXP ans, width, ans_width;

	skip_code0 = (unsigned char) INTEGER(skip_code)[0];
	ncodes = LENGTH(dna_codes);
	if (ncodes * ncodes * ncodes != LENGTH(lkup))
		error("Biostrings internal error in "
		      "DNAStringSet_translate(): length of 'lkup' "
		      "must equal length of 'dna_codes' power 3");
	if (ncodes == 4) {
		teb = _new_TwobitEncodingBuffer(dna_codes, 3, 0);
	} else {
		_init_byte2offset_with_INTEGER(&byte2offset, dna_codes, 1);
		s1 = CHAR(STRING_ELT(if_non_ambig, 0));
		if (strcmp(s1, "error") == 0)
			if_non_ambig0 = TRANSLATE_ERROR;
		else if (strcmp(s1, "solve") == 0)
			if_non_ambig0 = TRANSLATE_SOLVE;
		else if (strcmp(s1, "X") == 0)
			if_non_ambig0 = TRANSLATE_TO_X;
		else
			error("Biostrings internal error in "
			      "DNAStringSet_translate(): "
			      "invalid 'if_non_ambig' argument");
		s2 = CHAR(STRING_ELT(if_ambig, 0));
		if (strcmp(s2, "error") == 0)
			if_ambig0 = TRANSLATE_ERROR;
		else if (strcmp(s2, "X") == 0)
			if_ambig0 = TRANSLATE_TO_X;
		else
			error("Biostrings internal error in "
			      "DNAStringSet_translate(): "
			      "invalid 'if_ambig' argument");
	}
	X = _hold_XStringSet(x);
	ans_length = _get_length_from_XStringSet_holder(&X);
	PROTECT(width = NEW_INTEGER(ans_length));
	for (i = 0; i < ans_length; i++) {
		X_elt = _get_elt_from_XStringSet_holder(&X, i);
		INTEGER(width)[i] = X_elt.length / 3;
	}
	PROTECT(ans = alloc_XRawList("AAStringSet", "AAString", width));
	Y = _hold_XStringSet(ans);
	ans_width = _get_XStringSet_width(ans);
	for (i = 0; i < ans_length; i++) {
		X_elt = _get_elt_from_XStringSet_holder(&X, i);
		Y_elt = _get_elt_from_XStringSet_holder(&Y, i);
		errcode = ncodes == 4 ?
			fast_translate(&X_elt, &Y_elt,
				       skip_code0, &teb, lkup) :
			translate(&X_elt, &Y_elt,
				  skip_code0, ncodes, &byte2offset, lkup,
				  if_non_ambig0, if_ambig0);
		if (errcode == -1) {
			UNPROTECT(2);
			if (ans_length == 1)
				error("%s", errmsg_buf);
			else
				error("in 'x[[%d]]': %s", i + 1, errmsg_buf);
		}
		if (errcode >= 1) {
			if (errcode == 1) {
				snprintf(errmsg_buf, sizeof(errmsg_buf),
					 "last base was ignored");
			} else {
				snprintf(errmsg_buf, sizeof(errmsg_buf),
					 "last %d bases were ignored", errcode);
			}
			if (ans_length == 1)
				warning("%s", errmsg_buf);
			else
				warning("in 'x[[%d]]': %s", i + 1, errmsg_buf);
		}
		INTEGER(ans_width)[i] = Y_elt.length;
	}
	UNPROTECT(2);
	return ans;
}

