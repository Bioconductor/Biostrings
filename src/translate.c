#include "Biostrings.h"
#include "IRanges_interface.h"
#include "XVector_interface.h"


static char errmsg_buf[200];

/* Returns 0 if OK, -1 if error, and 1 if warning. */
static int translate(const Chars_holder *dna, Chars_holder *aa,
		TwobitEncodingBuffer *teb, SEXP lkup, char skipcode0)
{
	int phase, i, lkup_key;
	const char *c;

	phase = 0;
	aa->length = 0;
	_reset_twobit_signature(teb);
	for (i = 0, c = dna->seq; i < dna->length; i++, c++) {
		if (*c == skipcode0)
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
		/* aa->seq is a const char * so we need to cast it to
		   char * before we can write to it */
		((char *) aa->seq)[aa->length++] =
			(char) INTEGER(lkup)[lkup_key];
		phase = 0;
	}
	if (phase == 1) {
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "last base was ignored");
		return 1;
	}
	if (phase == 2) {
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "last %d bases were ignored", phase);
		return 1;
	}
	return 0;
}

/*
 * --- .Call ENTRY POINT ---
 * Return an AAStringSet object.
 */
SEXP DNAStringSet_translate(SEXP x, SEXP base_codes, SEXP lkup, SEXP skipcode)
{
	XStringSet_holder X, Y;
	Chars_holder X_elt, Y_elt;
	char skipcode0;
	int ans_length, i, errcode;
	SEXP ans, width, ans_width;
	TwobitEncodingBuffer teb;

	X = _hold_XStringSet(x);
	skipcode0 = (unsigned char) INTEGER(skipcode)[0];
	ans_length = _get_length_from_XStringSet_holder(&X);
	PROTECT(width = NEW_INTEGER(ans_length));
	for (i = 0; i < ans_length; i++) {
		X_elt = _get_elt_from_XStringSet_holder(&X, i);
		INTEGER(width)[i] = X_elt.length / 3;
	}
	PROTECT(ans = alloc_XRawList("AAStringSet", "AAString", width));
	Y = _hold_XStringSet(ans);
	teb = _new_TwobitEncodingBuffer(base_codes, 3, 0);
	ans_width = _get_XStringSet_width(ans);
	for (i = 0; i < ans_length; i++) {
		X_elt = _get_elt_from_XStringSet_holder(&X, i);
		Y_elt = _get_elt_from_XStringSet_holder(&Y, i);
		errcode = translate(&X_elt, &Y_elt, &teb, lkup, skipcode0);
		if (errcode == -1) {
			UNPROTECT(2);
			if (ans_length == 1)
				error("%s", errmsg_buf);
			else
				error("in 'x[[%d]]': %s", i + 1, errmsg_buf);
		}
		if (errcode == 1) {
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

