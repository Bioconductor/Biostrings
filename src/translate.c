#include "Biostrings.h"
#include "IRanges_interface.h"


static char errmsg_buf[200];

/* Returns 0 if OK, -1 if error, and 1 if warning. */
static int translate(const cachedCharSeq *dna, cachedCharSeq *aa,
		TwobitEncodingBuffer *teb, SEXP lkup, SEXP skipcode)
{
	int phase, i, lkup_key;
	const char *c;

	phase = 0;
	aa->length = 0;
	_reset_twobit_signature(teb);
	for (i = 0, c = dna->seq; i < dna->length; i++, c++) {
		lkup_key = _shift_twobit_signature(teb, *c);
		if (phase < 2) {
			phase++;
			continue;
		}
		if (lkup_key == NA_INTEGER) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "invalid codon starting at pos %d", i + 1);
			return -1;
		}
		/* aa->seq is a const char * so we need to cast it to
		   char * before we can write to it */
		((char *) aa->seq)[aa->length++] = (char) INTEGER(lkup)[lkup_key];
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
 * Return an AAString object.
 */
SEXP DNAString_translate(SEXP x, SEXP base_codes, SEXP lkup, SEXP skipcode)
{
	cachedCharSeq X, Y;
	int ans_length;
	SEXP ans;
	TwobitEncodingBuffer teb;
	int errcode;

	X = cache_XRaw(x);
	ans_length = X.length / 3;
	PROTECT(ans = alloc_XRaw("AAString", ans_length));
	Y = cache_XRaw(ans);
	teb = _new_TwobitEncodingBuffer(base_codes, 3, 0);
	errcode = translate(&X, &Y, &teb, lkup, skipcode);
	if (errcode == -1) {
		UNPROTECT(1);
		error("%s", errmsg_buf);
	}
	if (errcode == 1)
		warning("%s", errmsg_buf);
	UNPROTECT(1);
	return ans;
}

/*
 * --- .Call ENTRY POINT ---
 * Return an AAStringSet object.
 */
SEXP DNAStringSet_translate(SEXP x, SEXP base_codes, SEXP lkup, SEXP skipcode)
{
	cachedXStringSet X, Y;
	cachedCharSeq X_elt, Y_elt;
	int ans_length, i;
	SEXP ans, ans_width;
	TwobitEncodingBuffer teb;
	int errcode;

	X = _cache_XStringSet(x);
	ans_length = _get_cachedXStringSet_length(&X);
	PROTECT(ans_width = NEW_INTEGER(ans_length));
	for (i = 0; i < ans_length; i++) {
		X_elt = _get_cachedXStringSet_elt(&X, i);
		INTEGER(ans_width)[i] = X_elt.length / 3;
	}
	PROTECT(ans = alloc_XRawList("AAStringSet", "AAString", ans_width));
	Y = _cache_XStringSet(ans);
	teb = _new_TwobitEncodingBuffer(base_codes, 3, 0);
	for (i = 0; i < ans_length; i++) {
		X_elt = _get_cachedXStringSet_elt(&X, i);
		Y_elt = _get_cachedXStringSet_elt(&Y, i);
		errcode = translate(&X_elt, &Y_elt, &teb, lkup, skipcode);
		if (errcode == -1) {
			UNPROTECT(2);
			error("in 'x' element %d: %s", i + 1, errmsg_buf);
		}
		if (errcode == 1)
			warning("in 'x' element %d: %s", i + 1, errmsg_buf);
	}
	UNPROTECT(2);
	return ans;
}

