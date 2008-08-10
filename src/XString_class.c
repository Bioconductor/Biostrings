/****************************************************************************
 *                  Basic manipulation of XString objects                   *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"

static int debug = 0;

SEXP debug_XString_class()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'XString_class.c'\n", debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'XString_class.c'\n");
#endif
	return R_NilValue;
}


/****************************************************************************
 * Encoding/decoding XString data.
 */

static int DNA_enc_chrtrtable[CHRTRTABLE_LENGTH],
	   DNA_dec_chrtrtable[CHRTRTABLE_LENGTH],
	   RNA_enc_chrtrtable[CHRTRTABLE_LENGTH],
	   RNA_dec_chrtrtable[CHRTRTABLE_LENGTH];

const int *get_enc_chrtrtable(const char *class)
{
	if (strcmp(class, "DNAString") == 0)
		return DNA_enc_chrtrtable;
	else if (strcmp(class, "RNAString") == 0)
		return RNA_enc_chrtrtable;
	return NULL;
}

const int *get_dec_chrtrtable(const char *class)
{
	if (strcmp(class, "DNAString") == 0)
		return DNA_dec_chrtrtable;
	else if (strcmp(class, "RNAString") == 0)
		return RNA_dec_chrtrtable;
	return NULL;
}

static void copy_lkup(const int *lkup1, int len1, int *lkup2, int len2)
{
	int i;

	if (len1 > len2)
		error("Biostrings internal error in copy_lkup(): len1 > len2");
	for (i = 0; i < len1; i++)
		lkup2[i] = lkup1[i];
	for ( ; i < len2; i++)
		lkup2[i] = NA_INTEGER;
	return;
}

/*
 * --- .Call ENTRY POINT ---
 */
SEXP init_DNAlkups(SEXP enc_lkup, SEXP dec_lkup)
{
	copy_lkup(INTEGER(enc_lkup), LENGTH(enc_lkup),
		  DNA_enc_chrtrtable, CHRTRTABLE_LENGTH);
	copy_lkup(INTEGER(dec_lkup), LENGTH(dec_lkup),
		  DNA_dec_chrtrtable, CHRTRTABLE_LENGTH);
	return R_NilValue;
}

char _DNAencode(char c)
{
	int code;

	code = DNA_enc_chrtrtable[(unsigned char) c];
	if (code == NA_INTEGER)
		error("_DNAencode(): key %d not in lookup table", (int) c);
	return code;
}

char _DNAdecode(char code)
{
	int c;

	c = DNA_dec_chrtrtable[(unsigned char) code];
	if (c == NA_INTEGER)
		error("_DNAdecode(): key %d not in lookup table", (int) code);
	return c;
}

/*
 * --- .Call ENTRY POINT ---
 */
SEXP init_RNAlkups(SEXP enc_lkup, SEXP dec_lkup)
{
	copy_lkup(INTEGER(enc_lkup), LENGTH(enc_lkup),
		  RNA_enc_chrtrtable, CHRTRTABLE_LENGTH);
	copy_lkup(INTEGER(dec_lkup), LENGTH(dec_lkup),
		  RNA_dec_chrtrtable, CHRTRTABLE_LENGTH);
	return R_NilValue;
}

char _RNAencode(char c)
{
	int code;

	code = RNA_enc_chrtrtable[(unsigned char) c];
	if (code == NA_INTEGER)
		error("_RNAencode(): key %d not in lookup table", (int) c);
	return code;
}

char _RNAdecode(char code)
{
	int c;

	c = RNA_dec_chrtrtable[(unsigned char) code];
	if (c == NA_INTEGER)
		error("_RNAdecode(): key %d not in lookup table", (int) code);
	return (char) c;
}


/****************************************************************************
 * Low-level manipulation of XString objects.
 */

SEXP _get_XString_xdata(SEXP x)
{
	return GET_SLOT(x, install("xdata"));
}

RoSeq _get_XString_asRoSeq(SEXP x)
{
	RoSeq seq;
	SEXP tag;
	int offset;

	tag = _get_XRaw_tag(_get_XString_xdata(x));
	offset = INTEGER(GET_SLOT(x, install("offset")))[0];
	seq.elts = (const char *) (RAW(tag) + offset);
	seq.nelt = INTEGER(GET_SLOT(x, install("length")))[0];
	return seq;
}

/*
 * Never try to make this a .Call() entry point!
 * Its arguments are NOT duplicated so it would be a disaster if they were
 * coming from the user space.
 */
SEXP _new_XString(const char *class, SEXP xdata, int offset, int length)
{
	SEXP class_def, ans;

	class_def = MAKE_CLASS(class);
	PROTECT(ans = NEW_OBJECT(class_def));
	SET_SLOT(ans, mkChar("xdata"), xdata);
	SET_SLOT(ans, mkChar("offset"), ScalarInteger(offset));
	SET_SLOT(ans, mkChar("length"), ScalarInteger(length));
	UNPROTECT(1);
	return ans;
}

SEXP _new_XString_from_RoSeqs(const char *class, const RoSeqs *seqs)
{
	const int *enc_lkup;
        SEXP lkup, xdata, ans;

	enc_lkup = get_enc_chrtrtable(class);
	if (enc_lkup == NULL) {
		lkup = R_NilValue;
	} else {
		PROTECT(lkup = NEW_INTEGER(CHRTRTABLE_LENGTH));
		copy_lkup(enc_lkup, CHRTRTABLE_LENGTH,
			  INTEGER(lkup), LENGTH(lkup));
	}
	PROTECT(xdata = _new_XRaw_from_RoSeqs(seqs, lkup));
	PROTECT(ans = _new_XString(class, xdata, 0, _get_XRaw_length(xdata)));
	if (enc_lkup == NULL)
		UNPROTECT(2);
	else
		UNPROTECT(3);
	return ans;
}


/****************************************************************************
 * Utilities for creating an XString instance in 2 steps: first create the
 * skeleton (with junk data in it), then fill it with some character data.
 */

/*
 * Allocate only. The 'xdata' slot is not initialized (it contains junk,
 * or zeros).
 */
SEXP _alloc_XString(const char *class, int length)
{
	SEXP tag, xdata, ans;

	PROTECT(tag = NEW_RAW(length));
	PROTECT(xdata = _new_XRaw(tag));
	PROTECT(ans = _new_XString(class, xdata, 0, length));
	UNPROTECT(3);
	return ans;
}

void _write_RoSeq_to_XString(SEXP x, int start, const RoSeq *seq, int encode)
{
	int offset;
	const int *enc_chrtrtable;

	offset = INTEGER(GET_SLOT(x, install("offset")))[0];
	enc_chrtrtable = encode ? get_enc_chrtrtable(_get_class(x)) : NULL;
	_write_RoSeq_to_XRaw(_get_XString_xdata(x), offset + start - 1, seq, enc_chrtrtable);
	return;
}

