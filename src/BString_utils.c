/*
 * Low-level manipulation of BString and BStringList objects.
 */
#include "Biostrings.h"

static int DNA_enc_lkup[256], DNA_dec_lkup[256];
static int RNA_enc_lkup[256], RNA_dec_lkup[256];

static void copy_lkup(const int *lkup1, int len1, int *lkup2, int len2)
{
	int i;

	if (len1 > len2)
		error("Biostrings internal error in copy_lkup(): len1 > len2");
	for (i = 0; i < len1; i++)
		lkup2[i] = lkup1[i];
	for ( ; i < len2; i++)
		lkup2[i] = NA_INTEGER;
}

static int class2code(const char *class)
{
	error("class2code() not ready yet");
	return 0;
}

static const char *code2class(int code)
{
	error("code2class() not ready yet");
	return 0;
}

const char *getBString_class(SEXP x)
{
	return CHAR(STRING_ELT(GET_CLASS(x), 0));
}

static SEXP getBString_data(SEXP x)
{
	return GET_SLOT(x, install("data"));
}

SEXP init_DNAlkups(SEXP enc_lkup, SEXP dec_lkup)
{
	copy_lkup(INTEGER(enc_lkup), LENGTH(enc_lkup),
		DNA_enc_lkup, sizeof(DNA_enc_lkup) / sizeof(int));
	copy_lkup(INTEGER(dec_lkup), LENGTH(dec_lkup),
		DNA_dec_lkup, sizeof(DNA_dec_lkup) / sizeof(int));
	return R_NilValue;
}

char _DNAencode(char c)
{
	int code;

	code = DNA_enc_lkup[(unsigned char) c];
	if (code == NA_INTEGER)
		error("_DNAencode: key %d not in lookup table", (int) c);
	return code;
}

char _DNAdecode(char code)
{
	int c;

	c = DNA_dec_lkup[(unsigned char) code];
	if (c == NA_INTEGER)
		error("_DNAdecode: key %d not in lookup table", (int) code);
	return c;
}

SEXP init_RNAlkups(SEXP enc_lkup, SEXP dec_lkup)
{
	copy_lkup(INTEGER(enc_lkup), LENGTH(enc_lkup),
		RNA_enc_lkup, sizeof(RNA_enc_lkup) / sizeof(int));
	copy_lkup(INTEGER(dec_lkup), LENGTH(dec_lkup),
		RNA_dec_lkup, sizeof(RNA_dec_lkup) / sizeof(int));
	return R_NilValue;
}

char _RNAencode(char c)
{
	int code;

	code = RNA_enc_lkup[(unsigned char) c];
	if (code == NA_INTEGER)
		error("_RNAencode: key %d not in lookup table", (int) c);
	return code;
}

char _RNAdecode(char code)
{
	int c;

	c = RNA_dec_lkup[(unsigned char) code];
	if (c == NA_INTEGER)
		error("_RNAdecode: key %d not in lookup table", (int) code);
	return (char) c;
}

const char *_get_BString_charseq(SEXP x, int *length)
{
	SEXP xp;
	int offset;

	xp = GET_SLOT(getBString_data(x), install("xp"));
	offset = INTEGER(GET_SLOT(x, install("offset")))[0];
	*length = INTEGER(GET_SLOT(x, install("length")))[0];
	return (const char *) (RAW(R_ExternalPtrTag(xp)) + offset);
}

/* UNTESTED */
SEXP mkBString(const char *class, SEXP data, int offset, int length)
{
	SEXP class_def, ans;

	class_def = MAKE_CLASS(class);
	PROTECT(ans = NEW_OBJECT(class_def));
	SET_SLOT(ans, mkChar("data"), data);
	SET_SLOT(ans, mkChar("offset"), ScalarInteger(offset));
	SET_SLOT(ans, mkChar("length"), ScalarInteger(length));
	UNPROTECT(1);
	return ans;
}

/* NOT a Call() entry point! */
SEXP CHARSXP_to_BString(SEXP seq, SEXP start, SEXP nchar, SEXP lkup, SEXP proto)
{
	SEXP ans, data, data_xp;
	int length;

	PROTECT(ans = duplicate(proto));
	PROTECT(data = CHARSXP_to_XRaw(seq, start, nchar, lkup));
	data_xp = GET_SLOT(data, install("xp"));
	length = LENGTH(R_ExternalPtrTag(data_xp));
	SET_SLOT(ans, mkChar("data"), data);
	SET_SLOT(ans, mkChar("offset"), ScalarInteger(0));
	SET_SLOT(ans, mkChar("length"), ScalarInteger(length));
	UNPROTECT(2);
	return ans;
}

SEXP charseq_to_BString(SEXP seq, SEXP start, SEXP nchar, SEXP lkup, SEXP proto)
{
	return CHARSXP_to_BString(STRING_ELT(seq, 0), start, nchar, lkup, proto);
}

/* Return the list, NOT the BStringList object! */
SEXP charseqs_to_BStrings(SEXP seqs, SEXP start, SEXP nchar, SEXP lkup, SEXP proto)
{
	SEXP ans, seq;
	int nseq, i;

	nseq = LENGTH(seqs);
	PROTECT(ans = NEW_LIST(nseq));
	for (i = 0; i < nseq; i++) {
		seq = STRING_ELT(seqs, i);
		SET_ELEMENT(ans, i, CHARSXP_to_BString(seq, start, nchar, lkup, proto));
	}
	UNPROTECT(1);
	return ans;
}

/* 'x_seqs' must be the list, NOT the BStringList object!
   Return the list, NOT the BStringList object! */
SEXP subBStrings(SEXP x_seqs, SEXP start, SEXP nchar, SEXP proto)
{
	SEXP ans, ans_elt, x_seq, data;
	int nseq, offset, length, seq_length, i;

	offset = _start2offset(INTEGER(start)[0]);
	nseq = LENGTH(x_seqs);
	PROTECT(ans = NEW_LIST(nseq));
	for (i = 0; i < nseq; i++) {
		x_seq = VECTOR_ELT(x_seqs, i);
		_get_BString_charseq(x_seq, &seq_length);
		length = _nchar2length(INTEGER(nchar)[0], offset, seq_length);
		if (proto == R_NilValue) {
			PROTECT(ans_elt = duplicate(x_seq));
		} else {
			PROTECT(ans_elt = duplicate(proto));
			PROTECT(data = duplicate(GET_SLOT(x_seq, install("data"))));
			SET_SLOT(ans_elt, mkChar("data"), data);
			UNPROTECT(1);
		}
		SET_SLOT(ans_elt, mkChar("offset"), ScalarInteger(offset));
		SET_SLOT(ans_elt, mkChar("length"), ScalarInteger(length));
		SET_ELEMENT(ans, i, ans_elt);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	return ans;
}

/* 'x_seqs' must be the list, NOT the BStringList object! */
SEXP BStrings_to_nchars(SEXP x_seqs)
{
	SEXP ans, x_seq;
	int nseq, i, *seq_length;

	nseq = LENGTH(x_seqs);
	PROTECT(ans = NEW_INTEGER(nseq));
	for (i = 0, seq_length = INTEGER(ans); i < nseq; i++, seq_length++) {
		x_seq = VECTOR_ELT(x_seqs, i);
		_get_BString_charseq(x_seq, seq_length);
	}
	UNPROTECT(1);
	return ans;
}

/*
 * --- .Call ENTRY POINT ---
 * SEN_to_locs() arguments are assumed to be:
 *   seq_nchars: vector of non-negative integers (no NAs either)
 *   start, end, nchar: single integer (possibly NA)
 * SEN_to_locs() converts user-specified values 'start', 'end', 'nchar' into
 * valid start/nchar locations.
 * Return a list with 2 elements named "start" and "nchar", each of them being
 * an integer vector of the same length as 'seq_nchars'.
 */
SEXP SEN_to_locs(SEXP seq_nchars, SEXP start, SEXP end, SEXP nchar)
{
	SEXP ans_start, ans_nchar, ans, ans_names;
	int start0, end0, nchar0, nseq, i, *seq_nchar, *start_elt, *nchar_elt;

	start0 = INTEGER(start)[0];
	end0 = INTEGER(end)[0];
	nchar0 = INTEGER(nchar)[0];
	// Checking user-specified values 'start', 'end', 'nchar'
	if (start0 == 0)
		error("'start' must be a single >= 1, <= -1 or NA integer");
	if (end0 == 0)
		error("'end' must be a single >= 1, <= -1 or NA integer");
	if (nchar0 == NA_INTEGER) {
		if (start0 == NA_INTEGER)
			start0 = 1;
		if (end0 == NA_INTEGER)
			end0 = -1;
		if ((end0 > 0 || start0 < 0) && end0 < start0)
			error("invalid ('start','end') combination");
	} else if (nchar0 < 0) {
		error("'nchar' must be a single >= 0 or NA integer");
	} else if ((start0 == NA_INTEGER) == (end0 == NA_INTEGER)) {
		error("either 'start' or 'end' (but not both) must be NA when 'nchar' is not NA");
	} else if (start0 == NA_INTEGER) {
		// end0 is not NA
		if (0 < end0 && end0 < nchar0)
			error("invalid ('end','nchar') combination");
		start0 = end0 - nchar0 + 1; // will be 0 iff end0 = -1 and nchar0 = 0
	} else {
		// end0 is NA
		if (start0 < 0 && -start0 < nchar0)
			error("invalid ('start','nchar') combination");
		end0 = start0 + nchar0 - 1; // will be 0 iff start0 = 1 and nchar0 = 0
	}
	// From here, start0 and end0 can't be NA anymore so we don't
	// need nchar0 anymore.

	nseq = LENGTH(seq_nchars);
	PROTECT(ans_start = NEW_INTEGER(nseq));
	PROTECT(ans_nchar = NEW_INTEGER(nseq));
	for (i = 0, seq_nchar = INTEGER(seq_nchars),
		    start_elt = INTEGER(ans_start),
		    nchar_elt = INTEGER(ans_nchar);
	     i < nseq;
	     i++, seq_nchar++, start_elt++, nchar_elt++)
	{
		if (start0 > 0)
			*start_elt = start0;
		else
			*start_elt = *seq_nchar + start0 + 1;
		if (end0 >= 0)
			*nchar_elt = end0 - *start_elt + 1;
		else
			*nchar_elt = *seq_nchar + end0 + 1 - *start_elt + 1;
		if (*start_elt < 1) {
			UNPROTECT(2);
			error("trying to read before the start of input sequence %d", i+1);
		}
		if (*nchar_elt < 0) {
			UNPROTECT(2);
			error("trying to read a negative number of letters from input sequence %d", i+1);
		}
		if (*start_elt + *nchar_elt - 1 > *seq_nchar) {
			UNPROTECT(2);
			error("trying to read after the end of input sequence %d", i+1);
		}
	}

	PROTECT(ans = NEW_LIST(2));
	PROTECT(ans_names = NEW_CHARACTER(2));
	SET_STRING_ELT(ans_names, 0, mkChar("start"));
	SET_STRING_ELT(ans_names, 1, mkChar("nchar"));
	SET_NAMES(ans, ans_names);
	SET_ELEMENT(ans, 0, ans_start);
	SET_ELEMENT(ans, 1, ans_nchar);
	UNPROTECT(4);
	return ans;
}

SEXP get_start_for_adjacent_seqs(SEXP seq_nchars)
{
	SEXP ans;
	int nseq, i, *seq_nchar, *ans_elt0, *ans_elt1;

	nseq = LENGTH(seq_nchars);
	PROTECT(ans = NEW_INTEGER(nseq));
	if (nseq >= 1)
		INTEGER(ans)[0] = 1;
	if (nseq >= 2)
		for (i = 1, seq_nchar = INTEGER(seq_nchars),
			    ans_elt0 = INTEGER(ans),
			    ans_elt1 = INTEGER(ans)+1;
		     i < nseq;
		     i++, seq_nchar++, ans_elt0++, ans_elt1++) {
			*ans_elt1 = *ans_elt0 + *seq_nchar;
		}
	UNPROTECT(1);
	return ans;
}

