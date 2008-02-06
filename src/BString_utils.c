/*
 * Low-level manipulation of BString and BStringList objects.
 */
#include "Biostrings.h"


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

const char *getBString_charseq(SEXP x, int *length)
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
		getBString_charseq(x_seq, &seq_length);
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
		getBString_charseq(x_seq, seq_length);
	}
	UNPROTECT(1);
	return ans;
}

