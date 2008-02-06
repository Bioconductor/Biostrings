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
SEXP charseqs_to_BStringList(SEXP seqs, SEXP start, SEXP nchar, SEXP lkup, SEXP proto)
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

