/****************************************************************************
 *              Turning a set of sequences into an XRaw object              *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"

static int debug = 0;

SEXP Biostrings_debug_seqs_to_XRaw()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'seqs_to_XRaw.c'\n", debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'seqs_to_XRaw.c'\n");
#endif
	return R_NilValue;
}

/* NOT a Call() entry point! */
SEXP mkXRaw(SEXP tag)
{
	SEXP ans;

	PROTECT(ans = NEW_OBJECT(MAKE_CLASS("XRaw")));
	SET_SLOT(ans, mkChar("xp"), R_MakeExternalPtr(NULL, tag, R_NilValue));
	UNPROTECT(1);
        return ans;
}

/* NOT a Call() entry point! */
SEXP CHARSXP_to_XRaw(SEXP x, SEXP start, SEXP nchar, SEXP lkup)
{
	SEXP ans, dest;
	int offset, length;
	const char *src;

	offset = _start2offset(INTEGER(start)[0]);
	length = _nchar2length(INTEGER(nchar)[0], offset, LENGTH(x));
	src = CHAR(x) + offset;

	PROTECT(dest = NEW_RAW(length));
	if (lkup == R_NilValue) {
		_Biostrings_memcpy_to_i1i2(0, length - 1,
			(char *) RAW(dest), length,
			src, length, sizeof(char));
	} else {
		_Biostrings_translate_charcpy_to_i1i2(0, length - 1,
			(char *) RAW(dest), length,
			src, length,
			INTEGER(lkup), LENGTH(lkup));
	}
	ans = mkXRaw(dest);
	UNPROTECT(1);
	return ans;
}

SEXP char_to_XRaw(SEXP x, SEXP start, SEXP nchar, SEXP lkup)
{
	return CHARSXP_to_XRaw(STRING_ELT(x, 0), start, nchar, lkup);
}

SEXP copy_subXRaw(SEXP x, SEXP start, SEXP nchar, SEXP lkup)
{
	SEXP ans;

	error("copy_subXRaw() not ready yet");
	return R_NilValue;
}

static SEXP charseqs_to_RAW(const CharSeq *seqs, int nseq, SEXP lkup)
{
	SEXP ans;
	int ans_length, i;
	const CharSeq *seq;
	char *dest;

	ans_length = 0;
	for (i = 0, seq = seqs; i < nseq; i++, seq++)
		ans_length += seq->length;
	PROTECT(ans = NEW_RAW(ans_length));
	dest = (char *) RAW(ans);
	for (i = 0, seq = seqs; i < nseq; i++, seq++) {
		if (lkup == R_NilValue) {
			_Biostrings_memcpy_to_i1i2(0, seq->length - 1,
				dest, seq->length,
				seq->data, seq->length, sizeof(char));
		} else {
			_Biostrings_translate_charcpy_to_i1i2(0, seq->length - 1,
				dest, seq->length,
				seq->data, seq->length,
				INTEGER(lkup), LENGTH(lkup));
		}
		dest += seq->length;
	}
	UNPROTECT(1);
	return ans;
}

SEXP STRSXP_to_XRaw(SEXP x, SEXP start, SEXP nchar, SEXP lkup)
{
	int nseq;
	const CharSeq *seqs;
	SEXP tag, ans;

	nseq = LENGTH(x);
	if (LENGTH(start) != nseq || LENGTH(nchar) != nseq)
		error("invalid length of 'start' or 'end'");
	seqs = STRSXP_to_charseqs(x, INTEGER(start), INTEGER(nchar), &nseq);
	PROTECT(tag = charseqs_to_RAW(seqs, nseq, lkup));
	ans = mkXRaw(tag);
	UNPROTECT(1);
	return ans;
}

SEXP BStringSet_to_XRaw(SEXP x, SEXP start, SEXP nchar, SEXP lkup)
{
	int nseq;
	const CharSeq *seqs;
	SEXP tag, ans;

	nseq = _get_BStringSet_length(x);
	if (LENGTH(start) != nseq || LENGTH(nchar) != nseq)
		error("invalid length of 'start' or 'end'");
	seqs = BStringSet_to_charseqs(x, INTEGER(start), INTEGER(nchar), &nseq);
	PROTECT(tag = charseqs_to_RAW(seqs, nseq, lkup));
	ans = mkXRaw(tag);
	UNPROTECT(1);
	return ans;
}

