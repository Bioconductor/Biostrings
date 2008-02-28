/****************************************************************************
 *     Turning a set of sequences into an XRaw or a BStringList object      *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"

static int debug = 0;

SEXP Biostrings_debug_seqs_to_SEXP()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'seqs_to_SEXP.c'\n", debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'seqs_to_SEXP.c'\n");
#endif
	return R_NilValue;
}


/****************************************************************************
 * ... into an XRaw object
 */

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
		error("invalid length of 'start' or 'nchar'");
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
		error("invalid length of 'start' or 'nchar'");
	seqs = BStringSet_to_charseqs(x, INTEGER(start), INTEGER(nchar), &nseq);
	PROTECT(tag = charseqs_to_RAW(seqs, nseq, lkup));
	ans = mkXRaw(tag);
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * ... into a BStringList object
 */

/* NOT a Call() entry point! */
SEXP mkBStringList(const char *class, SEXP seqs)
{
	SEXP class_def, ans;

	class_def = MAKE_CLASS(class);
	PROTECT(ans = NEW_OBJECT(class_def));
	SET_SLOT(ans, mkChar("seqs"), seqs);
	UNPROTECT(1);
	return ans;
}

SEXP XRaw_to_BStringList(SEXP x, SEXP start, SEXP nchar, SEXP proto)
{
	int nseq, x_length, i, offset, length;
	SEXP ans, ans_seqs, ans_seq, data;
	const int *start_p, *nchar_p;
	char classbuf[14]; // longest string will be "DNAStringList"

	if (LENGTH(start) != LENGTH(nchar))
		error("'start' and 'nchar' must have the same length");
	nseq = LENGTH(start);
	x_length = LENGTH(R_ExternalPtrTag(GET_SLOT(x, install("xp"))));
	PROTECT(ans_seqs = NEW_LIST(nseq));
	for (i = 0, start_p = INTEGER(start), nchar_p = INTEGER(nchar);
	     i < nseq;
	     i++, start_p++, nchar_p++) {
		PROTECT(ans_seq = duplicate(proto));
		PROTECT(data = duplicate(x));
		offset = _start2offset(*start_p);
		length = _nchar2length(*nchar_p, offset, x_length);
		SET_SLOT(ans_seq, mkChar("data"), data);
		SET_SLOT(ans_seq, mkChar("offset"), ScalarInteger(offset));
		SET_SLOT(ans_seq, mkChar("length"), ScalarInteger(length));
		SET_ELEMENT(ans_seqs, i, ans_seq);
		UNPROTECT(2);
	}
	snprintf(classbuf, sizeof(classbuf), "%sList", get_class(proto));
	ans = mkBStringList(classbuf, ans_seqs);
	UNPROTECT(1);
	return ans;
}

SEXP narrow_BStringList(SEXP x, SEXP start, SEXP nchar, SEXP proto)
{
	int nseq, i, seq_length, offset, length;
	SEXP x_seqs, x_seq, ans, ans_seqs, ans_seq, data;
	const int *start_p, *nchar_p;
	char classbuf[14]; // longest string will be "DNAStringList"
	const char *class;

	nseq = _get_BStringList_length(x);
	if (LENGTH(start) != nseq || LENGTH(nchar) != nseq)
		error("invalid length of 'start' or 'nchar'");
	x_seqs = GET_SLOT(x, install("seqs"));
	PROTECT(ans_seqs = NEW_LIST(nseq));
	for (i = 0, start_p = INTEGER(start), nchar_p = INTEGER(nchar);
	     i < nseq;
	     i++, start_p++, nchar_p++) {
	        x_seq = VECTOR_ELT(x_seqs, i);
		_get_BString_charseq(x_seq, &seq_length);
		if (proto == R_NilValue) {
			PROTECT(ans_seq = duplicate(x_seq));
		} else {
			PROTECT(ans_seq = duplicate(proto));
			PROTECT(data = duplicate(GET_SLOT(x_seq, install("data"))));
			SET_SLOT(ans_seq, mkChar("data"), data);
			UNPROTECT(1);
		}
		offset = _start2offset(*start_p);
		length = _nchar2length(*nchar_p, offset, seq_length);
		SET_SLOT(ans_seq, mkChar("offset"), ScalarInteger(offset));
		SET_SLOT(ans_seq, mkChar("length"), ScalarInteger(length));
		SET_ELEMENT(ans_seqs, i, ans_seq);
		UNPROTECT(1);
	}
	if (proto == R_NilValue) {
		class = get_class(x);
	} else {
		snprintf(classbuf, sizeof(classbuf), "%sList", get_class(proto));
		class = classbuf;
	}
	ans = mkBStringList(class, ans_seqs);
	UNPROTECT(1);
	return ans;
}

