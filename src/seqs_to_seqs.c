/****************************************************************************
 *      Converting a set of sequences from one internal representation      *
 *                               into another                               *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"
#include <S.h> /* for Salloc() and Srealloc() */

static int debug = 0;

SEXP Biostrings_debug_seqs_to_seqs()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'seqs_to_seqs.c'\n", debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'seqs_to_seqs.c'\n");
#endif
	return R_NilValue;
}


/****************************************************************************
 * The SEN (Start/End/Nchar) interface.
 */

typedef struct startend {
	int start;
	int end;
} StartEnd;

/*
 * Check and simplify user-specified values 'start', 'end', 'nchar'.
 */
static StartEnd SEN_to_StartEnd(int start, int end, int nchar)
{
	StartEnd startend;

	if (start == 0)
		error("'start' must be a single >= 1, <= -1 or NA integer");
	if (end == 0)
		error("'end' must be a single >= 1, <= -1 or NA integer");
	if (nchar == NA_INTEGER) {
		if (start == NA_INTEGER)
			start = 1;
		if (end == NA_INTEGER)
			end = -1;
		if ((end > 0 || start < 0) && end < start)
			error("invalid ('start','end') combination");
	} else if (nchar < 0) {
		error("'nchar' must be a single >= 0 or NA integer");
	} else if ((start == NA_INTEGER) == (end == NA_INTEGER)) {
		error("either 'start' or 'end' (but not both) must be NA when 'nchar' is not NA");
	} else if (start == NA_INTEGER) {
		// 'end' is not NA
		if (0 < end && end < nchar)
			error("invalid ('end','nchar') combination");
		start = end - nchar + 1; // will be 0 iff 'end' = -1 and 'nchar' = 0
	} else {
		// 'end' is NA
		if (start < 0 && -start < nchar)
			error("invalid ('start','nchar') combination");
		end = start + nchar - 1; // will be 0 iff 'start' = 1 and 'nchar' = 0
	}
	// 'start' and 'end' cannot be NA anymore!
	startend.start = start;
	startend.end = end;
	return startend;
}

/*
 * --- .Call ENTRY POINT ---
 * SEN_to_locs() arguments are assumed to be:
 *   start, end, nchar: single integer (possibly NA)
 *   seq_nchars: vector of non-negative integers (no NAs either)
 * SEN_to_locs() converts user-specified values 'start', 'end', 'nchar' into
 * valid start/nchar locations.
 * Return a list with 2 elements named "start" and "nchar", each of them being
 * an integer vector of the same length as 'seq_nchars'.
 */
SEXP SEN_to_locs(SEXP start, SEXP end, SEXP nchar, SEXP seq_nchars)
{
	SEXP ans_start, ans_nchar, ans, ans_names;
	StartEnd startend;
	int nseq, i, *seq_nchar, *ans_start_elt, *ans_nchar_elt;

	startend = SEN_to_StartEnd(INTEGER(start)[0], INTEGER(end)[0], INTEGER(nchar)[0]);
	nseq = LENGTH(seq_nchars);
	PROTECT(ans_start = NEW_INTEGER(nseq));
	PROTECT(ans_nchar = NEW_INTEGER(nseq));
	for (i = 0, seq_nchar = INTEGER(seq_nchars),
		    ans_start_elt = INTEGER(ans_start),
		    ans_nchar_elt = INTEGER(ans_nchar);
	     i < nseq;
	     i++, seq_nchar++, ans_start_elt++, ans_nchar_elt++)
	{
		*ans_start_elt = startend.start;
		if (startend.start <= 0) // yes, <= 0
			*ans_start_elt += *seq_nchar + 1;
		*ans_nchar_elt = startend.end - *ans_start_elt + 1;
		if (startend.end < 0) // yes, < 0
			*ans_nchar_elt += *seq_nchar + 1;
		if (*ans_start_elt < 1) {
			UNPROTECT(2);
			error("trying to read before the start of input sequence %d", i+1);
		}
		if (*ans_nchar_elt < 0) {
			UNPROTECT(2);
			error("trying to read a negative number of letters from input sequence %d", i+1);
		}
		if (*ans_start_elt + *ans_nchar_elt - 1 > *seq_nchar) {
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

/*
 * --- .Call ENTRY POINT ---
 */
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


/****************************************************************************
 * 2 helper functions
 */

static int start2offset(int start)
{
	if (start < 1)
		error("'start' must be >= 1");
	return --start;
}

static int nchar2length(int nchar, int offset, int seq_length)
{
	if (nchar == NA_INTEGER) {
		nchar = seq_length - offset;
		if (nchar < 0L)
			error("cannot read a negative number of letters");
	} else {
		if (nchar < 0L)
			error("cannot read a negative number of letters");
		if (offset + nchar > seq_length)
			error("cannot read beyond the end of the input sequence");
	}
	return nchar;
}


/****************************************************************************
 * Converting a set of sequences into an array of CharSeq structs.
 *
 * The 3 functions below assume that 'start' and 'nchar' have the same
 * length as 'x'.
 */

const CharSeq *STRSXP_to_charseqs(SEXP x, const int *start, const int *nchar, int *nseq)
{
	CharSeq *seqs, *seq;
	int i, offset;
	const int *start_p, *nchar_p;
	SEXP x_elt;

	*nseq = LENGTH(x);
	seqs = Salloc((long) *nseq, CharSeq);
	for (i = 0, start_p = start, nchar_p = nchar, seq = seqs;
	     i < *nseq;
	     i++, start_p++, nchar_p++, seq++) {
		x_elt = STRING_ELT(x, i);
		if (x_elt == NA_STRING)
			error("input sequence %d is NA", i+1);
		seq->data = CHAR(x_elt);
		seq->data += offset = start2offset(*start_p);
		seq->length = nchar2length(*nchar_p, offset, LENGTH(x_elt));
	}
	return seqs;
}

const CharSeq *BStringList_to_charseqs(SEXP x, const int *start, const int *nchar, int *nseq)
{
	CharSeq *seqs, *seq;
	int i, offset;
	const int *start_p, *nchar_p;

	*nseq = _get_BStringList_length(x);
	seqs = Salloc((long) *nseq, CharSeq);
	for (i = 0, start_p = start, nchar_p = nchar, seq = seqs;
	     i < *nseq;
	     i++, start_p++, nchar_p++, seq++) {
		seq->data = _get_BStringList_charseq(x, i, &(seq->length));
		seq->data += offset = start2offset(*start_p);
		seq->length = nchar2length(*nchar_p, offset, seq->length);
	}
	return seqs;
}

const CharSeq *BStringSet_to_charseqs(SEXP x, const int *start, const int *nchar, int *nseq)
{
	CharSeq *seqs, *seq;
	int i, offset;
	const int *start_p, *nchar_p;

	*nseq = _get_BStringSet_length(x);
	seqs = Salloc((long) *nseq, CharSeq);
	for (i = 0, start_p = start, nchar_p = nchar, seq = seqs;
	     i < *nseq;
	     i++, start_p++, nchar_p++, seq++) {
		seq->data = _get_BStringSet_charseq(x, i, &(seq->length));
		seq->data += offset = start2offset(*start_p);
		seq->length = nchar2length(*nchar_p, offset, seq->length);
	}
	return seqs;
}


/****************************************************************************
 * Converting a set of sequences into an XRaw object.
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

/*
 * --- .Call ENTRY POINT ---
 */
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

/*
 * --- .Call ENTRY POINT ---
 * TODO: Support the 'collapse' argument
 */
SEXP STRSXP_to_XRaw(SEXP x, SEXP start, SEXP nchar, SEXP collapse, SEXP lkup)
{
	int nseq;
	const CharSeq *seqs;
	SEXP tag, ans;

	nseq = LENGTH(x);
	if (LENGTH(start) != nseq || LENGTH(nchar) != nseq)
		error("invalid length of 'start' or 'nchar'");
	if (collapse == R_NilValue) {
		if (nseq != 1)
			error("'collapse' must be specified when the number of input sequences is not exactly 1");
	} else {
		if (LENGTH(collapse) != 1 || LENGTH(STRING_ELT(collapse, 0)) != 0)
			error("'collapse' can only be NULL or the empty string for now");
	}
	seqs = STRSXP_to_charseqs(x, INTEGER(start), INTEGER(nchar), &nseq);
	PROTECT(tag = charseqs_to_RAW(seqs, nseq, lkup));
	ans = mkXRaw(tag);
	UNPROTECT(1);
	return ans;
}

/*
 * --- .Call ENTRY POINT ---
 */
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
 * Converting a set of sequences into a BStringList object.
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

/*
 * --- .Call ENTRY POINT ---
 */
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
		offset = start2offset(*start_p);
		length = nchar2length(*nchar_p, offset, x_length);
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

/*
 * --- .Call ENTRY POINT ---
 */
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
		offset = start2offset(*start_p);
		length = nchar2length(*nchar_p, offset, seq_length);
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

