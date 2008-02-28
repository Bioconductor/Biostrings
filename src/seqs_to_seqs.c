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
 * SEN_to_safelocs() arguments are assumed to be:
 *   start, end, nchar: single integer (possibly NA)
 *   seq_nchars: vector of non-negative integers (no NAs either)
 * Create a set of valid locations in 'x' from user-specified Start/End/Nchar
 * values ('start', 'end', 'nchar').
 * Return a list with 2 elements named "start" and "nchar", each of
 * them being an integer vector of the same length as 'seq_nchars'. Those
 * vectors are _safe_ i.e. they describe a set of valid locations in the
 * sequences whose lengths are passed to 'seq_nchars'.
 */
SEXP SEN_to_safelocs(SEXP start, SEXP end, SEXP nchar, SEXP seq_nchars)
{
	SEXP safe_starts, safe_nchars, ans, ans_names;
	StartEnd startend;
	int nseq, i, *seq_nchar, *start_p, *nchar_p;

	startend = SEN_to_StartEnd(INTEGER(start)[0], INTEGER(end)[0], INTEGER(nchar)[0]);
	nseq = LENGTH(seq_nchars);
	PROTECT(safe_starts = NEW_INTEGER(nseq));
	PROTECT(safe_nchars = NEW_INTEGER(nseq));
	for (i = 0, seq_nchar = INTEGER(seq_nchars),
		    start_p = INTEGER(safe_starts),
		    nchar_p = INTEGER(safe_nchars);
	     i < nseq;
	     i++, seq_nchar++, start_p++, nchar_p++)
	{
		*start_p = startend.start;
		if (startend.start <= 0) // yes, <= 0
			*start_p += *seq_nchar + 1;
		*nchar_p = startend.end - *start_p + 1;
		if (startend.end < 0) // yes, < 0
			*nchar_p += *seq_nchar + 1;
		if (*start_p < 1) {
			UNPROTECT(2);
			error("trying to read before the start of input sequence %d", i+1);
		}
		if (*nchar_p < 0) {
			UNPROTECT(2);
			error("trying to read a negative number of letters from input sequence %d", i+1);
		}
		if (*start_p + *nchar_p - 1 > *seq_nchar) {
			UNPROTECT(2);
			error("trying to read after the end of input sequence %d", i+1);
		}
	}
	PROTECT(ans = NEW_LIST(2));
	PROTECT(ans_names = NEW_CHARACTER(2));
	SET_STRING_ELT(ans_names, 0, mkChar("start"));
	SET_STRING_ELT(ans_names, 1, mkChar("nchar"));
	SET_NAMES(ans, ans_names);
	SET_ELEMENT(ans, 0, safe_starts);
	SET_ELEMENT(ans, 1, safe_nchars);
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
 * This helper function should never raise an error.
 *
 * All the functions in the rest of this file assume that their 'safe_starts'
 * and 'safe_nchars' arguments have been obtained by processing the
 * user-specified Start/End/Nchar values thru SEN_to_safelocs().
 * Hence those arguments are assumed to be safe i.e. to describe a set of
 * valid locations in 'x'.
 */

static int start2offset(int safe_start)
{
	if (safe_start < 1)
		error("Biostrings internal error in start2offset(): safe_start < 1");
	return --safe_start;
}


/****************************************************************************
 * Converting a set of sequences into an array of CharSeq structs.
 */

const CharSeq *STRSXP_to_charseqs(SEXP x,
		const int *safe_starts, const int *safe_nchars, int *nseq)
{
	CharSeq *seqs, *seq;
	int i;
	const int *start_p, *nchar_p;
	SEXP x_elt;

	*nseq = LENGTH(x);
	seqs = Salloc((long) *nseq, CharSeq);
	for (i = 0, start_p = safe_starts, nchar_p = safe_nchars, seq = seqs;
	     i < *nseq;
	     i++, start_p++, nchar_p++, seq++) {
		x_elt = STRING_ELT(x, i);
		if (x_elt == NA_STRING)
			error("input sequence %d is NA", i+1);
		seq->data = CHAR(x_elt);
		seq->data += start2offset(*start_p);
		seq->length = *nchar_p;
	}
	return seqs;
}

const CharSeq *BStringList_to_charseqs(SEXP x,
		const int *safe_starts, const int *safe_nchars, int *nseq)
{
	CharSeq *seqs, *seq;
	int i;
	const int *start_p, *nchar_p;

	*nseq = _get_BStringList_length(x);
	seqs = Salloc((long) *nseq, CharSeq);
	for (i = 0, start_p = safe_starts, nchar_p = safe_nchars, seq = seqs;
	     i < *nseq;
	     i++, start_p++, nchar_p++, seq++) {
		seq->data = _get_BStringList_charseq(x, i, &(seq->length));
		seq->data += start2offset(*start_p);
		seq->length = *nchar_p;
	}
	return seqs;
}

const CharSeq *BStringSet_to_charseqs(SEXP x,
		const int *safe_starts, const int *safe_nchars, int *nseq)
{
	CharSeq *seqs, *seq;
	int i;
	const int *start_p, *nchar_p;

	*nseq = _get_BStringSet_length(x);
	seqs = Salloc((long) *nseq, CharSeq);
	for (i = 0, start_p = safe_starts, nchar_p = safe_nchars, seq = seqs;
	     i < *nseq;
	     i++, start_p++, nchar_p++, seq++) {
		seq->data = _get_BStringSet_charseq(x, i, &(seq->length));
		seq->data += start2offset(*start_p);
		seq->length = *nchar_p;
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
SEXP STRSXP_to_XRaw(SEXP x, SEXP safe_starts, SEXP safe_nchars, SEXP collapse, SEXP lkup)
{
	int nseq;
	const CharSeq *seqs;
	SEXP tag, ans;

	nseq = LENGTH(x);
	if (LENGTH(safe_starts) != nseq || LENGTH(safe_nchars) != nseq)
		error("invalid length of 'safe_starts' or 'safe_nchars'");
	if (collapse == R_NilValue) {
		if (nseq != 1)
			error("'collapse' must be specified when the number of input sequences is not exactly 1");
	} else {
		if (LENGTH(collapse) != 1 || LENGTH(STRING_ELT(collapse, 0)) != 0)
			error("'collapse' can only be NULL or the empty string for now");
	}
	seqs = STRSXP_to_charseqs(x, INTEGER(safe_starts), INTEGER(safe_nchars), &nseq);
	PROTECT(tag = charseqs_to_RAW(seqs, nseq, lkup));
	ans = mkXRaw(tag);
	UNPROTECT(1);
	return ans;
}

/*
 * --- .Call ENTRY POINT ---
 */
SEXP BStringSet_to_XRaw(SEXP x, SEXP safe_starts, SEXP safe_nchars, SEXP lkup)
{
	int nseq;
	const CharSeq *seqs;
	SEXP tag, ans;

	nseq = _get_BStringSet_length(x);
	if (LENGTH(safe_starts) != nseq || LENGTH(safe_nchars) != nseq)
		error("invalid length of 'safe_starts' or 'safe_nchars'");
	seqs = BStringSet_to_charseqs(x, INTEGER(safe_starts), INTEGER(safe_nchars), &nseq);
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
SEXP XRaw_to_BStringList(SEXP x, SEXP safe_starts, SEXP safe_nchars, SEXP proto)
{
	int nseq, x_length, i;
	SEXP ans, ans_seqs, ans_seq, data;
	const int *start_p, *nchar_p;
	char classbuf[14]; // longest string will be "DNAStringList"

	if (LENGTH(safe_starts) != LENGTH(safe_nchars))
		error("'safe_starts' and 'safe_nchars' must have the same length");
	nseq = LENGTH(safe_starts);
	x_length = LENGTH(R_ExternalPtrTag(GET_SLOT(x, install("xp"))));
	PROTECT(ans_seqs = NEW_LIST(nseq));
	for (i = 0, start_p = INTEGER(safe_starts), nchar_p = INTEGER(safe_nchars);
	     i < nseq;
	     i++, start_p++, nchar_p++) {
		PROTECT(ans_seq = duplicate(proto));
		PROTECT(data = duplicate(x));
		SET_SLOT(ans_seq, mkChar("data"), data);
		SET_SLOT(ans_seq, mkChar("offset"), ScalarInteger(start2offset(*start_p)));
		SET_SLOT(ans_seq, mkChar("length"), ScalarInteger(*nchar_p));
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
SEXP narrow_BStringList(SEXP x, SEXP safe_starts, SEXP safe_nchars, SEXP proto)
{
	int nseq, i, seq_length;
	SEXP x_seqs, x_seq, ans, ans_seqs, ans_seq, data;
	const int *start_p, *nchar_p;
	char classbuf[14]; // longest string will be "DNAStringList"
	const char *class;

	nseq = _get_BStringList_length(x);
	if (LENGTH(safe_starts) != nseq || LENGTH(safe_nchars) != nseq)
		error("invalid length of 'safe_starts' or 'safe_nchars'");
	x_seqs = GET_SLOT(x, install("seqs"));
	PROTECT(ans_seqs = NEW_LIST(nseq));
	for (i = 0, start_p = INTEGER(safe_starts), nchar_p = INTEGER(safe_nchars);
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
		SET_SLOT(ans_seq, mkChar("offset"), ScalarInteger(start2offset(*start_p)));
		SET_SLOT(ans_seq, mkChar("length"), ScalarInteger(*nchar_p));
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

