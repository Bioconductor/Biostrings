/****************************************************************************
 *             Turning an SEXP object containing input sequences            *
 *                        into an array of CharSeq                          *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"
#include <S.h> /* for Salloc() and Srealloc() */

static int debug = 0;

SEXP Biostrings_debug_SEXP_to_charseqs()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'SEXP_to_charseqs.c'\n", debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'SEXP_to_charseqs.c'\n");
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
 * Turning a character vector ...
 *
 * 'start' and 'nchar' must have the same length as 'x'.
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
		seq->data += offset = _start2offset(*start_p);
		seq->length = _nchar2length(*nchar_p, offset, LENGTH(x_elt));
	}
	return seqs;
}


/****************************************************************************
 * Turning a BStringList object ...
 *
 * 'start' and 'nchar' must have the same length as 'x'.
 */

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
		seq->data += offset = _start2offset(*start_p);
		seq->length = _nchar2length(*nchar_p, offset, seq->length);
	}
	return seqs;
}


/****************************************************************************
 * Turning a BStringSet object ...
 *
 * 'start' and 'nchar' must have the same length as 'x'.
 */

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
		seq->data += offset = _start2offset(*start_p);
		seq->length = _nchar2length(*nchar_p, offset, seq->length);
	}
	return seqs;
}

