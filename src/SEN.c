/****************************************************************************
 *                   The SEN (Start/End/Nchar) interface                    *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"

typedef struct startend {
	int start;
	int end;
} StartEnd;

static int debug = 0;

SEXP Biostrings_debug_SEN()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'SEN.c'\n", debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'SEN.c'\n");
#endif
	return R_NilValue;
}

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

