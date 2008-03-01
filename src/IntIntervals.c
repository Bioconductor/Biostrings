/****************************************************************************
 *                  Fast IntIntervals objects manipulation                  *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"

typedef struct startend {
	int start;
	int end;
} StartEnd;

static int debug = 0;

SEXP Biostrings_debug_IntIntervals()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'IntIntervals.c'\n", debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'IntIntervals.c'\n");
#endif
	return R_NilValue;
}

/*
 * Check and simplify user-specified values 'start', 'end', 'width'.
 */
static StartEnd uSEW_to_StartEnd(int start, int end, int width)
{
	StartEnd startend;

	if (start == 0)
		error("'start' must be a single >= 1, <= -1 or NA integer");
	if (end == 0)
		error("'end' must be a single >= 1, <= -1 or NA integer");
	if (width == NA_INTEGER) {
		if (start == NA_INTEGER)
			start = 1;
		if (end == NA_INTEGER)
			end = -1;
		if ((end > 0 || start < 0) && end < start)
			error("invalid ('start','end') combination");
	} else if (width < 0) {
		error("'width' must be a single >= 0 or NA integer");
	} else if ((start == NA_INTEGER) == (end == NA_INTEGER)) {
		error("either 'start' or 'end' (but not both) must be NA when 'width' is not NA");
	} else if (start == NA_INTEGER) {
		// 'end' is not NA
		if (0 < end && end < width)
			error("invalid ('end','width') combination");
		start = end - width + 1; // will be 0 iff 'end' = -1 and 'width' = 0
	} else {
		// 'end' is NA
		if (start < 0 && -start < width)
			error("invalid ('start','width') combination");
		end = start + width - 1; // will be 0 iff 'start' = 1 and 'width' = 0
	}
	// 'start' and 'end' cannot be NA anymore!
	startend.start = start;
	startend.end = end;
	return startend;
}

int _get_IntIntervals_length(SEXP x)
{
	SEXP inters;

	inters = GET_SLOT(x, install("inters"));
	return LENGTH(VECTOR_ELT(inters, 0));
}

const int *_get_IntIntervals_start(SEXP x)
{
	SEXP inters;

	inters = GET_SLOT(x, install("inters"));
	return INTEGER(VECTOR_ELT(inters, 0));
}

const int *_get_IntIntervals_width(SEXP x)
{
	SEXP inters;

	inters = GET_SLOT(x, install("inters"));
	return INTEGER(VECTOR_ELT(inters, 1));
}

/*
 * --- .Call ENTRY POINT ---
 */
SEXP narrow_IntIntervals(SEXP x, SEXP start, SEXP end, SEXP width)
{
	StartEnd startend;
	const int *old_start, *old_width;
	int x_length, i, *new_start, *new_width, shift1, shift2;
	SEXP ans_start, ans_width, ans, ans_names;

	startend = uSEW_to_StartEnd(INTEGER(start)[0], INTEGER(end)[0], INTEGER(width)[0]);
	x_length = _get_IntIntervals_length(x);
	PROTECT(ans_start = NEW_INTEGER(x_length));
	PROTECT(ans_width = NEW_INTEGER(x_length));
	for (i = 0, old_start = _get_IntIntervals_start(x),
		    old_width = _get_IntIntervals_width(x),
		    new_start = INTEGER(ans_start),
		    new_width = INTEGER(ans_width);
	     i < x_length;
	     i++, old_start++, old_width++, new_start++, new_width++)
	{
		if (startend.start > 0)
			shift1 = startend.start - 1;
		else
			shift1 = startend.start + *old_width;
		if (shift1 < 0) {
			UNPROTECT(2);
			error("cannot narrow interval %d, this would require moving "
			      "its 'start' (%d) to the left", i + 1, *old_start);
		}
		if (startend.end < 0)
			shift2 = startend.end + 1;
		else
			shift2 = startend.end - *old_width;
		if (shift2 > 0) {
			UNPROTECT(2);
			error("cannot narrow interval %d, this would require moving "
			      "its 'end' (%d) to the right", i + 1, *old_start + *old_width - 1);
		}
		*new_width = *old_width - shift1 + shift2;
		if (*new_width < 0) {
			UNPROTECT(2);
			error("cannot narrow interval %d, its 'width' (%d) is too small",
			      i + 1, *old_width);
		}
		*new_start = *old_start + shift1;
	}
	PROTECT(ans = NEW_LIST(2));
	PROTECT(ans_names = NEW_CHARACTER(2));
	SET_STRING_ELT(ans_names, 0, mkChar("start"));
	SET_STRING_ELT(ans_names, 1, mkChar("width"));
	SET_NAMES(ans, ans_names);
	SET_ELEMENT(ans, 0, ans_start);
	SET_ELEMENT(ans, 1, ans_width);
	UNPROTECT(4);
	return ans;
}

/*
 * --- .Call ENTRY POINT ---
 */
SEXP int_to_adjacent_intervals(SEXP x)
{
	SEXP ans;
	int nseq, i, *x_elt, *ans_elt0, *ans_elt1;

	nseq = LENGTH(x);
	PROTECT(ans = NEW_INTEGER(nseq));
	if (nseq >= 1)
		INTEGER(ans)[0] = 1;
	if (nseq >= 2)
		for (i = 1, x_elt = INTEGER(x),
			    ans_elt0 = INTEGER(ans),
			    ans_elt1 = INTEGER(ans)+1;
		     i < nseq;
		     i++, x_elt++, ans_elt0++, ans_elt1++) {
			*ans_elt1 = *ans_elt0 + *x_elt;
		}
	UNPROTECT(1);
	return ans;
}

