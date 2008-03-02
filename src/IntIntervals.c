/****************************************************************************
 *                  Fast IntIntervals objects manipulation                  *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"

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


/****************************************************************************
 * Check and simplify user-specified values 'start', 'end', 'width'.
 */

typedef struct startend {
	int start;
	int end;
} StartEnd;

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
	int x_length, i, *x_elt, *ans_elt0, *ans_elt1;

	x_length = LENGTH(x);
	PROTECT(ans = NEW_INTEGER(x_length));
	if (x_length >= 1)
		INTEGER(ans)[0] = 1;
	if (x_length >= 2)
		for (i = 1, x_elt = INTEGER(x),
			    ans_elt0 = INTEGER(ans),
			    ans_elt1 = INTEGER(ans)+1;
		     i < x_length;
		     i++, x_elt++, ans_elt0++, ans_elt1++) {
			*ans_elt1 = *ans_elt0 + *x_elt;
		}
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * Normalization
 */

static InterBuf normalized_intervals;
static int max_end, inframe_offset;

static void add_to_normalized_intervals(int start, int width)
{
	int buf_length, end, gap;

	buf_length = normalized_intervals.start.count;
	end = start + width - 1;
	if (buf_length == 0 || (gap = start - max_end - 1) > 0) {
		_InterBuf_insert_at(&normalized_intervals, buf_length, start, width);
		if (buf_length == 0)
			inframe_offset = start - 1;
		else
			inframe_offset += gap;
		max_end = end;
		return;
	}
	if (end <= max_end)
		return;
	normalized_intervals.width.vals[buf_length - 1] += end - max_end;
	max_end = end;
	return;
}

static void normalize_intervals(int length, const int *start, const int *width, int *inframe_start)
{
	int i, j;
	IBuf start_order;

	_IBuf_init(&start_order, length, 0);
	get_intorder(length, start, start_order.vals);
	_InterBuf_init(&normalized_intervals, 0, 0);
	for (i = 0; i < length; i++) {
		j = start_order.vals[i];
		add_to_normalized_intervals(start[j], width[j]);
		if (inframe_start != NULL)
			inframe_start[j] = start[j] - inframe_offset;
	}
	return;
}

/*
 * --- .Call ENTRY POINT ---
 */
SEXP normalize_IntIntervals(SEXP x, SEXP with_inframe_start)
{
	int x_length;
	const int *x_start, *x_width;
	SEXP ans, ans_names, ans_inframe_start;
	int *inframe_start;

	x_length = _get_IntIntervals_length(x);
	x_start = _get_IntIntervals_start(x);
	x_width = _get_IntIntervals_width(x);
	if (LOGICAL(with_inframe_start)[0]) {
		PROTECT(ans_inframe_start = NEW_INTEGER(x_length));
		inframe_start = INTEGER(ans_inframe_start);
	} else {
		inframe_start = NULL;
	}
	normalize_intervals(x_length, x_start, x_width, inframe_start);

	PROTECT(ans = NEW_LIST(3));
	PROTECT(ans_names = NEW_CHARACTER(3));
	SET_STRING_ELT(ans_names, 0, mkChar("start"));
	SET_STRING_ELT(ans_names, 1, mkChar("width"));
	SET_STRING_ELT(ans_names, 2, mkChar("inframe.start"));
	SET_NAMES(ans, ans_names);
	UNPROTECT(1);
	SET_ELEMENT(ans, 0, _IBuf_asINTEGER(&(normalized_intervals.start)));
	SET_ELEMENT(ans, 1, _IBuf_asINTEGER(&(normalized_intervals.width)));
	if (inframe_start != NULL) {
		SET_ELEMENT(ans, 2, ans_inframe_start);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	return ans;
}

