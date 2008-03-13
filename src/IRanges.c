/****************************************************************************
 *                    Fast IRanges objects manipulation                     *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"

static int debug = 0;

SEXP Biostrings_debug_IRanges()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'IRanges.c'\n", debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'IRanges.c'\n");
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

int _get_IRanges_length(SEXP x)
{
	SEXP ranges;

	ranges = GET_SLOT(x, install("ranges"));
	return LENGTH(VECTOR_ELT(ranges, 0));
}

const int *_get_IRanges_start(SEXP x)
{
	SEXP ranges;

	ranges = GET_SLOT(x, install("ranges"));
	return INTEGER(VECTOR_ELT(ranges, 0));
}

const int *_get_IRanges_width(SEXP x)
{
	SEXP ranges;

	ranges = GET_SLOT(x, install("ranges"));
	return INTEGER(VECTOR_ELT(ranges, 1));
}

static SEXP mk_ranges_slot(SEXP start, SEXP width, SEXP names)
{
	SEXP ans, ans_names, ans_row_names, ans_attribs, ans_attrib_names;
	int i;

	PROTECT(ans_attribs = NEW_LIST(3));

	/* set the attrib names */
	PROTECT(ans_attrib_names = NEW_CHARACTER(3));
	SET_STRING_ELT(ans_attrib_names, 0, mkChar("names"));
	SET_STRING_ELT(ans_attrib_names, 1, mkChar("row.names"));
	SET_STRING_ELT(ans_attrib_names, 2, mkChar("class"));
	SET_NAMES(ans_attribs, ans_attrib_names);
	UNPROTECT(1);

	/* set the class */
	SET_ELEMENT(ans_attribs, 2, mkString("data.frame"));

	/* set the row names */
	PROTECT(ans_row_names = NEW_INTEGER(LENGTH(start)));
	for (i = 0; i < LENGTH(ans_row_names); i++)
		INTEGER(ans_row_names)[i] = i + 1;
	SET_ELEMENT(ans_attribs, 1, ans_row_names);
	UNPROTECT(1);

	/* set the elements and their names */
	if (names == R_NilValue) {
		PROTECT(ans = NEW_LIST(2));
		PROTECT(ans_names = NEW_CHARACTER(2));
	} else {
		PROTECT(ans = NEW_LIST(3));
		PROTECT(ans_names = NEW_CHARACTER(3));
	}
	SET_ELEMENT(ans, 0, start);
	SET_STRING_ELT(ans_names, 0, mkChar("start"));
	SET_ELEMENT(ans, 1, width);
	SET_STRING_ELT(ans_names, 1, mkChar("width"));
	if (names != R_NilValue) {
		SET_ELEMENT(ans, 2, names);
		SET_STRING_ELT(ans_names, 2, mkChar("names"));
	}
	SET_ELEMENT(ans_attribs, 0, ans_names);
	UNPROTECT(1);

	SET_ATTRIB(ans, ans_attribs);
	UNPROTECT(2);
	return ans;
}

/*
 * Do NOT make this a .Call() entry point!
 * Its arguments are NOT duplicated so it would be a disaster if they were
 * coming from the user space.
 */
SEXP _new_IRanges(SEXP start, SEXP width, SEXP names)
{
	SEXP class_def, ans, ans_slot;

	class_def = MAKE_CLASS(".IRanges");
	PROTECT(ans = NEW_OBJECT(class_def));
	PROTECT(ans_slot = mk_ranges_slot(start, width, names));
	SET_SLOT(ans, mkChar("ranges"), ans_slot);
	UNPROTECT(2);
	return ans;
}

/*
 * --- .Call ENTRY POINT ---
 */
SEXP narrow_IRanges(SEXP x, SEXP start, SEXP end, SEXP width)
{
	StartEnd startend;
	const int *old_start, *old_width;
	int x_length, i, *new_start, *new_width, shift1, shift2;
	SEXP ans_start, ans_width, ans, ans_names;

	startend = uSEW_to_StartEnd(INTEGER(start)[0], INTEGER(end)[0], INTEGER(width)[0]);
	x_length = _get_IRanges_length(x);
	PROTECT(ans_start = NEW_INTEGER(x_length));
	PROTECT(ans_width = NEW_INTEGER(x_length));
	for (i = 0, old_start = _get_IRanges_start(x),
		    old_width = _get_IRanges_width(x),
		    new_start = INTEGER(ans_start),
		    new_width = INTEGER(ans_width);
	     i < x_length;
	     i++, old_start++, old_width++, new_start++, new_width++)
	{
		if (startend.start > 0)
			shift1 = startend.start - 1;
		else
			shift1 = startend.start + *old_width;
		if (startend.end < 0)
			shift2 = startend.end + 1;
		else
			shift2 = startend.end - *old_width;
		*new_width = *old_width - shift1 + shift2;
		if (shift1 < 0 || shift2 > 0 || *new_width < 0) {
			UNPROTECT(2);
			error("width of range %d is too small (%d) for this narrowing",
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
SEXP int_to_adjacent_ranges(SEXP x)
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
 * Reduction (aka extracting the frame)
 */

static RangeBuf reduced_ranges;
static int max_end, inframe_offset;

static void add_to_reduced_ranges(int start, int width)
{
	int buf_length, end, gap;

	buf_length = reduced_ranges.start.nelt;
	end = start + width - 1;
	if (buf_length == 0 || (gap = start - max_end - 1) > 0) {
		_RangeBuf_insert_at(&reduced_ranges, buf_length, start, width);
		if (buf_length == 0)
			inframe_offset = start - 1;
		else
			inframe_offset += gap;
		max_end = end;
		return;
	}
	if (end <= max_end)
		return;
	reduced_ranges.width.elts[buf_length - 1] += end - max_end;
	max_end = end;
	return;
}

static void reduce_ranges(int length, const int *start, const int *width, int *inframe_start)
{
	int i, j;
	IntBuf start_order;

	start_order = _new_IntBuf(length, 0);
	get_intorder(length, start, start_order.elts);
	reduced_ranges = _new_RangeBuf(0, 0);
	for (i = 0; i < length; i++) {
		j = start_order.elts[i];
		add_to_reduced_ranges(start[j], width[j]);
		if (inframe_start != NULL)
			inframe_start[j] = start[j] - inframe_offset;
	}
	return;
}

/*
 * --- .Call ENTRY POINT ---
 */
SEXP reduce_IRanges(SEXP x, SEXP with_inframe_start)
{
	int x_length;
	const int *x_start, *x_width;
	SEXP ans, ans_names, ans_inframe_start;
	int *inframe_start;

	x_length = _get_IRanges_length(x);
	x_start = _get_IRanges_start(x);
	x_width = _get_IRanges_width(x);
	if (LOGICAL(with_inframe_start)[0]) {
		PROTECT(ans_inframe_start = NEW_INTEGER(x_length));
		inframe_start = INTEGER(ans_inframe_start);
	} else {
		inframe_start = NULL;
	}
	reduce_ranges(x_length, x_start, x_width, inframe_start);

	PROTECT(ans = NEW_LIST(3));
	PROTECT(ans_names = NEW_CHARACTER(3));
	SET_STRING_ELT(ans_names, 0, mkChar("start"));
	SET_STRING_ELT(ans_names, 1, mkChar("width"));
	SET_STRING_ELT(ans_names, 2, mkChar("inframe.start"));
	SET_NAMES(ans, ans_names);
	UNPROTECT(1);
	SET_ELEMENT(ans, 0, _IntBuf_asINTEGER(&(reduced_ranges.start)));
	SET_ELEMENT(ans, 1, _IntBuf_asINTEGER(&(reduced_ranges.width)));
	if (inframe_start != NULL) {
		SET_ELEMENT(ans, 2, ans_inframe_start);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	return ans;
}

