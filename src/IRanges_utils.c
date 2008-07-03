/****************************************************************************
 *                          Fast IRanges utilities                          *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"

static int debug = 0;

SEXP debug_IRanges_utils()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'IRanges_utils.c'\n",
	        debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'IRanges_utils.c'\n");
#endif
	return R_NilValue;
}

/*
 * --- .Call ENTRY POINT ---
 */
SEXP narrow_IRanges(SEXP x, SEXP start, SEXP end, SEXP width)
{
	int start0, end0, x_length, i, *new_start, *new_width, shift1, shift2, start_length;
	const int *old_start, *old_width;
	SEXP ans_start, ans_width, ans, ans_names;

	x_length = _get_IRanges_length(x);
	PROTECT(ans_start = NEW_INTEGER(x_length));
	PROTECT(ans_width = NEW_INTEGER(x_length));
	start_length = LENGTH(start);
	if (start_length == 1) {
		start0 = INTEGER(start)[0];
		end0 = INTEGER(end)[0];
		_normargs_startend(&start0, &end0, INTEGER(width)[0], "");
		for (i = 0, old_start = _get_IRanges_start0(x),
			    old_width = _get_IRanges_width0(x),
			    new_start = INTEGER(ans_start),
			    new_width = INTEGER(ans_width);
		     i < x_length;
		     i++, old_start++, old_width++, new_start++, new_width++)
		{
			if (start0 > 0)
				shift1 = start0 - 1;
			else
				shift1 = start0 + *old_width;
			if (end0 < 0)
				shift2 = end0 + 1;
			else
				shift2 = end0 - *old_width;
			*new_width = *old_width - shift1 + shift2;
			if (shift1 < 0 || shift2 > 0 || *new_width < 0) {
				UNPROTECT(2);
				error("width of range %d is too small (%d) for this narrowing",
				      i + 1, *old_width);
			}
			*new_start = *old_start + shift1;
		}
	} else {
		for (i = 0, old_start = _get_IRanges_start0(x),
			    old_width = _get_IRanges_width0(x),
			    new_start = INTEGER(ans_start),
			    new_width = INTEGER(ans_width);
		     i < x_length;
		     i++, old_start++, old_width++, new_start++, new_width++)
		{
			start0 = INTEGER(start)[i];
			end0 = INTEGER(end)[i];
			_normargs_startend(&start0, &end0, INTEGER(width)[i], "");
			if (start0 > 0)
				shift1 = start0 - 1;
			else
				shift1 = start0 + *old_width;
			if (end0 < 0)
				shift2 = end0 + 1;
			else
				shift2 = end0 - *old_width;
			*new_width = *old_width - shift1 + shift2;
			if (shift1 < 0 || shift2 > 0 || *new_width < 0) {
				UNPROTECT(2);
				error("width of range %d is too small (%d) for this narrowing",
				      i + 1, *old_width);
			}
			*new_start = *old_start + shift1;
		}
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

	start_order = _new_IntBuf(length, 0, 0);
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
	x_start = _get_IRanges_start0(x);
	x_width = _get_IRanges_width0(x);
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

/*
 * --- .Call ENTRY POINT ---
 */
SEXP IRanges_coverage(SEXP x, SEXP ans_length)
{
	int x_len, ans_len, *ans_elt, i, j;
	const int *x_start, *x_width;
	SEXP ans;

	ans_len = INTEGER(ans_length)[0];
	PROTECT(ans = NEW_INTEGER(ans_len));
	memset(INTEGER(ans), 0, ans_len * sizeof(int));
	x_len = _get_IRanges_length(x);
	for (i = 0, x_start = _get_IRanges_start0(x),
	            x_width = _get_IRanges_width0(x);
	     i < x_len;
	     i++, x_start++, x_width++)
	{
		for (j = 0, ans_elt = INTEGER(ans) + *x_start - 1;
		     j < *x_width;
		     j++, ans_elt++)
		{
			*ans_elt += 1;
		}
	}
	UNPROTECT(1);
	return ans;
}

