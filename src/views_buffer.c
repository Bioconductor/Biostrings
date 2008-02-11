/*
 * Manipulation of the "views buffer".
 *
 * The "views buffer" is used for temporarily storing a set of views or
 * matches found in a sequence. E.g. it is used by the string searching
 * functions like find_palindromes(), find_repeats(), the pattern matching
 * functions, etc...
 * Except for Biostrings_debug_views_buffer(), the functions defined in
 * this file are NOT .Call methods (but they are used by .Call methods
 * defined in other .c files) so THEY DON'T NEED TO BE REGISTERED in
 * R_init_Biostrings.c. They are prefixed with "_Biostrings_" to minimize
 * the risk of a clash with symbols defined elsewhere (e.g. in libc).
 */
#include "Biostrings.h"
#include <S.h> /* for Srealloc() */


static int debug = 0;

SEXP Biostrings_debug_views_buffer()
{
#ifdef DEBUG_BIOSTRINGS
        debug = !debug;
        Rprintf("Debug mode turned %s in 'views_buffer.c'\n", debug ? "on" : "off");
#else
        Rprintf("Debug mode not available in 'views_buffer.c'\n");
#endif
        return R_NilValue;
}

static int viewsbuf_reporting_mode;

static int *viewsbuf_start, *viewsbuf_end;
static char **viewsbuf_desc;
static int viewsbuf_maxcount, viewsbuf_count; /* viewsbuf_maxcount >= viewsbuf_count */

static int new_view()
{
	long new_maxcount;

	if (viewsbuf_count >= viewsbuf_maxcount) {
		/* Buffer is full */
		if (viewsbuf_maxcount == 0)
			new_maxcount = 1024;
		else
			new_maxcount = 2 * viewsbuf_maxcount;
		viewsbuf_start = Srealloc((char *) viewsbuf_start, new_maxcount,
						(long) viewsbuf_maxcount, int);
		viewsbuf_end = Srealloc((char *) viewsbuf_end, new_maxcount,
						(long) viewsbuf_maxcount, int);
		if (viewsbuf_reporting_mode == 0)
			viewsbuf_desc = Srealloc((char *) viewsbuf_desc, new_maxcount,
						(long) viewsbuf_maxcount, char *);
		viewsbuf_maxcount = new_maxcount;
	}
	return viewsbuf_count++;
}

/*
 * Must be used in reporting mode >= 2 (see _Biostrings_reset_viewsbuf() below).
 */
static void insert_view_at(int start, int end, int insert_at)
{
	int i, j;

	j = new_view(); // viewsbuf_count - 1
	i = j - 1;
	while (insert_at <= i) {
		viewsbuf_start[j] = viewsbuf_start[i];
		viewsbuf_end[j] = viewsbuf_end[i];
		i--;
		j--;
	}
	viewsbuf_start[j] = start;
	viewsbuf_end[j] = end;
	return;
}

/*
 * Used by reporting mode 3 (see _Biostrings_reset_viewsbuf() below).
 */
static void insert_view_if_new(int start, int end)
{
	int i;

	i = viewsbuf_count - 1;
	while (0 <= i && start < viewsbuf_start[i]) i--;
	while (0 <= i && start == viewsbuf_start[i] && end < viewsbuf_end[i]) i--;
	if (0 <= i && start == viewsbuf_start[i] && end == viewsbuf_end[i])
		return;
	i++;
	insert_view_at(start, end, i);
	return;
}

/*
 * Used by reporting mode 4 (see _Biostrings_reset_viewsbuf() below).
 */
static void merge_with_last_view(int start, int end)
{
	int end1;

	if (viewsbuf_count == 0) {
		insert_view_at(start, end, viewsbuf_count);
		return;
	}
	end1 = viewsbuf_end[viewsbuf_count - 1] + 1;
	if (start > end1) {
		insert_view_at(start, end, viewsbuf_count);
		return;
	}
	if (end >= end1)
		viewsbuf_end[viewsbuf_count - 1] = end;
	return;
}

/*
 * Used by reporting mode 5 (see _Biostrings_reset_viewsbuf() below).
 */
static void merge_view(int start, int end)
{
	int start1, end1, i, j;

	start1 = start - 1;
	end1 = end + 1;
	j = viewsbuf_count - 1;
	while (0 <= j && end1 < viewsbuf_start[j]) j--;
	i = j;
	while (0 <= i && start1 <= viewsbuf_end[i]) i--;
	if (i == j) {
		insert_view_at(start, end, i + 1);
		return;
	}
	i++;
	if (viewsbuf_start[i] < start)
		start = viewsbuf_start[i];
	if (end < viewsbuf_end[j])
		end = viewsbuf_end[j];
	viewsbuf_start[i] = start;
	viewsbuf_end[i] = end;
	if (i == j)
		return;
	i++;
	j++;
	while (j < viewsbuf_count) {
		viewsbuf_start[i] = viewsbuf_start[j];
		viewsbuf_end[i] = viewsbuf_end[j];
		i++;
		j++;
	}
	viewsbuf_count = i;
	return;
}

/*
 * Reset views buffer
 *
 * 6 valid "reporting modes":
 *
 *   0: Views must be reported thru _Biostrings_append_view(start, end, desc).
 *      They are not reordered, nor checked for duplicated.
 *
 *   1: Matches must be reported thru _Biostrings_report_match(Lpos, Rpos).
 *      They are counted only (no need to allocate memory to store them).
 *
 *   2: Matches must be reported thru _Biostrings_report_match(Lpos, Rpos).
 *      They are not reordered, nor checked for duplicated.
 *
 *   3: Matches must be reported thru _Biostrings_report_match(Lpos, Rpos).
 *      They are reordered and the duplicated are ignored.
 *
 *   4: Matches must be reported thru _Biostrings_report_match(Lpos, Rpos)
 *      "in ascending Lpos order" i.e. for each new call to
 *      _Biostrings_report_match(), the value passed to the 'Lpos' argument
 *      must greater or equal than for the previous call.
 *      Each new match is merged with the previous match if overlapping or
 *      adjacent so the current set of matches is kept normalized.
 *
 *   5: Matches must be reported thru _Biostrings_report_match(Lpos, Rpos).
 *      They are reordered and merged when overlapping or adjacent
 *      (normalization). Hence a call to _Biostrings_report_match() can
 *      actually reduce the current set of matches if several matches are
 *      replaced by a single one.
 *      IMPORTANT: If the matches are reported "in ascending Lpos order",
 *      mode 4 and 5 will be semantically equivalent: both will normalize
 *      the set of matches but 4 will be slightly faster (20-30%) and
 *      should be preferred. Use mode 5 only if there is no guarantee that
 *      the matches will be found (and reported) "in ascending Lpos order".
 *      Note that the reordering is performed every time a new match is
 *      reported by using a simple bubble sort approach. This is clearly
 *      not optimal since the bubble sort approach is not efficient when
 *      the input needs a lot of reordering. However the bubble sort approach
 *      is simple, has minimal memory footprint and its performance is not
 *      too bad when the number of matches is small (< 10^5) or doesn't need
 *      too much reodering.
 */
void _Biostrings_reset_viewsbuf(int reporting_mode)
{
	viewsbuf_reporting_mode = reporting_mode;
	/* No memory leak here, because we use transient storage allocation */
	viewsbuf_start = viewsbuf_end = NULL;
	viewsbuf_desc = NULL;
	viewsbuf_maxcount = viewsbuf_count = 0;
	return;
}

/*
 * Can only be used in reporting mode 0.
 * Return the new number of views.
 */
int _Biostrings_append_view(int start, int end, const char *desc)
{
	int i;
	size_t desc_size;

	if (viewsbuf_reporting_mode != 0)
		error("_Biostrings_append_view(): viewsbuf_reporting_mode != 0");
	i = new_view();
	viewsbuf_start[i] = start;
	viewsbuf_end[i] = end;
	if (desc != NULL) {
		desc_size = strlen(desc) + 1; /* + 1 for the terminating '\0' character */
		viewsbuf_desc[i] = Salloc((long) desc_size, char);
		memcpy(viewsbuf_desc[i], desc, desc_size);
	} else {
		viewsbuf_desc[i] = NULL;
	}
	return viewsbuf_count;
}

/*
 * Can only be used in reporting mode != 0.
 */
int _Biostrings_report_match(int Lpos, int Rpos)
{
	int start, end;

	start = ++Lpos;
	end = ++Rpos;
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _Biostrings_report_match(): ");
		Rprintf("match found at start=%d end=%d --> ", start, end);
	}
#endif
	switch (viewsbuf_reporting_mode) {
		case 0:
			error("_Biostrings_report_match(): viewsbuf_reporting_mode == 0");
		break;
		case 1:
			viewsbuf_count++;
		break;
		case 2:
			insert_view_at(start, end, viewsbuf_count);
		break;
		case 3:
			insert_view_if_new(start, end);
		break;
		case 4:
			merge_with_last_view(start, end);
		break;
		case 5:
			merge_view(start, end);
		break;
	}
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("viewsbuf_count=%d\n", viewsbuf_count);
	}
#endif
	return viewsbuf_count;
}

SEXP _Biostrings_viewsbuf_count_asINTEGER()
{
	SEXP ans;

	PROTECT(ans = NEW_INTEGER(1));
	INTEGER(ans)[0] = viewsbuf_count;
	UNPROTECT(1);
	return ans;
}

SEXP _Biostrings_viewsbuf_start_asINTEGER()
{
	SEXP ans;

	if (viewsbuf_reporting_mode == 1)
		error("_Biostrings_viewsbuf_start_asINTEGER(): viewsbuf_reporting_mode == 1");
	PROTECT(ans = NEW_INTEGER(viewsbuf_count));
	memcpy(INTEGER(ans), viewsbuf_start, sizeof(int) * viewsbuf_count);
	UNPROTECT(1);
	return ans;
}

SEXP _Biostrings_viewsbuf_end_asINTEGER()
{
	SEXP ans;

	if (viewsbuf_reporting_mode == 1)
		error("_Biostrings_viewsbuf_end_asINTEGER(): viewsbuf_reporting_mode == 1");
	PROTECT(ans = NEW_INTEGER(viewsbuf_count));
	memcpy(INTEGER(ans), viewsbuf_end, sizeof(int) * viewsbuf_count);
	UNPROTECT(1);
	return ans;
}

SEXP _Biostrings_viewsbuf_desc_asCHARACTER()
{
	SEXP ans;
	int i;

	if (viewsbuf_reporting_mode != 0)
		error("_Biostrings_viewsbuf_desc_asCHARACTER(): viewsbuf_reporting_mode != 0");
	PROTECT(ans = NEW_CHARACTER(viewsbuf_count));
	for (i = 0; i < viewsbuf_count; i++)
		SET_STRING_ELT(ans, i, viewsbuf_desc[i] == NULL ? NA_STRING : mkChar(viewsbuf_desc[i]));
	UNPROTECT(1);
	return ans;
}

SEXP _Biostrings_viewsbuf_asLIST()
{
	SEXP ans, ans_names, ans_elt;

	if (viewsbuf_reporting_mode == 1)
		error("_Biostrings_viewsbuf_asLIST(): viewsbuf_reporting_mode == 1");
	PROTECT(ans = NEW_LIST(2));
	/* set the names */
	PROTECT(ans_names = NEW_CHARACTER(2));
	SET_STRING_ELT(ans_names, 0, mkChar("start"));
	SET_STRING_ELT(ans_names, 1, mkChar("end"));
	SET_NAMES(ans, ans_names);
	UNPROTECT(1);
	/* set the "start" element */
	PROTECT(ans_elt = _Biostrings_viewsbuf_start_asINTEGER());
	SET_ELEMENT(ans, 0, ans_elt);
	UNPROTECT(1);
	/* set the "end" element */
	PROTECT(ans_elt = _Biostrings_viewsbuf_end_asINTEGER());
	SET_ELEMENT(ans, 1, ans_elt);
	UNPROTECT(1);
	/* ans is ready */
	UNPROTECT(1);
	return ans;
}

SEXP _reported_matches_asSEXP()
{
	if (viewsbuf_reporting_mode == 1)
		return _Biostrings_viewsbuf_count_asINTEGER();
	if (viewsbuf_reporting_mode == 2)
		return _Biostrings_viewsbuf_start_asINTEGER();
	error("_reported_matches_asSEXP(): unsupported reporting mode %d", viewsbuf_reporting_mode);
	return R_NilValue;
}

