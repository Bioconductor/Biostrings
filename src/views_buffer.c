/*
 * Manipulation of the "views buffer".
 *
 * The "views buffer" is used for temporarily storing a set of views or
 * matches found in a sequence (e.g. the matches found by the matching algos).
 * Except for Biostrings_debug_views_buffer(), the functions defined in this
 * file are NOT .Call methods (but they are used by .Call methods defined in
 * other .c files) so THEY DON'T NEED TO BE REGISTERED (in R_init_Biostrings.c).
 * They are prefixed with "_Biostrings_" to minimize the risk of a clash
 * with symbols defined elsewhere (e.g. in libc).
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

/* 4 valid modes:
 *   0: views must be reported thru _Biostrings_report_view(start, end, desc),
 *      they are not reordered, not checked for duplicated
 *   1: matches must be reported thru _Biostrings_report_match(Lpos, Rpos),
 *      they are counted only
 *   2: matches must be reported thru _Biostrings_report_match(Lpos, Rpos),
 *      they are not reordered, not checked for duplicated
 *   3: matches must be reported thru _Biostrings_report_match(Lpos, Rpos),
 *      they are reordered and the duplicated are ignored
 */
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

static void set_view(int i, int start, int end, const char *desc)
{
	size_t desc_size;

	viewsbuf_start[i] = start;
	viewsbuf_end[i] = end;
	if (viewsbuf_reporting_mode != 0)
		return;
	if (desc != NULL) {
		desc_size = strlen(desc) + 1; /* + 1 for the terminating '\0' character */
		viewsbuf_desc[i] = Salloc((long) desc_size, char);
		memcpy(viewsbuf_desc[i], desc, desc_size);
	} else {
		viewsbuf_desc[i] = NULL;
	}
	return;
}

/* Reset views buffer */
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
int _Biostrings_report_view(int start, int end, const char *desc)
{
	if (viewsbuf_reporting_mode != 0)
		error("_Biostrings_report_view(): viewsbuf_reporting_mode != 0");
	set_view(new_view(), start, end, desc);
	return viewsbuf_count;
}

/*
 * Can only be used in reporting mode 1, 2 and 3.
 * Return 1 if the match is counted, 0 if it's ignored.
 */
int _Biostrings_report_match(int Lpos, int Rpos)
{
	int start, end, i, j1, j2;

	if (viewsbuf_reporting_mode == 0)
		error("_Biostrings_report_match(): viewsbuf_reporting_mode == 0");
	start = ++Lpos;
	end = ++Rpos;
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _Biostrings_report_match(): ");
		Rprintf("match found at start=%d end=%d --> ", start, end);
	}
#endif
	if (viewsbuf_reporting_mode == 1) {
		viewsbuf_count++;
#ifdef DEBUG_BIOSTRINGS
		if (debug) Rprintf("COUNTED (viewsbuf_count=%d)\n", viewsbuf_count);
#endif
		return 1;
	}
	if (viewsbuf_reporting_mode == 2) {
		set_view(new_view(), start, end, NULL);
#ifdef DEBUG_BIOSTRINGS
		if (debug) Rprintf("COUNTED (viewsbuf_count=%d)\n", viewsbuf_count);
#endif
		return 1;
	}
	i = viewsbuf_count - 1;
	while (0 <= i && start < viewsbuf_start[i]) i--;
	while (0 <= i && start == viewsbuf_start[i] && end < viewsbuf_end[i]) i--;
	if (0 <= i && start == viewsbuf_start[i] && end == viewsbuf_end[i]) {
#ifdef DEBUG_BIOSTRINGS
		if (debug) Rprintf("IGNORED (duplicated)\n");
#endif
		return 0;
	}
	j2 = new_view();
	j1 = j2 - 1;
	while (j1 > i) {
		viewsbuf_start[j2] = viewsbuf_start[j1];
		viewsbuf_end[j2] = viewsbuf_end[j1];
		j1--;
		j2--;
	}
	set_view(j2, start, end, NULL);
#ifdef DEBUG_BIOSTRINGS
	if (debug) Rprintf("COUNTED (viewsbuf_count=%d)\n", viewsbuf_count);
#endif
	return 1;
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

