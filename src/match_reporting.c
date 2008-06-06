/*
 * Match (and view) reporting facilities.
 *
 * The "match buffer" is used for temporarily storing a set of matches (or
 * views) found in (or defined on) a sequence. E.g. it is used by the string
 * searching functions like find_palindromes() or the pattern matching
 * functions.
 * Except for debug_match_reporting(), the functions defined in this
 * file are NOT .Call entry points so THEY DON'T NEED TO BE REGISTERED
 * in R_init_Biostrings.c.
 */
#include "Biostrings.h"
#include <S.h> /* for Srealloc() */


static int debug = 0;

SEXP debug_match_reporting()
{
#ifdef DEBUG_BIOSTRINGS
        debug = !debug;
        Rprintf("Debug mode turned %s in 'match_reporting.c'\n",
		debug ? "on" : "off");
#else
        Rprintf("Debug mode not available in 'match_reporting.c'\n");
#endif
        return R_NilValue;
}

static int matchbuf_mrmode;

static int *matchbuf_start, *matchbuf_end;
static char **matchbuf_name;
static int matchbuf_maxcount, matchbuf_count; /* matchbuf_maxcount >= matchbuf_count */
static int match_shift;

static int new_match()
{
	long new_maxcount;

	if (matchbuf_count >= matchbuf_maxcount) {
		/* Buffer is full */
		if (matchbuf_maxcount == 0)
			new_maxcount = 1024;
		else
			new_maxcount = 2 * matchbuf_maxcount;
		matchbuf_start = Srealloc((char *) matchbuf_start, new_maxcount,
						(long) matchbuf_maxcount, int);
		matchbuf_end = Srealloc((char *) matchbuf_end, new_maxcount,
						(long) matchbuf_maxcount, int);
		if (matchbuf_mrmode == 0)
			matchbuf_name = Srealloc((char *) matchbuf_name, new_maxcount,
						(long) matchbuf_maxcount, char *);
		matchbuf_maxcount = new_maxcount;
	}
	return matchbuf_count++;
}

/*
 * Must be used in mrmode >= 2 (see _init_match_reporting() below).
 */
static void insert_match_at(int start, int end, int insert_at)
{
	int i, j;

	j = new_match(); // matchbuf_count - 1
	i = j - 1;
	while (insert_at <= i) {
		matchbuf_start[j] = matchbuf_start[i];
		matchbuf_end[j] = matchbuf_end[i];
		i--;
		j--;
	}
	matchbuf_start[j] = start;
	matchbuf_end[j] = end;
	return;
}

/*
 * Used by mrmode 3 (see _init_match_reporting() below).
 */
static void insert_match_if_new(int start, int end)
{
	int i;

	i = matchbuf_count - 1;
	while (0 <= i && start < matchbuf_start[i]) i--;
	while (0 <= i && start == matchbuf_start[i] && end < matchbuf_end[i]) i--;
	if (0 <= i && start == matchbuf_start[i] && end == matchbuf_end[i])
		return;
	i++;
	insert_match_at(start, end, i);
	return;
}

/*
 * Used by mrmode 4 (see _init_match_reporting() below).
 */
static void merge_with_last_match(int start, int end)
{
	int end1;

	if (matchbuf_count == 0) {
		insert_match_at(start, end, matchbuf_count);
		return;
	}
	end1 = matchbuf_end[matchbuf_count - 1] + 1;
	if (start > end1) {
		insert_match_at(start, end, matchbuf_count);
		return;
	}
	if (end >= end1)
		matchbuf_end[matchbuf_count - 1] = end;
	return;
}

/*
 * Used by mrmode 5 (see _init_match_reporting() below).
 */
static void merge_match(int start, int end)
{
	int start1, end1, i, j;

	start1 = start - 1;
	end1 = end + 1;
	j = matchbuf_count - 1;
	while (0 <= j && end1 < matchbuf_start[j]) j--;
	i = j;
	while (0 <= i && start1 <= matchbuf_end[i]) i--;
	if (i == j) {
		insert_match_at(start, end, i + 1);
		return;
	}
	i++;
	if (matchbuf_start[i] < start)
		start = matchbuf_start[i];
	if (end < matchbuf_end[j])
		end = matchbuf_end[j];
	matchbuf_start[i] = start;
	matchbuf_end[i] = end;
	if (i == j)
		return;
	i++;
	j++;
	while (j < matchbuf_count) {
		matchbuf_start[i] = matchbuf_start[j];
		matchbuf_end[i] = matchbuf_end[j];
		i++;
		j++;
	}
	matchbuf_count = i;
	return;
}

/*
 * _init_match_reporting() suports 6 valid "match reporting modes":
 *
 *   0: Views must be reported thru _report_view(start, end, name).
 *      They are not reordered, nor checked for duplicated.
 *
 *   1: Matches must be reported thru _report_match(start, end).
 *      They are counted only (no need to allocate memory to store them).
 *
 *   2: Matches must be reported thru _report_match(start, end).
 *      They are not reordered, nor checked for duplicated.
 *
 *   3: Matches must be reported thru _report_match(start, end).
 *      They are reordered and the duplicated are ignored.
 *
 *   4: Matches must be reported thru _report_match(start, end)
 *      "in ascending start order" i.e. for each new call to
 *      _report_match(), the value passed to the 'start' argument
 *      must be greater or equal than for the previous call.
 *      Each new match is merged with the previous match if overlapping or
 *      adjacent so the current set of matches is kept normalized.
 *
 *   5: Matches must be reported thru _report_match(start, end).
 *      They are reordered and merged when overlapping or adjacent
 *      (normalization). Hence a call to _report_match() can
 *      actually reduce the current set of matches if several matches are
 *      replaced by a single one.
 *      IMPORTANT: If the matches are reported "in ascending start order",
 *      mrmode 4 and 5 will be semantically equivalent: both will normalize
 *      the set of matches but 4 will be slightly faster (20-30%) and
 *      should be preferred. Use mrmode 5 only if there is no guarantee that
 *      the matches will be found (and reported) "in ascending start order".
 *      Note that the reordering is performed every time a new match is
 *      reported by using a simple bubble sort approach. This is clearly
 *      not optimal since the bubble sort approach is not efficient when
 *      the input needs a lot of reordering. However the bubble sort approach
 *      is simple, has minimal memory footprint and its performance is not
 *      too bad when the number of matches is small (< 10^5) or doesn't need
 *      too much reodering.
 */

void _init_match_reporting(int mrmode)
{
	matchbuf_mrmode = mrmode;
	/* No memory leak here, because we use transient storage allocation */
	matchbuf_start = matchbuf_end = NULL;
	matchbuf_name = NULL;
	matchbuf_maxcount = matchbuf_count = 0;
	match_shift = 0;
	return;
}

void _drop_current_matches()
{
	matchbuf_count = 0;
}

void _set_match_shift(int shift)
{
	match_shift = shift;
}

/*
 * Can only be used in mrmode 0.
 * Return the new number of views.
 */
int _report_view(int start, int end, const char *name)
{
	int i;
	size_t name_size;

	if (matchbuf_mrmode != 0)
		error("_report_view(): matchbuf_mrmode != 0");
	i = new_match();
	matchbuf_start[i] = start;
	matchbuf_end[i] = end;
	if (name != NULL) {
		name_size = strlen(name) + 1; /* + 1 for the terminating '\0' character */
		matchbuf_name[i] = Salloc((long) name_size, char);
		memcpy(matchbuf_name[i], name, name_size);
	} else {
		matchbuf_name[i] = NULL;
	}
	return matchbuf_count;
}

/*
 * Can only be used in mrmode != 0.
 */
int _report_match(int start, int end)
{
	start += match_shift;
	end += match_shift;
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _report_match(): ");
		Rprintf("match found at start=%d end=%d --> ", start, end);
	}
#endif
	switch (matchbuf_mrmode) {
		case 0:
			error("_report_match(): matchbuf_mrmode == 0");
		break;
		case 1:
			matchbuf_count++;
		break;
		case 2:
			insert_match_at(start, end, matchbuf_count);
		break;
		case 3:
			insert_match_if_new(start, end);
		break;
		case 4:
			merge_with_last_match(start, end);
		break;
		case 5:
			merge_match(start, end);
		break;
	}
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("matchbuf_count=%d\n", matchbuf_count);
	}
#endif
	return matchbuf_count;
}

SEXP _reported_match_count_asINTEGER()
{
	SEXP ans;

	PROTECT(ans = NEW_INTEGER(1));
	INTEGER(ans)[0] = matchbuf_count;
	UNPROTECT(1);
	return ans;
}

SEXP _reported_match_starts_asINTEGER()
{
	SEXP ans;

	if (matchbuf_mrmode == 1)
		error("_reported_match_starts_asINTEGER(): matchbuf_mrmode == 1");
	PROTECT(ans = NEW_INTEGER(matchbuf_count));
	memcpy(INTEGER(ans), matchbuf_start, sizeof(int) * matchbuf_count);
	UNPROTECT(1);
	return ans;
}

SEXP _reported_match_ends_asINTEGER()
{
	SEXP ans;

	if (matchbuf_mrmode == 1)
		error("_reported_match_ends_asINTEGER(): matchbuf_mrmode == 1");
	PROTECT(ans = NEW_INTEGER(matchbuf_count));
	memcpy(INTEGER(ans), matchbuf_end, sizeof(int) * matchbuf_count);
	UNPROTECT(1);
	return ans;
}

SEXP _reported_view_names_asCHARACTER()
{
	SEXP ans;
	int i;

	if (matchbuf_mrmode != 0)
		error("_reported_view_names_asCHARACTER(): matchbuf_mrmode != 0");
	PROTECT(ans = NEW_CHARACTER(matchbuf_count));
	for (i = 0; i < matchbuf_count; i++)
		SET_STRING_ELT(ans, i, matchbuf_name[i] == NULL ? NA_STRING : mkChar(matchbuf_name[i]));
	UNPROTECT(1);
	return ans;
}

SEXP _reported_matches_asLIST()
{
	SEXP ans, ans_names, ans_elt;

	if (matchbuf_mrmode == 1)
		error("_reported_matches_asLIST(): matchbuf_mrmode == 1");
	PROTECT(ans = NEW_LIST(2));
	/* set the names */
	PROTECT(ans_names = NEW_CHARACTER(2));
	SET_STRING_ELT(ans_names, 0, mkChar("start"));
	SET_STRING_ELT(ans_names, 1, mkChar("end"));
	SET_NAMES(ans, ans_names);
	UNPROTECT(1);
	/* set the "start" element */
	PROTECT(ans_elt = _reported_match_starts_asINTEGER());
	SET_ELEMENT(ans, 0, ans_elt);
	UNPROTECT(1);
	/* set the "end" element */
	PROTECT(ans_elt = _reported_match_ends_asINTEGER());
	SET_ELEMENT(ans, 1, ans_elt);
	UNPROTECT(1);
	/* ans is ready */
	UNPROTECT(1);
	return ans;
}

SEXP _reported_matches_asSEXP()
{
	if (matchbuf_mrmode == 1)
		return _reported_match_count_asINTEGER();
	if (matchbuf_mrmode == 2)
		return _reported_match_starts_asINTEGER();
	error("_reported_matches_asSEXP(): invalid mrmode %d", matchbuf_mrmode);
	return R_NilValue;
}

