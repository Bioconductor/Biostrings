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

static int *views_startbuf, *views_endbuf;
static char **views_descbuf;
static int views_bufsize, views_count; /* views_bufsize >= views_count */
static int count_only_mode; /* 1 for yes, 0 for no */

static int new_view()
{
	long new_size;

	if (views_count >= views_bufsize) {
		/* Buffer is full */
		if (views_bufsize == 0)
			new_size = 1024;
		else
			new_size = 2 * views_bufsize;
		views_startbuf = Srealloc((char *) views_startbuf, new_size,
						(long) views_bufsize, int);
		views_endbuf = Srealloc((char *) views_endbuf, new_size,
						(long) views_bufsize, int);
		views_descbuf = Srealloc((char *) views_descbuf, new_size,
						(long) views_bufsize, char *);
		views_bufsize = new_size;
	}
	return views_count++;
}

static void set_view(int i, int start, int end, const char *desc)
{
	size_t desc_size;

	views_startbuf[i] = start;
	views_endbuf[i] = end;
	if (desc == NULL) {
		views_descbuf[i] = NULL;
	} else {
		desc_size = strlen(desc) + 1; /* + 1 for the terminating '\0' character */
		views_descbuf[i] = Salloc((long) desc_size, char);
		memcpy(views_descbuf[i], desc, desc_size);
	}
	return;
}

/* Reset views buffer */
void _Biostrings_reset_views_buffer(int count_only)
{
	/* No memory leak here, because we use transient storage allocation */
	views_startbuf = views_endbuf = NULL;
	views_descbuf = NULL;
	views_bufsize = views_count = 0;
	count_only_mode = count_only;
	return;
}

/* Return the new number of views */
int _Biostrings_report_view(int start, int end, const char *desc)
{
	if (count_only_mode)
		views_count++;
	else
		set_view(new_view(), start, end, desc);
	return views_count;
}

/* Return 0 if the match is ignored (because it's a duplicated) or 1 if it's counted */
int _Biostrings_report_match(int Lpos, int Rpos)
{
	int start, end, i, j1, j2;

	start = ++Lpos;
	end = ++Rpos;
	i = views_count - 1;
	while (0 <= i && start < views_startbuf[i]) i--;
	while (0 <= i && start == views_startbuf[i] && end < views_endbuf[i]) i--;
	if (0 <= i && start == views_startbuf[i] && end == views_endbuf[i])
		return 0;
	if (count_only_mode) {
		views_count++;
	} else {
		j2 = new_view();
		j1 = j2 - 1;
		while (j1 > i) {
			views_startbuf[j2] = views_startbuf[j1];
			views_endbuf[j2] = views_endbuf[j1];
			views_descbuf[j2] = views_descbuf[j1];
			j1--;
			j2--;
		}
		set_view(j2, start, end, NULL);
	}
	return 1;
}

SEXP _Biostrings_get_views_count_INTEGER()
{
	SEXP ans;

	PROTECT(ans = NEW_INTEGER(1));
	INTEGER(ans)[0] = views_count;
	UNPROTECT(1);
	return ans;
}

SEXP _Biostrings_get_views_start_INTEGER()
{
	SEXP ans;

	if (count_only_mode)
		error("Biostrings internals: _Biostrings_get_views_start_INTEGER() called while in \"count only\" mode");
	PROTECT(ans = NEW_INTEGER(views_count));
	memcpy(INTEGER(ans), views_startbuf, sizeof(int) * views_count);
	UNPROTECT(1);
	return ans;
}

SEXP _Biostrings_get_views_end_INTEGER()
{
	SEXP ans;

	if (count_only_mode)
		error("Biostrings internals: _Biostrings_get_views_end_INTEGER() called while in \"count only\" mode");
	PROTECT(ans = NEW_INTEGER(views_count));
	memcpy(INTEGER(ans), views_endbuf, sizeof(int) * views_count);
	UNPROTECT(1);
	return ans;
}

SEXP _Biostrings_get_views_desc_CHARACTER()
{
	SEXP ans;
	int i;

	if (count_only_mode)
		error("Biostrings internals: _Biostrings_get_views_desc_CHARACTER() called while in \"count only\" mode");
	PROTECT(ans = NEW_CHARACTER(views_count));
	for (i = 0; i < views_count; i++)
		SET_STRING_ELT(ans, i, mkChar(views_descbuf[i]));
	UNPROTECT(1);
	return ans;
}

SEXP _Biostrings_get_views_LIST()
{
	SEXP ans, ans_names, ans_elt;

	if (count_only_mode)
		error("Biostrings internals: _Biostrings_get_views_LIST() called while in \"count only\" mode");
	PROTECT(ans = NEW_LIST(2));
	/* set the names */
	PROTECT(ans_names = NEW_CHARACTER(2));
	SET_STRING_ELT(ans_names, 0, mkChar("start"));
	SET_STRING_ELT(ans_names, 1, mkChar("end"));
	SET_NAMES(ans, ans_names);
	UNPROTECT(1);
	/* set the "start" element */
	PROTECT(ans_elt = _Biostrings_get_views_start_INTEGER());
	SET_ELEMENT(ans, 0, ans_elt);
	UNPROTECT(1);
	/* set the "end" element */
	PROTECT(ans_elt = _Biostrings_get_views_end_INTEGER());
	SET_ELEMENT(ans, 1, ans_elt);
	UNPROTECT(1);
	/* ans is ready */
	UNPROTECT(1);
	return ans;
}

