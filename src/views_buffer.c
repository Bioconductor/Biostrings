#include "Biostrings.h"
#include <S.h> /* for Srealloc() */


/* ==========================================================================
 * Helper functions used for storing views in a temporary buffer like the
 * matches found by the matching algos.
 * --------------------------------------------------------------------------
 */

static int *views_startbuf, *views_endbuf;
static char **views_descbuf;
static int views_bufsize, views_count; /* views_bufsize >= views_count */

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
void _Biostrings_reset_views_buffer()
{
	/* No memory leak here, because we use transient storage allocation */
	views_startbuf = views_endbuf = NULL;
	views_descbuf = NULL;
	views_bufsize = views_count = 0;
	return;
}

int *_Biostrings_get_views_start()
{
	return views_startbuf;
}

int *_Biostrings_get_views_end()
{
	return views_endbuf;
}

char **_Biostrings_get_views_desc()
{
	return views_descbuf;
}

/* Return the new number of views */
int _Biostrings_report_view(int start, int end, const char *desc)
{
	set_view(new_view(), start, end, desc);
	return views_count;
}

/* Return 0 if the match is ignored (because it's a duplicated) or 1 if it's inserted */
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
	return 1;
}

SEXP _Biostrings_get_views_LIST()
{
	SEXP ans, ans_names, ans_elt;

	PROTECT(ans = NEW_LIST(2));
	/* set the names */
	PROTECT(ans_names = allocVector(STRSXP, 2));
	SET_STRING_ELT(ans_names, 0, mkChar("start"));
	SET_STRING_ELT(ans_names, 1, mkChar("end"));
	SET_NAMES(ans, ans_names);
	UNPROTECT(1);
	/* set the "start" element */
	PROTECT(ans_elt = allocVector(INTSXP, views_count));
	memcpy(INTEGER(ans_elt), _Biostrings_get_views_start(), sizeof(int) * views_count);
	SET_ELEMENT(ans, 0, ans_elt);
	UNPROTECT(1);
	/* set the "end" element */
	PROTECT(ans_elt = allocVector(INTSXP, views_count));
	memcpy(INTEGER(ans_elt), _Biostrings_get_views_end(), sizeof(int) * views_count);
	SET_ELEMENT(ans, 1, ans_elt);
	UNPROTECT(1);
	/* ans is ready */
	UNPROTECT(1);
	return ans;
}

