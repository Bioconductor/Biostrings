#include "Biostrings.h"

#include <stdio.h>


/*
 * Assumes that 'start[i] <= end[i]' for all 0 <= i < nviews
 * and that 'start' is sorted in increasing order.
 */
static void normalize_views(const int *start, const int *end, int nviews)
{
	int i;

	for (i = 0; i < nviews; i++)
		_Biostrings_report_match(*(start++) - 1, *(end++) - 1);
	return;
}

/*
 * 'start' and 'end': the set of views i.e. 2 NA-free integer vectors of the
 *                    same length and such that 'start <= end'
 */
SEXP Biostrings_normalize_views(SEXP start, SEXP end)
{
	_Biostrings_reset_viewsbuf(4);
	normalize_views(INTEGER(start), INTEGER(end), LENGTH(start));
	return _Biostrings_viewsbuf_asLIST();
}

