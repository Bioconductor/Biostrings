#include "Biostrings.h"

#include <stdio.h>


/*
 * Assumes that 'start[i] <= end[i]' for all 0 <= i < nviews
 * and that 'start' is sorted in increasing order.
 */
static void normalize_views(const int *start, const int *end, int nviews)
{
	int start0, end0, start1, end1, i;

	if (nviews == 0)
		return;
	start0 = start[0];
	end0 = end[0];
	for (i = 1; i < nviews; i++) {
		start1 = start[i];
		end1 = end[i];
		if (start1 <= end0 + 1) {
			if (end1 > end0)
				end0 = end1;
			continue;
		}
		_Biostrings_report_view(start0, end0, NULL);
		start0 = start1;
		end0 = end1;
	}
	_Biostrings_report_view(start0, end0, NULL);
	return;
}

/*
 * 'start' and 'end: integer vectors of same length and such
 * that 'start <= end' and 'start' sorted in increasing order.
 */
SEXP Biostrings_normalize_views(SEXP start, SEXP end)
{
	_Biostrings_reset_viewsbuf(0);
	normalize_views(INTEGER(start), INTEGER(end), LENGTH(start));
	return _Biostrings_viewsbuf_asLIST();
}

