#include "Biostrings.h"

#include <stdio.h>


/*
 * Assumes that 'start[i] <= end[i]' for all 0 <= i < nviews
 * and that 'start' is sorted in increasing order.
 */
static int normalize_views(const int *start, const int *end, int nviews)
{
	int start0, end0, start1, end1, i;

	if (nviews == 0)
		return 0;
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
		_Biostrings_report_view(start0, end0, "");
		start0 = start1;
		end0 = end1;
	}
	return _Biostrings_report_view(start0, end0, "");
}

/*
 * 'start' and 'end: integer vectors of same length and such
 * that 'start <= end' and 'start' sorted in increasing order.
 */
SEXP Biostrings_normalize_views(SEXP start, SEXP end)
{
	int count;
	SEXP ans, ans_names, ans_elt;

	_Biostrings_reset_views_buffer();
	count = normalize_views(INTEGER(start), INTEGER(end), LENGTH(start));

	PROTECT(ans = NEW_LIST(2));
	/* set the names */
	PROTECT(ans_names = allocVector(STRSXP, 2));
	SET_STRING_ELT(ans_names, 0, mkChar("start"));
	SET_STRING_ELT(ans_names, 1, mkChar("end"));
	SET_NAMES(ans, ans_names);
	UNPROTECT(1);
	/* set the "start" element */
	PROTECT(ans_elt = allocVector(INTSXP, count));
	memcpy(INTEGER(ans_elt), _Biostrings_get_views_start(), sizeof(int) * count);
	SET_ELEMENT(ans, 0, ans_elt);
	UNPROTECT(1);
	/* set the "end" element */
	PROTECT(ans_elt = allocVector(INTSXP, count));
	memcpy(INTEGER(ans_elt), _Biostrings_get_views_end(), sizeof(int) * count);
	SET_ELEMENT(ans, 1, ans_elt);
	UNPROTECT(1);
	/* ans is ready */
	UNPROTECT(1);
	return ans;
}

