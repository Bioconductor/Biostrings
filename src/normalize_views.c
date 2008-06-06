/****************************************************************************
 *                        Normalizing a set of views                        *
 *                           Author: Herve Pages                            *
 *                           -------------------                            *
 *                                                                          *
 * See R/mask.R for a quick intro to views normalization.                   *
 *                                                                          *
 ****************************************************************************/
#include "Biostrings.h"

#include <stdio.h>


/*
 * Assumes that 'start[i] <= end[i]' for all 0 <= i < nviews
 * and that 'start' is sorted in increasing order.
 *
 * Some benchmarks:
 *   > set.seed(23)
 *   > start <- sort(as.integer(runif(10^7, min=1, max=10^8)))
 *   > sum(duplicated(start))
 *   [1] 483908
 *   > mean(diff(start))
 *   [1] 10
 *   > end <- start + abs(as.integer(rnorm(10^7, mean=8, sd=4)))
 *   > sum(start == end)
 *   [1] 278805
 *   > mean(end - start)
 *   [1] 7.571926
 *   > max(end - start)
 *   [1] 29
 *   > system.time(views0 <- .Call("Biostrings_normalize_views", start, end, PACKAGE="Biostrings"))
 *      user  system elapsed
 *     0.216   0.000   0.217
 *   > length(views0$start)
 *   [1] 4038960
 *   > as.data.frame(views0)[1:7, ]
 *     start end
 *   1     3   7
 *   2    12  26
 *   3    60  70
 *   4    91 114
 *   5   133 147
 *   6   177 194
 *   7   200 216
 *
 */
static void normalize_orderedbystartviews(const int *start, const int *end, int nviews)
{
        int i;

	_init_match_reporting(4);
	for (i = 0; i < nviews; i++)
		_Biostrings_report_match(*(start++) - 1, *(end++) - 1);
        return;
}

/*
 * Only assumes that 'start[i] <= end[i]' for all 0 <= i < nviews. Nothing
 * else. Because when used in reporting mode 5, _Biostrings_report_match()
 * takes care of keeping the views buffer normalized.
 * Only 30% slower than normalize_orderedbystartviews() when the views are
 * already ordered by start (i.e. 'start' is sorted in increasing order).
 */
static void normalize_views(const int *start, const int *end, int nviews)
{
	int i;

	_init_match_reporting(5);
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
	//normalize_orderedbystartviews(INTEGER(start), INTEGER(end), LENGTH(start));
	normalize_views(INTEGER(start), INTEGER(end), LENGTH(start));
	return _reported_matches_asLIST();
}

