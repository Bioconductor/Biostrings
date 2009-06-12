/* Driver routine to call RCONT2 from R, B times.
   Calculates the log-likelihood ratio for each generated table.
   Largely a cut'n'paste translation of chisqsim()
   Append this file to []/R-x.x.x/src/library/ctest/chisqsim.c
   Pete Hurd - Sept 29 2001
   */

void
gtestsim(int *nrow, int *ncol, int *nrowt, int *ncolt, int *n,
	 int *b, double *expected, int *observed, double *fact,
	 int *jwork, double *results)
{
    /* Local variables */
    int i, j, iter;
    double g, e, o, x;

    /* Calculate log-factorials */
    x = 0.;
    fact[0] = 0.;
    for (i = 1; i <= *n; ++i) {
	x += log((double) i);
	fact[i] = x;
    }

    GetRNGstate();

    for (iter = 0; iter < *b; ++iter) {
	rcont2(nrow, ncol, nrowt, ncolt, n, fact, jwork, observed);
	/* Calculate G value from the random table: */
	g = 0.;
	for (i = 0; i < *nrow; ++i) {
	    for (j = 0; j < *ncol; ++j) {
		e = expected[i + j * *nrow];
		o = observed[i + j * *nrow];
		if (o!=0) g += o * log(o / e);
	    }
	}
	g = 2 * g;
	results[iter] = g;
    }

    PutRNGstate();

    return;
}
