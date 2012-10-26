#include <Rmath.h>
#include <R_ext/Random.h>
#include <R_ext/Error.h>
#include <R_ext/Utils.h>

/* rcont2 function copy'n'pasted from R/src/library/stats/src/rcont.c */
static void
rcont2(int *nrow, int *ncol,
       /* vectors of row and column totals, and their sum ntotal: */
       int *nrowt, int *ncolt, int *ntotal,
       double *fact, int *jwork, int *matrix)
{
    int j, l, m, ia, ib, ic, jc, id, ie, ii, nll, nlm, nr_1, nc_1;
    double x, y, dummy, sumprb;
    Rboolean lsm, lsp;

    nr_1 = *nrow - 1;
    nc_1 = *ncol - 1;

    ib = 0; /* -Wall */

    /* Construct random matrix */
    for (j = 0; j < nc_1; ++j)
        jwork[j] = ncolt[j];

    jc = *ntotal;

    for (l = 0; l < nr_1; ++l) { /* -----  matrix[ l, * ] ----- */
        ia = nrowt[l];
        ic = jc;
        jc -= ia;/* = n_tot - sum(nr[0:l]) */

        for (m = 0; m < nc_1; ++m) {
            id = jwork[m];
            ie = ic;
            ic -= id;
            ib = ie - ia;
            ii = ib - id;

            if (ie == 0) { /* Row [l,] is full, fill rest with zero entries */
                for (j = m; j < nc_1; ++j)
                    matrix[l + j * *nrow] = 0;
                ia = 0;
                break;
            }

            /* Generate pseudo-random number */
            dummy = unif_rand();

            do {/* Outer Loop */

                /* Compute conditional expected value of MATRIX(L, M) */

                nlm = (int)(ia * (id / (double) ie) + 0.5);
                x = exp(fact[ia] + fact[ib] + fact[ic] + fact[id]
                        - fact[ie] - fact[nlm]
                        - fact[id - nlm] - fact[ia - nlm] - fact[ii + nlm]);
                if (x >= dummy)
                    break;
                if (x == 0.)/* MM: I haven't seen this anymore */
                    error("rcont2 [%d,%d]: exp underflow to 0; algorithm failure", l, m);

                sumprb = x;
                y = x;
                nll = nlm;

                do {
                    /* Increment entry in row L, column M */
                    j = (int)((id - nlm) * (double)(ia - nlm));
                    lsp = (j == 0);
                    if (!lsp) {
                        ++nlm;
                        x = x * j / ((double) nlm * (ii + nlm));
                        sumprb += x;
                        if (sumprb >= dummy)
                            goto L160;
                    }

                    do {
                        R_CheckUserInterrupt();

                        /* Decrement entry in row L, column M */
                        j = (int)(nll * (double)(ii + nll));
                        lsm = (j == 0);
                        if (!lsm) {
                            --nll;
                            y = y * j / ((double) (id - nll) * (ia - nll));
                            sumprb += y;
                            if (sumprb >= dummy) {
                                nlm = nll;
                                goto L160;
                            }
                            /* else */
                            if (!lsp)
                                break;/* to while (!lsp) */
                        }
                    } while (!lsm);

                } while (!lsp);

                dummy = sumprb * unif_rand();

            } while (1);

L160:
            matrix[l + m * *nrow] = nlm;
            ia -= nlm;
            jwork[m] -= nlm;
        }
        matrix[l + nc_1 * *nrow] = ia;/* last column in row l */
    }

    /* Compute entries in last row of MATRIX */
    for (m = 0; m < nc_1; ++m)
        matrix[nr_1 + m * *nrow] = jwork[m];

    matrix[nr_1 + nc_1 * *nrow] = ib - matrix[nr_1 + (nc_1-1) * *nrow];

    return;
}

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
