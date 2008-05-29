/****************************************************************************
 *                        Fast SparseList utilities                         *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"

static int debug = 0;

SEXP debug_SparseList_utils()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'SparseList_utils.c'\n",
		debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'SparseList_utils.c'\n");
#endif
	return R_NilValue;
}

/* 'symbol' must be a CHARSXP */
SEXP getSymbolVal(SEXP symbol, SEXP envir, int error_on_unbound_value)
{
	SEXP ans;

	/* The following code was inspired by R's do_get() code.
	 * Note that do_get() doesn't use PROTECT at all and so do we...
	 */
	ans = findVar(install(translateChar(symbol)), envir);
	if (ans == R_UnboundValue) {
		if (error_on_unbound_value)
			error("Biostrings internal error in getSymbolVal(): unbound value");
		return R_UnboundValue;
	}
	if (TYPEOF(ans) == PROMSXP)
		ans = eval(ans, envir);
	if (ans != R_NilValue && NAMED(ans) == 0)
		SET_NAMED(ans, 1);
	return ans;
}

