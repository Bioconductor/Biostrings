/****************************************************************************
 *                        Fast SparseList utilities                         *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"
#include "IRanges_interface.h"

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

SEXP _SparseList_int2symb(int symb_as_int)
{
	char symbbuf[11];

	snprintf(symbbuf, sizeof(symbbuf), "%010d", symb_as_int);
	return mkChar(symbbuf); /* UNPROTECTED! */
}

int _SparseList_symb2int(SEXP symbol)
{
	int symb_as_int;

	sscanf(CHAR(symbol), "%d", &symb_as_int);
	return symb_as_int;
}

/* 'symbol' must be a CHARSXP */
SEXP _get_val_from_env(SEXP symbol, SEXP env, int error_on_unbound_value)
{
	SEXP ans;

	/* The following code was inspired by R's do_get() code.
	 * Note that do_get() doesn't use PROTECT at all and so do we...
	 */
	ans = findVar(install(translateChar(symbol)), env);
	if (ans == R_UnboundValue) {
		if (error_on_unbound_value)
			error("Biostrings internal error in _get_val_from_env(): "
			      "unbound value");
		return R_UnboundValue;
	}
	if (TYPEOF(ans) == PROMSXP)
		ans = eval(ans, env);
	if (ans != R_NilValue && NAMED(ans) == 0)
		SET_NAMED(ans, 1);
	return ans;
}

SEXP _get_val_from_SparseList(int symb_as_int, SEXP env, int error_on_unbound_value)
{
	SEXP symbol, ans;

	PROTECT(symbol = _SparseList_int2symb(symb_as_int));
	ans = _get_val_from_env(symbol, env, error_on_unbound_value);
	UNPROTECT(1);
	return ans;
}

int _get_int_from_SparseList(int symb_as_int, SEXP env)
{
	SEXP value;
	int val;

	value = _get_val_from_SparseList(symb_as_int, env, 0);
	if (value == R_UnboundValue)
		return NA_INTEGER;
	if (LENGTH(value) != 1)
		error("Biostrings internal error in _get_int_from_SparseList(): "
		      "value is not a single integer");
	val = INTEGER(value)[0];
	if (val == NA_INTEGER)
		error("Biostrings internal error in _get_int_from_SparseList(): "
		      "value is NA");
	return val;
}

void _set_env_from_IntAE(SEXP env, IntAE *int_ae)
{
	int symb_as_int, *elt;
	SEXP symbol, value;

	for (symb_as_int = 1, elt = int_ae->elts;
	     symb_as_int <= int_ae->nelt;
	     symb_as_int++, elt++)
	{
		if (*elt == NA_INTEGER)
			continue;
		PROTECT(symbol = _SparseList_int2symb(symb_as_int));
		PROTECT(value = ScalarInteger(*elt));
		defineVar(install(translateChar(symbol)), value, env);
		UNPROTECT(2);
	}
	return;
}

void _set_env_from_IntAEAE(SEXP env, IntAEAE *int_aeae)
{
	int symb_as_int;
	IntAE *elt;
	SEXP symbol, value;

	for (symb_as_int = 1, elt = int_aeae->elts;
	     symb_as_int <= int_aeae->nelt;
	     symb_as_int++, elt++)
	{
		if (elt->nelt == 0)
			continue;
		PROTECT(symbol = _SparseList_int2symb(symb_as_int));
		PROTECT(value = IntAE_asINTEGER(elt));
		defineVar(install(translateChar(symbol)), value, env);
		UNPROTECT(2);
	}
	return;
}

