/****************************************************************************
 *                     PDict object utility functions                       *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"
#include "IRanges_interface.h"

static int debug = 0;

SEXP debug_PDict_utils()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'PDict_utils.c'\n",
		debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'PDict_utils.c'\n");
#endif
	return R_NilValue;
}


SEXP Dups_diff(SEXP x_unq2dup, SEXP y_dup2unq)
{
	SEXP ans, ans_elt, dups;
	int ans_length, i, j;
	const int *dup;
	IntAE new_dups;

	new_dups = new_IntAE(0, 0, 0);
	ans_length = LENGTH(x_unq2dup);
	PROTECT(ans = NEW_LIST(ans_length));
	for (i = 0; i < ans_length; i++) {
		dups = VECTOR_ELT(x_unq2dup, i);
		if (dups == R_NilValue)
			continue;
		new_dups.nelt = 0;
		for (j = 0, dup = INTEGER(dups); j < LENGTH(dups); j++, dup++) {
			if (INTEGER(y_dup2unq)[*dup - 1] == NA_INTEGER)
				IntAE_insert_at(&new_dups, new_dups.nelt, *dup);
		}
		PROTECT(ans_elt = IntAE_asINTEGER(&new_dups));
		SET_ELEMENT(ans, i, ans_elt);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * Buffer of duplicates.
 */

static IntAE dup2unq_buf;

void _init_dup2unq_buf(int length)
{
	dup2unq_buf = new_IntAE(length, length, NA_INTEGER);
	return;
}

void _report_dup(int poffset, int P_id)
{
	dup2unq_buf.elts[poffset] = P_id;
	return;
}


/****************************************************************************
 * Turning our local data structures into an R list (SEXP).
 *
 * Note that none of the functions below is a .Call() entry point!
 */

/* NOT a .Call() entry point! */
SEXP _dup2unq_asINTEGER()
{
	return IntAE_asINTEGER(&dup2unq_buf);
}

