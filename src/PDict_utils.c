/****************************************************************************
 *                     PDict object utility functions                       *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"
#include <S.h> /* for Salloc()/Srealloc() */

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


/****************************************************************************
 * Buffer of duplicates.
 */

static IntBuf dup2unq_buf;

void _init_dup2unq_buf(int length)
{
	dup2unq_buf = _new_IntBuf(length, length, NA_INTEGER);
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
	return _IntBuf_asINTEGER(&dup2unq_buf);
}

