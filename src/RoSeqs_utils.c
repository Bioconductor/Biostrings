/****************************************************************************
 *                       RoSeqs low-level utilities                         *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"
#include "IRanges_interface.h"
#include <S.h> /* for Salloc() */

static int debug = 0;

SEXP debug_RoSeqs_utils()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'RoSeqs_utils.c'\n",
		debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'RoSeqs_utils.c'\n");
#endif
	return R_NilValue;
}

RoSeqs _alloc_RoSeqs(int nelt)
{
	RoSeqs seqs;

	seqs.elts = Salloc((long) nelt, cachedCharSeq);
	seqs.nelt = nelt;
	return seqs;
}

