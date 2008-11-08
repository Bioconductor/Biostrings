/****************************************************************************
               A SIMPLE MATCHING ALGO WITH SUPPORT FOR INDELS
	                     Author: Herve Pages
 ****************************************************************************/
#include "Biostrings.h"

static int debug = 0;

SEXP debug_match_pattern_indels()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in file %s\n",
		debug ? "on" : "off", __FILE__);
#else
	Rprintf("Debug mode not available in file %s\n", __FILE__);
#endif
	return R_NilValue;
}

static ByteTrTable byte2offset;

void _match_pattern_indels(const RoSeq *P, const RoSeq *S,
		int max_mm, int fixedP, int fixedS)
{
	if (P->nelt <= 0)
		error("empty pattern");
	_init_byte2offset_with_RoSeq(byte2offset, P, 0);
	return;
}

