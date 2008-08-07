/****************************************************************************
 *                        RoSeq low-level utilities                         *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"

static int debug = 0;

SEXP debug_RoSeq_utils()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'RoSeq_utils.c'\n",
		debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'RoSeq_utils.c'\n");
#endif
	return R_NilValue;
}

SEXP _new_IRanges_from_RoSeqs(const char *class, RoSeqs seqs)
{
	const RoSeq *seq;
	SEXP start, width, ans;
	int *start_elt, *width_elt, *start_prev_elt, i;

#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _new_IRanges_from_RoSeqs(): BEGIN\n");
	}
#endif
	seq = seqs.elts;
	PROTECT(start = NEW_INTEGER(seqs.nelt));
	PROTECT(width = NEW_INTEGER(seqs.nelt));
	start_elt = INTEGER(start);
	width_elt = INTEGER(width);
	if (seqs.nelt >= 1) {
		*(start_elt++) = 1;
		*(width_elt++) = seq->nelt;
	}
	if (seqs.nelt >= 2)
		for (i = 1, start_prev_elt = INTEGER(start); i < seqs.nelt; i++) {
			*(start_elt++) = *(start_prev_elt++) + (seq++)->nelt;
			*(width_elt++) = seq->nelt;
		}
	PROTECT(ans = new_IRanges(class, start, width, R_NilValue));
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _new_IRanges_from_RoSeqs(): END\n");
	}
#endif
	UNPROTECT(3);
	return ans;
}

