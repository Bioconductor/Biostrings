/****************************************************************************
 *             Turning an SEXP object containing input sequences            *
 *                        into an array of CharSeq                          *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"
#include <S.h> /* for Salloc() and Srealloc() */

static int debug = 0;

SEXP Biostrings_debug_SEXP_to_charseqs()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'SEXP_to_charseqs.c'\n", debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'SEXP_to_charseqs.c'\n");
#endif
	return R_NilValue;
}

/*
 * 'start' and 'nchar' must have the same length as 'x'.
 */
const CharSeq *STRSXP_to_charseqs(SEXP x, const int *start, const int *nchar, int *nseq)
{
	CharSeq *seqs, *seq;
	int i, offset;
	const int *start_p, *nchar_p;
	SEXP x_elt;

	*nseq = LENGTH(x);
	seqs = Salloc((long) *nseq, CharSeq);
	for (i = 0, start_p = start, nchar_p = nchar, seq = seqs;
	     i < *nseq;
	     i++, start_p++, nchar_p++, seq++) {
		x_elt = STRING_ELT(x, i);
		if (x_elt == NA_STRING)
			error("input sequence %d is NA", i+1);
		seq->data = CHAR(x_elt);
		seq->data += offset = _start2offset(*start_p);
		seq->length = _nchar2length(*nchar_p, offset, LENGTH(x_elt));
	}
	return seqs;
}

const CharSeq *BStringSet_to_charseqs(SEXP x, const int *start, const int *nchar, int *nseq)
{
	CharSeq *seqs, *seq;
	int i, offset;
	const int *start_p, *nchar_p;

	*nseq = _get_BStringSet_length(x);
	seqs = Salloc((long) *nseq, CharSeq);
	for (i = 0, start_p = start, nchar_p = nchar, seq = seqs;
	     i < *nseq;
	     i++, start_p++, nchar_p++, seq++) {
		seq->data = _get_BStringSet_charseq(x, i, &(seq->length));
		seq->data += offset = _start2offset(*start_p);
		seq->length = _nchar2length(*nchar_p, offset, seq->length);
	}
	return seqs;
}

