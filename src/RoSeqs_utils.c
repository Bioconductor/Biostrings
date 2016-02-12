/****************************************************************************
 *                        RoSeqs low-level utilities                        *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "Biostrings.h"
#include "IRanges_interface.h"
#include <S.h> /* for Salloc() */

RoSeqs _alloc_RoSeqs(int nelt)
{
	RoSeqs seqs;

	seqs.elts = Salloc((long) nelt, Chars_holder);
	seqs.nelt = nelt;
	return seqs;
}

