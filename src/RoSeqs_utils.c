/****************************************************************************
 *                        RoSeqs low-level utilities                        *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "Biostrings.h"
#include "IRanges_interface.h"

#define Salloc(n,t) (t*)S_alloc(n, sizeof(t))  /* from old <S.h> */

RoSeqs _alloc_RoSeqs(int nelt)
{
	RoSeqs seqs;

	seqs.elts = Salloc((long) nelt, Chars_holder);
	seqs.nelt = nelt;
	return seqs;
}

