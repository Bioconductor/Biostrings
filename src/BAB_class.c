/****************************************************************************
 *                    Basic manipulation of BAB objects                     *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "Biostrings.h"
#include "IRanges_interface.h"


SEXP IntegerBAB_new(SEXP max_nblock)
{
	SEXP blocks, prot, xp, classdef, ans;

	PROTECT(blocks = NEW_LIST(INTEGER(max_nblock)[0]));
	PROTECT(prot = NEW_INTEGER(2));
	INTEGER(prot)[0] = 0;  // nblock
	PROTECT(xp = R_MakeExternalPtr(NULL, blocks, prot));
	PROTECT(classdef = MAKE_CLASS("IntegerBAB"));
	PROTECT(ans = NEW_OBJECT(classdef));
	SET_SLOT(ans, mkChar("xp"), xp);
	UNPROTECT(5);
	return ans;
}

int *_get_BAB_nblock_ptr(SEXP x)
{
	SEXP xp, prot;

	xp = GET_SLOT(x, install("xp"));
	prot = R_ExternalPtrProtected(xp);
	return INTEGER(prot);
}

int *_get_BAB_lastblock_nelt_ptr(SEXP x)
{
	return _get_BAB_nblock_ptr(x) + 1;
}

SEXP _get_BAB_blocks(SEXP x)
{
	SEXP xp;

	xp = GET_SLOT(x, install("xp"));
	return R_ExternalPtrTag(xp);
}

SEXP _IntegerBAB_addblock(SEXP x, int block_length)
{
	SEXP xp, blocks, prot, block;
	int max_nblock, nblock; 

	xp = GET_SLOT(x, install("xp"));
	blocks = R_ExternalPtrTag(xp);
	max_nblock = LENGTH(blocks);
	prot = R_ExternalPtrProtected(xp);
	nblock = INTEGER(prot)[0];
	if (nblock >= max_nblock)
		error("_IntegerBAB_addblock(): reached max buffer size");
	PROTECT(block = NEW_INTEGER(block_length));
	SET_ELEMENT(blocks, nblock, block);
	UNPROTECT(1);
	nblock++;
	INTEGER(prot)[0] = nblock;
	INTEGER(prot)[1] = 0;  // lastblock_nelt
	return block;
}

