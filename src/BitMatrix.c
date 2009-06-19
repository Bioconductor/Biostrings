/****************************************************************************
 *                                                                          *
 *                   Routines for BitMatrix manipulation                    *
 *                           Author: Herve Pages                            *
 *                                                                          *
 ****************************************************************************/
#include "Biostrings.h"
#include "IRanges_interface.h"
#include <S.h> /* for Salloc() */

#include <stdio.h>
#include <limits.h> /* for CHAR_BIT */
#include <stdlib.h> /* for div() */


static int debug = 0;

#define BITMATBYROW_NCOL (sizeof(int) * CHAR_BIT)

typedef IntAE BitMatByRow;

void _BitCol_set_val(const BitCol *bitcol, BitWord val)
{
	int i;
	BitWord *word;

	for (i = 0, word = bitcol->words; i < bitcol->nword; i++, word++)
		*word = val;
	return;
}

BitCol _new_BitCol(int nbit, BitWord val)
{
	BitCol bitcol;
	div_t q;
	int nword;

	if (nbit <= 0)
		error("_new_BitCol(): nbit <= 0");
	q = div(nbit, NBIT_PER_BITWORD);
	nword = q.quot;
	if (q.rem != 0)
		nword++;
	bitcol.words = Salloc((long) nword, BitWord);
	bitcol.nword = nword;
	bitcol.nbit = nbit;
	_BitCol_set_val(&bitcol, val);
	return bitcol;
}

void _BitMatrix_set_val(const BitMatrix *bitmat, BitWord val)
{
	int nword, i;
	BitWord *word;

	nword = bitmat->nword_per_col * bitmat->ncol;
	for (i = 0, word = bitmat->words; i < nword; i++, word++)
		*word = val;
	return;
}

BitMatrix _new_BitMatrix(int nrow, int ncol, BitWord val)
{
	BitMatrix bitmat;
	div_t q;
	int nword_per_col, nword;

	if (nrow <= 0 || ncol <= 0)
		error("_new_BitMatrix(): nrow <= 0 || ncol <= 0");
	q = div(nrow, NBIT_PER_BITWORD);
	nword_per_col = q.quot;
	if (q.rem != 0)
		nword_per_col++;
	nword = nword_per_col * ncol;
	bitmat.words = Salloc((long) nword, BitWord);
	bitmat.nword_per_col = nword_per_col;
	bitmat.nrow = nrow;
	bitmat.ncol = ncol;
	_BitMatrix_set_val(&bitmat, val);
	return bitmat;
}

static void BitMatrix_tr(BitMatrix *in, BitMatByRow *out)
{
	BitWord rbit, *word;
	int i1, i2, i, j, cbit;

	if (in->nrow != out->nelt)
		error("BitMatrix_tr(): in and out are incompatible");
	if (in->ncol >= BITMATBYROW_NCOL)
		error("BitMatrix_tr(): in has too many columns");
	for (i1 = i = 0; i1 < in->nword_per_col; i1++) {
		for (i2 = 0, rbit = 1UL;
		     i2 < NBIT_PER_BITWORD;
		     i2++, i++, rbit <<= 1)
		{
			if (i >= in->nrow)
				return;
			out->elts[i] = 0;
			word = in->words + i1;
			for (j = 0, word = in->words + i1, cbit = 1;
			     j < in->ncol;
			     j++, word += in->nword_per_col, cbit <<= 1) {
				if (*word & rbit)
					out->elts[i] += cbit;
			}
		}
	}
	return;
}

static void BitMatrix_print(BitMatrix *bitmat)
{
	BitMatByRow bitmat_byrow;
	int i, *row, j, cbit, bit;

	bitmat_byrow = new_IntAE(bitmat->nrow, bitmat->nrow, 0);
	BitMatrix_tr(bitmat, &bitmat_byrow);
	for (i = 0, row = bitmat_byrow.elts; i < bitmat_byrow.nelt; i++, row++) {
		Rprintf("%4d: %3d (", i, *row);
		for (j = 0, cbit = 1; j < bitmat->ncol; j++, cbit <<= 1) {
			bit = (*row & cbit) != 0;
			Rprintf("%d", bit);
		}
		Rprintf(")\n");
	}
	return;
}

void _BitMatrix_grow1rows(BitMatrix *bitmat, BitCol *bitcol)
{
	BitWord *Lword, Rword, ret;
	int i1, j;

	if (bitmat->nrow != bitcol->nbit)
		error("BitMatrix_addcol(): bitmat and bitcol are incompatible");
	for (i1 = 0; i1 < bitmat->nword_per_col; i1++) {
		Lword = bitmat->words + i1;
		Rword = bitcol->words[i1];
		for (j = 0; j < bitmat->ncol; j++) {
			ret = *Lword & Rword; // and
			*Lword |= Rword; // or
			Rword = ret;
			Lword += bitmat->nword_per_col;
		}
	}
	return;
}

/*
static void BitMatrix_addcol(BitMatrix *bitmat, BitCol *bitcol)
{
	BitWord *Lword, Rword, ret;
	int i1, j;

	if (bitmat->nrow != bitcol->nbit)
		error("BitMatrix_addcol(): bitmat and bitcol are incompatible");
	for (i1 = 0; i1 < bitmat->nword_per_col; i1++) {
		Lword = bitmat->words + i1;
		Rword = bitcol->words[i1];
		for (j = 0; j < bitmat->ncol; j++) {
			ret = *Lword & Rword; // and
			*Lword ^= Rword; // xor
			Rword = ret;
			Lword += bitmat->nword_per_col;
		}
		if (Rword)
			warning("integer overflow in BitMatrix_addcol()");
	}
	return;
}

static void testing1(BitMatrix *bitmat, BitCol *bitcol)
{
	int n, p;

	for (n = 0; n < 1000000; n++) {
		_BitMatrix_set_val(bitmat, 0UL);
		for (p = 0; p < 25; p++)
			BitMatrix_addcol(bitmat, bitcol);
	}
	//BitMatrix_print(&bitmat);
	return;
}

static void testing2(BitMatByRow *bitmat_byrow)
{
	int n, p, i, *row;

	for (n = 0; n < 1000000; n++) {
		IntAE_set_val(bitmat_byrow, 0);
		for (i = 0, row = bitmat_byrow->elts; i < bitmat_byrow->nelt; i++)
			for (p = 0; p < 25; p++)
				(*row)++;
	}
	return;
}
*/

SEXP debug_BitMatrix()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in file %s\n",
		debug ? "on" : "off", __FILE__);
	if (debug) {
		BitMatrix bitmat0;
		BitCol bitcol0;
		//BitMatByRow bitmat_byrow0;

		bitmat0 = _new_BitMatrix(40, 15, 0UL);
		bitcol0 = _new_BitCol(40, 33UL + (1UL << 39));
		
		BitMatrix_print(&bitmat0);
		_BitMatrix_grow1rows(&bitmat0, &bitcol0);
		BitMatrix_print(&bitmat0);
		_BitMatrix_grow1rows(&bitmat0, &bitcol0);
		BitMatrix_print(&bitmat0);
		_BitMatrix_grow1rows(&bitmat0, &bitcol0);
		BitMatrix_print(&bitmat0);
		_BitMatrix_grow1rows(&bitmat0, &bitcol0);
		BitMatrix_print(&bitmat0);
/*
		BitMatrix_print(&bitmat0);
		BitMatrix_addcol(&bitmat0, &bitcol0);
		BitMatrix_print(&bitmat0);
		BitMatrix_addcol(&bitmat0, &bitcol0);
		BitMatrix_print(&bitmat0);
		BitMatrix_addcol(&bitmat0, &bitcol0);
		BitMatrix_print(&bitmat0);
		BitMatrix_addcol(&bitmat0, &bitcol0);
		BitMatrix_print(&bitmat0);
*/

/*
		bitmat0 = _new_BitMatrix(3000, 5, 0UL);
		bitcol0 = _new_BitCol(3000, 0UL);
		bitcol0.words[0] = 33UL;
		bitcol0.words[4] = 1UL << 43;
		bitmat_byrow0 = new_IntAE(3000, 3000, 0);
		testing1(&bitmat0, &bitcol0);
		//testing2(&bitmat_byrow0);
*/
	}
#else
	Rprintf("Debug mode not available in file %s\n", __FILE__);
#endif
	return R_NilValue;
}

