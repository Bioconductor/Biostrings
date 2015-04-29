/****************************************************************************
 *                                                                          *
 *                   Routines for BitMatrix manipulation                    *
 *                           Author: Herve Pages                            *
 *                                                                          *
 ****************************************************************************/
#include "Biostrings.h"
#include "S4Vectors_interface.h"
#include <S.h> /* for Salloc() */

#include <stdio.h>
#include <limits.h> /* for CHAR_BIT and ULONG_MAX */
#include <stdlib.h> /* for div() */


static int debug = 0;

#define BITMATBYROW_NCOL (sizeof(int) * CHAR_BIT)

typedef IntAE BitMatByRow;

void _BitCol_set_val(BitCol *bitcol, BitWord val)
{
	div_t q;
	BitWord *bitword;
	int i1;

	q = div(bitcol->nbit, NBIT_PER_BITWORD);
	if (q.rem != 0)
		q.quot++;
	bitword = bitcol->bitword0;
	for (i1 = 0; i1 < q.quot; i1++)
		*(bitword++) = val;
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
	bitcol.bitword0 = Salloc((long) nword, BitWord);
	bitcol.nword = nword;
	bitcol.nbit = nbit;
	_BitCol_set_val(&bitcol, val);
	return bitcol;
}

int _BitCol_get_bit(const BitCol *bitcol, int i)
{
	div_t q;
	BitWord *bitword;

	q = div(i, NBIT_PER_BITWORD);
	bitword = bitcol->bitword0 + q.quot;
	return (*bitword >> q.rem) & 1UL;
}

void _BitCol_set_bit(BitCol *bitcol, int i, int bit)
{
	div_t q;
	BitWord *bitword, mask;

	q = div(i, NBIT_PER_BITWORD);
	bitword = bitcol->bitword0 + q.quot;
	mask = 1UL << q.rem;
	if (bit)
		*bitword |= mask;
	else
		*bitword &= ~mask;
	return;
}

void _BitCol_A_gets_BimpliesA(BitCol *A, const BitCol *B)
{
	div_t q;
	BitWord *Abitword;
	const BitWord *Bbitword;
	int i1;

	if (A->nbit != B->nbit)
		error("_BitCol_A_gets_BimpliesA(): "
		      "'A' and 'B' are incompatible");
	q = div(A->nbit, NBIT_PER_BITWORD);
	if (q.rem != 0)
		q.quot++;
	Abitword = A->bitword0;
	Bbitword = B->bitword0;
	for (i1 = 0; i1 < q.quot; i1++)
		*(Abitword++) |= ~(*(Bbitword++));
	return;
}

/* WARNING: This is a 0-copy column extraction! */
BitCol _BitMatrix_get_col(const BitMatrix *bitmat, int j)
{
	BitCol bitcol;

	bitcol.bitword0 = bitmat->bitword00 + j * bitmat->nword_per_col;
	bitcol.nword = bitmat->nword_per_col;
	bitcol.nbit = bitmat->nrow;
	return bitcol;
}

void _BitMatrix_set_col(BitMatrix *bitmat, int j, const BitCol *bitcol)
{
	div_t q;
	BitWord *Lbitword;
	const BitWord *Rbitword;
	int i1;

	if (bitmat->nrow != bitcol->nbit)
		error("_BitMatrix_set_col(): "
		      "'bitmat' and 'bitcol' are incompatible");
	q = div(bitmat->nrow, NBIT_PER_BITWORD);
	if (q.rem != 0)
		q.quot++;
	Lbitword = bitmat->bitword00 + j * bitmat->nword_per_col;
	Rbitword = bitcol->bitword0;
	for (i1 = 0; i1 < q.quot; i1++)
		*(Lbitword++) = *(Rbitword++);
	return;
}


/*
 * _BitMatrix_set_val() could also be implemented as:
 *   for (j = 0; j < bitmat->ncol; j++) {
 *     bitcol = _BitMatrix_get_col(bitmat, j);
 *     _BitCol_set_val(&bitcol, val);
 *   }
 * but the implementation below is faster (only 1 call to div()).
 */
void _BitMatrix_set_val(BitMatrix *bitmat, BitWord val)
{
	div_t q;
	int i1, j;
	BitWord *bitword0, *bitword;

	q = div(bitmat->nrow, NBIT_PER_BITWORD);
	if (q.rem != 0)
		q.quot++;
	for (j = 0, bitword0 = bitmat->bitword00;
	     j < bitmat->ncol;
	     j++, bitword0 += bitmat->nword_per_col)
	{
		bitword = bitword0;
		for (i1 = 0; i1 < q.quot; i1++)
			*(bitword++) = val;
	}
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
	bitmat.bitword00 = Salloc((long) nword, BitWord);
	bitmat.nword_per_col = nword_per_col;
	bitmat.nrow = nrow;
	bitmat.ncol = ncol;
	_BitMatrix_set_val(&bitmat, val);
	return bitmat;
}

/*
 * _BitMatrix_get_bit() could also be implemented as:
 *   bitcol = _BitMatrix_get_col(bitmat, j);
 *   return _BitCol_get_bit(&bitcol, i);
 */
int _BitMatrix_get_bit(const BitMatrix *bitmat, int i, int j)
{
	div_t q;
	BitWord *bitword;

	q = div(i, NBIT_PER_BITWORD);
	bitword = bitmat->bitword00 + j * bitmat->nword_per_col + q.quot;
	return (*bitword >> q.rem) & 1UL;
}

/*
 * _BitMatrix_set_bit() could also be implemented as:
 *   bitcol = _BitMatrix_get_col(bitmat, j);
 *   _BitCol_set_bit(&bitcol, i, bit);
 */
void _BitMatrix_set_bit(BitMatrix *bitmat, int i, int j, int bit)
{
	div_t q;
	BitWord *bitword, mask;

	q = div(i, NBIT_PER_BITWORD);
	bitword = bitmat->bitword00 + j * bitmat->nword_per_col + q.quot;
	mask = 1UL << q.rem;
	if (bit)
		*bitword |= mask;
	else
		*bitword &= ~mask;
	return;
}

void _BitMatrix_Rrot1(BitMatrix *bitmat)
{
	div_t q;
	BitWord *bitwordNE, *Lbitword, *Rbitword;
	int i1, j;

	if (bitmat->ncol == 0)
		error("_BitMatrix_Rrot1(): bitmat->ncol == 0");
	q = div(bitmat->nrow, NBIT_PER_BITWORD);
	if (q.rem != 0)
		q.quot++;
	bitwordNE = bitmat->bitword00 + (bitmat->ncol - 1) * bitmat->nword_per_col;
	for (i1 = 0; i1 < q.quot; i1++) {
		Rbitword = bitwordNE + i1;
		Lbitword = Rbitword - bitmat->nword_per_col;
		for (j = 1; j < bitmat->ncol; j++) {
			*Rbitword = *Lbitword;
			Rbitword = Lbitword;
			Lbitword -= bitmat->nword_per_col;
		}
		*Rbitword = ULONG_MAX;
	}
	return;
}

void _BitMatrix_grow1rows(BitMatrix *bitmat, const BitCol *bitcol)
{
	div_t q;
	BitWord *Lbitword, Rbitword, ret;
	int i1, j;

	if (bitmat->nrow != bitcol->nbit)
		error("_BitMatrix_grow1rows(): "
		      "'bitmat' and 'bitcol' are incompatible");
	q = div(bitmat->nrow, NBIT_PER_BITWORD);
	if (q.rem != 0)
		q.quot++;
	for (i1 = 0; i1 < q.quot; i1++) {
		Lbitword = bitmat->bitword00 + i1;
		Rbitword = bitcol->bitword0[i1];
		for (j = 0; j < bitmat->ncol; j++) {
			ret = *Lbitword & Rbitword; // and
			*Lbitword |= Rbitword; // or
			Rbitword = ret;
			Lbitword += bitmat->nword_per_col;
		}
	}
	return;
}



/*****************************************************************************
 * Testing and debugging stuff
 */

static void BitMatrix_tr(BitMatrix *in, BitMatByRow *out)
{
	BitWord rbit, *bitword;
	int i1, i2, i, j, cbit;

	if (in->nrow != IntAE_get_nelt(out))
		error("BitMatrix_tr(): 'in' and 'out' are incompatible");
	if (in->ncol >= BITMATBYROW_NCOL)
		error("BitMatrix_tr(): 'in' has too many columns");
	for (i1 = i = 0; i1 < in->nword_per_col; i1++) {
		for (i2 = 0, rbit = 1UL;
		     i2 < NBIT_PER_BITWORD;
		     i2++, i++, rbit <<= 1)
		{
			if (i >= in->nrow)
				return;
			out->elts[i] = 0;
			for (j = 0, bitword = in->bitword00 + i1, cbit = 1;
			     j < in->ncol;
			     j++, bitword += in->nword_per_col, cbit <<= 1) {
				if (*bitword & rbit)
					out->elts[i] += cbit;
			}
		}
	}
	return;
}

static void BitMatrix_print(BitMatrix *bitmat)
{
	BitMatByRow *bitmat_byrow;
	int i, *row, j, cbit, bit;

	bitmat_byrow = new_IntAE(bitmat->nrow, bitmat->nrow, 0);
	BitMatrix_tr(bitmat, bitmat_byrow);
	for (i = 0, row = bitmat_byrow->elts; i < bitmat->nrow; i++, row++) {
		Rprintf("%4d: ", i);
		for (j = 0, cbit = 1; j < bitmat->ncol; j++, cbit <<= 1) {
			bit = (*row & cbit) != 0;
			Rprintf("%d", bit);
		}
		Rprintf(" (%d)\n", *row);
	}
	return;
}

/*
static void BitMatrix_addcol(BitMatrix *bitmat, const BitCol *bitcol)
{
	BitWord *Lbitword, Rbitword, ret;
	int i1, j;

	if (bitmat->nrow != bitcol->nbit)
		error("BitMatrix_addcol(): "
		      "'bitmat' and 'bitcol' are incompatible");
	for (i1 = 0; i1 < bitmat->nword_per_col; i1++) {
		Lbitword = bitmat->bitword00 + i1;
		Rbitword = bitcol->bitword0[i1];
		for (j = 0; j < bitmat->ncol; j++) {
			ret = *Lbitword & Rbitword; // and
			*Lbitword ^= Rbitword; // xor
			Rbitword = ret;
			Lbitword += bitmat->nword_per_col;
		}
		if (Rbitword)
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

		bitmat0 = _new_BitMatrix(40, 5, 0UL);
		bitcol0 = _new_BitCol(40, 33UL + (1UL << 39));
		
		BitMatrix_print(&bitmat0);
		_BitMatrix_set_bit(&bitmat0, 0, 0, 1);
		BitMatrix_print(&bitmat0);
		_BitMatrix_set_bit(&bitmat0, 39, 4, 1);
		BitMatrix_print(&bitmat0);
		_BitMatrix_set_bit(&bitmat0, 39, 4, 0);
		BitMatrix_print(&bitmat0);
		_BitMatrix_grow1rows(&bitmat0, &bitcol0);
		BitMatrix_print(&bitmat0);
		_BitMatrix_grow1rows(&bitmat0, &bitcol0);
		BitMatrix_print(&bitmat0);
		_BitMatrix_grow1rows(&bitmat0, &bitcol0);
		BitMatrix_print(&bitmat0);
		_BitMatrix_grow1rows(&bitmat0, &bitcol0);
		BitMatrix_print(&bitmat0);
		_BitMatrix_Rrot1(&bitmat0);
		BitMatrix_print(&bitmat0);
		_BitMatrix_Rrot1(&bitmat0);
		BitMatrix_print(&bitmat0);
		_BitMatrix_Rrot1(&bitmat0);
		BitMatrix_print(&bitmat0);
		_BitMatrix_Rrot1(&bitmat0);
		BitMatrix_print(&bitmat0);
		_BitMatrix_Rrot1(&bitmat0);
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
		bitcol0.bitword0[0] = 33UL;
		bitcol0.bitword0[4] = 1UL << 43;
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

