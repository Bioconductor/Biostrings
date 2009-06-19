#include <stdio.h>
#include <limits.h> /* for CHAR_BIT */
#include <stdlib.h> /* for div() */

typedef long unsigned int BitWord;

#define NBIT_PER_BITWORD (sizeof(BitWord) * CHAR_BIT)
#define BITMATBYROW_NCOL (sizeof(int) * CHAR_BIT)

typedef struct bitcol {
	BitWord *words;
	int nword;
	int nbit; // <= nword * NBIT_PER_BITWORD
} BitCol;
	
typedef struct bitmat_bycol {
	BitWord *words;
	int nword_per_col;
	int nrow; // <= nword_per_col * NBIT_PER_BITWORD
	int ncol;
} BitMatByCol;

typedef struct bitmat_byrow {
	int *rows;
	int nrow;
} BitMatByRow;

BitCol new_BitCol(int nbit)
{
	BitCol bitcol;
	div_t q;
	int nword, j;
	BitWord *word;

	q = div(nbit, NBIT_PER_BITWORD);
	nword = q.quot;
	if (q.rem != 0)
		nword++;
	bitcol.words = (BitWord *) malloc(nword * sizeof(BitWord));
	if (bitcol.words == NULL) {
		printf("malloc() error in new_BitCol()\n");
		exit(-1);
	}
	/* Initialization */
	for (j = 0, word = bitcol.words; j < nword; j++, word++)
		*word = 0UL;
	bitcol.nword = nword;
	bitcol.nbit = nbit;
	return bitcol;
}

void BitMatByCol_reset(BitMatByCol *bitmat)
{
	int nword, i1;
	BitWord *word;

	nword = bitmat->nword_per_col * bitmat->ncol;
	for (i1 = 0, word = bitmat->words; i1 < nword; i1++, word++)
		*word = 0UL;
	return;
}

BitMatByCol new_BitMatByCol(int nrow, int ncol)
{
	BitMatByCol bitmat;
	div_t q;
	int nword_per_col, nword;

	q = div(nrow, NBIT_PER_BITWORD);
	nword_per_col = q.quot;
	if (q.rem != 0)
		nword_per_col++;
	nword = nword_per_col * ncol;
	bitmat.words = (BitWord *) malloc(nword * sizeof(BitWord));
	if (bitmat.words == NULL) {
		printf("malloc() error in new_BitMatByCol()\n");
		exit(-1);
	}
	bitmat.nword_per_col = nword_per_col;
	bitmat.nrow = nrow;
	bitmat.ncol = ncol;
	BitMatByCol_reset(&bitmat);
	return bitmat;
}

void BitMatByRow_reset(BitMatByRow *bitmat)
{
	int i, *row;

	for (i = 0, row = bitmat->rows; i < bitmat->nrow; i++, row++)
		*row = 0;
	return;
}

BitMatByRow new_BitMatByRow(int nrow)
{
	BitMatByRow bitmat;

	bitmat.rows = (int *) malloc(nrow * sizeof(int));
	if (bitmat.rows == NULL) {
		printf("malloc() error in new_BitMatByRow()\n");
		exit(-1);
	}
	bitmat.nrow = nrow;
	BitMatByRow_reset(&bitmat);
	return bitmat;
}

void BitMatByCol_tr(BitMatByCol *in, BitMatByRow *out)
{
	BitWord rbit, *word;
	int i1, i2, i, j, cbit;

	if (in->nrow != out->nrow) {
		printf("BitMatByCol_tr(): in and out are incompatible\n");
		exit(-1);
	}
	if (in->ncol >= BITMATBYROW_NCOL) {
		printf("BitMatByCol_tr(): in has too many columns\n");
		exit(-1);
	}
	for (i1 = i = 0; i1 < in->nword_per_col; i1++) {
		for (i2 = 0, rbit = 1UL;
		     i2 < NBIT_PER_BITWORD;
		     i2++, i++, rbit <<= 1)
		{
			if (i >= in->nrow)
				return;
			out->rows[i] = 0;
			word = in->words + i1;
			for (j = 0, word = in->words + i1, cbit = 1;
			     j < in->ncol;
			     j++, word += in->nword_per_col, cbit <<= 1) {
				if (*word & rbit)
					out->rows[i] += cbit;
			}
		}
	}
	return;
}

void BitMatByCol_print(BitMatByCol *bitmat)
{
	BitMatByRow bitmat_byrow;
	int i, *row, j, cbit, bit;

	bitmat_byrow = new_BitMatByRow(bitmat->nrow);
	BitMatByCol_tr(bitmat, &bitmat_byrow);
	for (i = 0, row = bitmat_byrow.rows; i < bitmat_byrow.nrow; i++, row++) {
		printf("%4d: %3d (", i, *row);
		for (j = 0, cbit = 1; j < bitmat->ncol; j++, cbit <<= 1) {
			bit = (*row & cbit) != 0;
			printf("%d", bit);
		}
		printf(")\n");
	}
	return;
}

void BitMatByCol_addcol(BitMatByCol *bitmat, BitCol *bitcol)
{
	BitWord *Lword, Rword, ret;
	int i1, j;

	if (bitmat->nrow != bitcol->nbit) {
		printf("BitMatByCol_addcol(): bitmat and bitcol are incompatible\n");
		exit(-1);
	}
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
			printf("warning: overflow in BitMatByCol_addcol()\n");
	}
	return;
}

/*
void testing1(BitMatByCol *bitmat, BitCol *bitcol)
{
	int n, p;

	for (n = 0; n < 1000000; n++) {
		BitMatByCol_reset(bitmat);
		for (p = 0; p < 25; p++)
			BitMatByCol_addcol(bitmat, bitcol);
	}
	//BitMatByCol_print(&bitmat);
	return;
}

void testing2(BitMatByRow *bitmat_byrow)
{
	int n, p, i, *row;

	for (n = 0; n < 1000000; n++) {
		BitMatByRow_reset(bitmat_byrow);
		for (i = 0, row = bitmat_byrow->rows; i < bitmat_byrow->nrow; i++)
			for (p = 0; p < 25; p++)
				(*row)++;
	}
	return;
}

int main()
{
	BitMatByCol bitmat0;
	BitCol bitcol0;
	BitMatByRow bitmat_byrow0;

	bitmat0 = new_BitMatByCol(3000, 5);
	bitcol0 = new_BitCol(3000);
	bitcol0.words[0] = 33UL;
	bitcol0.words[4] = 1UL << 43;
	bitmat_byrow0 = new_BitMatByRow(3000);
	//BitMatByCol_print(&bitmat0);
	//BitMatByCol_addcol(&bitmat0, &bitcol0);
	//BitMatByCol_print(&bitmat0);
	//BitMatByCol_addcol(&bitmat0, &bitcol0);
	//BitMatByCol_print(&bitmat0);
	//BitMatByCol_addcol(&bitmat0, &bitcol0);
	//BitMatByCol_print(&bitmat0);
	testing1(&bitmat0, &bitcol0);
	//testing2(&bitmat_byrow0);
	return 0;
}
*/

