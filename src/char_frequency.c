#include "Biostrings.h"
#include <string.h>


static int chrtrtable[CHRTRTABLE_LENGTH];

static void add_char_freqs(RoSeq X, int *freqs)
{
	int i;

	for (i = 0; i < X.nelt; i++, X.elts++)
		freqs[(unsigned char) *X.elts]++;
	return;
}

static void add_code_freqs(RoSeq X, int *freqs)
{
	int i, offset;

	for (i = 0; i < X.nelt; i++, X.elts++) {
		offset = chrtrtable[(unsigned char) *X.elts];
		if (offset == -1)
			continue;
		freqs[(unsigned char) offset]++;
	}
	return;
}

/*
 * --- .Call ENTRY POINT ---
 * Return a vector of 256 integers.
 */
SEXP XString_char_frequency(SEXP x)
{
	SEXP ans;
	RoSeq X;

	PROTECT(ans = NEW_INTEGER(256));
	memset(INTEGER(ans), 0, LENGTH(ans) * sizeof(int));
	X = _get_XString_asRoSeq(x);
	add_char_freqs(X, INTEGER(ans));
	UNPROTECT(1);
	return ans;
}

/*
 * --- .Call ENTRY POINT ---
 * Return a vector of LENGTH(codes) integers.
 */
SEXP XString_code_frequency(SEXP x, SEXP codes)
{
	SEXP ans;
	RoSeq X;

	_init_chrtrtable(INTEGER(codes), LENGTH(codes), chrtrtable);
	PROTECT(ans = NEW_INTEGER(LENGTH(codes)));
	memset(INTEGER(ans), 0, LENGTH(ans) * sizeof(int));
	X = _get_XString_asRoSeq(x);
	add_code_freqs(X, INTEGER(ans));
	
	UNPROTECT(1);
	return ans;
}

/*
 * --- .Call ENTRY POINT ---
 */
SEXP XStringSet_char_frequency(SEXP x, SEXP collapse)
{
	SEXP ans;
	int x_length, *freqs, i;
	RoSeq X;

	x_length = _get_XStringSet_length(x);
	if (LOGICAL(collapse)[0])
		PROTECT(ans = NEW_INTEGER(256));
	else
		PROTECT(ans = NEW_INTEGER(256 * x_length));
	freqs = INTEGER(ans);
	memset(freqs, 0, LENGTH(ans) * sizeof(int));
	for (i = 0; i < x_length; i++) {
		X = _get_XStringSet_elt_asRoSeq(x, i);
		add_char_freqs(X, freqs);
		if (!LOGICAL(collapse)[0])
			freqs += 256;
	}
	UNPROTECT(1);
	return ans;
}

/*
 * --- .Call ENTRY POINT ---
 */
SEXP XStringSet_code_frequency(SEXP x, SEXP collapse, SEXP codes)
{
	SEXP ans;
	int x_length, *freqs, i;
	RoSeq X;

	_init_chrtrtable(INTEGER(codes), LENGTH(codes), chrtrtable);
	x_length = _get_XStringSet_length(x);
	if (LOGICAL(collapse)[0])
		PROTECT(ans = NEW_INTEGER(LENGTH(codes)));
	else
		PROTECT(ans = NEW_INTEGER(LENGTH(codes) * x_length));
	freqs = INTEGER(ans);
	memset(freqs, 0, LENGTH(ans) * sizeof(int));
	for (i = 0; i < x_length; i++) {
		X = _get_XStringSet_elt_asRoSeq(x, i);
		add_code_freqs(X, freqs);
		if (!LOGICAL(collapse)[0])
			freqs += LENGTH(codes);
	}
	UNPROTECT(1);
	return ans;
}

/*
 * --- .Call ENTRY POINT ---
 */
SEXP oligonucleotide_frequency(SEXP x_xp, SEXP x_offset, SEXP x_length,
		SEXP base_codes, SEXP width, SEXP fast_moving_side)
{
	SEXP ans;
	const Rbyte *x;
	int ans_len, ans_offset_bitmask, ans_offset, x_len, width0, nbit_in_mask,
            right_moves_fastest, i, nb_valid_left_char, twobit;
	static int eightbit2twobit_lkup[256];

	x = RAW(R_ExternalPtrTag(x_xp)) + INTEGER(x_offset)[0];
	x_len = INTEGER(x_length)[0];
	if (LENGTH(base_codes) != 4)
		error("'base_codes' must be of length 4");
	_init_chrtrtable(INTEGER(base_codes), LENGTH(base_codes), eightbit2twobit_lkup);
	width0 = INTEGER(width)[0];
	if (width0 < 1 || width0 > 12)
		error("'width' must be >=1 and <= 12");
	right_moves_fastest = strcmp(CHAR(STRING_ELT(fast_moving_side, 0)), "right") == 0;
	ans_len = 1 << (width0 * 2);
	nbit_in_mask = (width0 - 1) * 2;
	ans_offset_bitmask = (1 << nbit_in_mask) - 1;
	if (!right_moves_fastest)
		ans_offset_bitmask <<= 2;
	PROTECT(ans = NEW_INTEGER(ans_len));
	memset(INTEGER(ans), 0, LENGTH(ans) * sizeof(int));
	nb_valid_left_char = 0;
	for (i = 0; i < x_len; i++, x++) {
		twobit = eightbit2twobit_lkup[(unsigned char) *x];
		if (twobit == -1) {
			nb_valid_left_char = 0;
		} else {
			nb_valid_left_char++;
			ans_offset &= ans_offset_bitmask;
			if (right_moves_fastest) {
				ans_offset <<= 2;
			} else {
				ans_offset >>= 2;
				twobit <<= nbit_in_mask;
			}
			ans_offset += twobit;
			if (nb_valid_left_char >= width0)
				INTEGER(ans)[ans_offset]++;
		}
	}
	UNPROTECT(1);
	return ans;
}

