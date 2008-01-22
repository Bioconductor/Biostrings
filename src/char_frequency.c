#include "Biostrings.h"
#include <string.h>


/*
 * Arguments:
 *   'x_xp': x@data@xp
 *   'x_offset': x@offset
 *   'x_length': x@length
 * Returns a vector of 256 integers.
 */
SEXP char_frequency(SEXP x_xp, SEXP x_offset, SEXP x_length)
{
	SEXP ans;
	const Rbyte *x;
	int x_len, i;

	x = RAW(R_ExternalPtrTag(x_xp)) + INTEGER(x_offset)[0];
	x_len = INTEGER(x_length)[0];
	PROTECT(ans = NEW_INTEGER(256));
	memset(INTEGER(ans), 0, LENGTH(ans) * sizeof(int));
	for (i = 0; i < x_len; i++, x++)
		INTEGER(ans)[(unsigned char) *x]++;
	UNPROTECT(1);
	return ans;
}

SEXP all_oligonucleotide_frequency(SEXP x_xp, SEXP x_offset, SEXP x_length,
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
	_init_code2offset_lkup(INTEGER(base_codes), LENGTH(base_codes), eightbit2twobit_lkup);
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

