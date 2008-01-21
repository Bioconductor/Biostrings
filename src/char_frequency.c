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
		SEXP base_codes, SEXP width)
{
	SEXP ans;
	const Rbyte *x;
	int ans_len, ans_offset_bitmask, ans_offset, x_len, i, nb_valid_left_char, twobit;
	static int eightbit2twobit_lkup[256];

	x = RAW(R_ExternalPtrTag(x_xp)) + INTEGER(x_offset)[0];
	x_len = INTEGER(x_length)[0];
	if (LENGTH(base_codes) != 4)
		error("'base_codes' must be of length 4");
	_init_code2offset_lkup(INTEGER(base_codes), LENGTH(base_codes), eightbit2twobit_lkup);
	if (INTEGER(width)[0] < 1 || INTEGER(width)[0] > 12)
		error("'width' must be >=1 and <= 12");
	ans_len = 1 << (INTEGER(width)[0] * 2);
	ans_offset_bitmask = ans_len - 1;
	PROTECT(ans = NEW_INTEGER(ans_len));
	memset(INTEGER(ans), 0, LENGTH(ans) * sizeof(int));
	nb_valid_left_char = 0;
	for (i = 0; i < x_len; i++, x++) {
		twobit = eightbit2twobit_lkup[(unsigned char) *x];
		if (twobit == -1) {
			nb_valid_left_char = 0;
		} else {
			nb_valid_left_char++;
			ans_offset <<= 2;
			ans_offset &= ans_offset_bitmask;
			ans_offset += twobit;
			if (nb_valid_left_char >= INTEGER(width)[0])
				INTEGER(ans)[ans_offset]++;
		}
	}
	UNPROTECT(1);
	return ans;
}

