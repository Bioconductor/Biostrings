#include "Biostrings.h"
#include <string.h>


static int code2offset[CHRTRTABLE_LENGTH];

static void init_code2offset(SEXP codes, int with_other, int nrow)
{
	int i;

	_init_chrtrtable(INTEGER(codes), LENGTH(codes), code2offset);
	for (i = 0; i < CHRTRTABLE_LENGTH; i++) {
		if (code2offset[i] == -1) {
			if (!with_other)
				continue;
			code2offset[i] = LENGTH(codes);
		}
		code2offset[i] *= nrow;
	}
	return;
}

static void add_freqs(RoSeq X, SEXP codes, int *freqs, int nrow)
{
	int i, offset;

	if (codes == R_NilValue)
		for (i = 0; i < X.nelt; i++, X.elts++)
			freqs[((unsigned char) *X.elts) * nrow]++;
	else 
		for (i = 0; i < X.nelt; i++, X.elts++) {
			offset = code2offset[(unsigned char) *X.elts];
			if (offset == -1) 
				continue;
			freqs[offset]++;
		}
	return;
}

static SEXP append_other_to_names(SEXP codes)
{
	SEXP names, name, codes_names;
	int i;

	PROTECT(names = NEW_CHARACTER(LENGTH(codes) + 1));
	codes_names = GET_NAMES(codes);
	for (i = 0; i < LENGTH(codes); i++) {
		if (codes_names == R_NilValue)
			PROTECT(name = mkChar(""));
		else
			PROTECT(name = duplicate(STRING_ELT(codes_names, i)));
		SET_STRING_ELT(names, i, name);
		UNPROTECT(1);
	}
	SET_STRING_ELT(names, i, mkChar("other"));
	UNPROTECT(1);
	return names;
}

static void set_names(SEXP x, SEXP codes, int with_other, int collapse)
{
	SEXP names, codes_names, dim_names;

	if (codes == R_NilValue)
		return;
	if (with_other) {
		PROTECT(names = append_other_to_names(codes));
	} else {
		codes_names = GET_NAMES(codes);
		PROTECT(names = duplicate(codes_names));
	}
	if (collapse) {
		SET_NAMES(x, names);
	} else {
		PROTECT(dim_names = NEW_LIST(2));
		SET_ELEMENT(dim_names, 0, R_NilValue);
		SET_ELEMENT(dim_names, 1, names);
		SET_DIMNAMES(x, dim_names);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	return;
}

/*
 * --- .Call ENTRY POINT ---
 */
SEXP XString_char_frequency(SEXP x, SEXP codes, SEXP with_other)
{
	int ans_length, *freqs;
	SEXP ans;
	RoSeq X;

	if (codes == R_NilValue) {
		ans_length = 256;
	} else {
		ans_length = LENGTH(codes);
		if (LOGICAL(with_other)[0])
			ans_length++;
		init_code2offset(codes, LOGICAL(with_other)[0], 1);
	}
	PROTECT(ans = NEW_INTEGER(ans_length));
	freqs = INTEGER(ans);
	memset(freqs, 0, ans_length * sizeof(int));
	X = _get_XString_asRoSeq(x);
	add_freqs(X, codes, freqs, 1);
	set_names(ans, codes, LOGICAL(with_other)[0], 1);
	UNPROTECT(1);
	return ans;
}

/*
 * --- .Call ENTRY POINT ---
 */
SEXP XStringSet_char_frequency(SEXP x, SEXP codes, SEXP with_other,
		SEXP collapse)
{
	int x_length, ans_nrow, ans_ncol, ans_length, *freqs, i;
	SEXP ans;
	RoSeq X;

	x_length = _get_XStringSet_length(x);
	ans_nrow = LOGICAL(collapse)[0] ? 1 : x_length;
	if (codes == R_NilValue) {
		ans_ncol = 256;
	} else {
		ans_ncol = LENGTH(codes);
		if (LOGICAL(with_other)[0])
			ans_ncol++;
		init_code2offset(codes, LOGICAL(with_other)[0], ans_nrow);
	}
	ans_length = ans_nrow * ans_ncol;
	if (LOGICAL(collapse)[0])
		PROTECT(ans = NEW_INTEGER(ans_length));
	else
		PROTECT(ans = allocMatrix(INTSXP, ans_nrow, ans_ncol));
	freqs = INTEGER(ans);
	memset(freqs, 0, ans_length * sizeof(int));
	for (i = 0; i < x_length; i++) {
		X = _get_XStringSet_elt_asRoSeq(x, i);
		add_freqs(X, codes, freqs, ans_nrow);
		if (!LOGICAL(collapse)[0])
			freqs++;
	}
	set_names(ans, codes, LOGICAL(with_other)[0], LOGICAL(collapse)[0]);
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

