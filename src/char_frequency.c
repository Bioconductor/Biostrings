#include "Biostrings.h"

static int code2offset[CHRTRTABLE_LENGTH];
static int rowbuf[256];

static int get_ans_width(SEXP codes, int with_other)
{
	int width, i;

	if (codes == R_NilValue)
		return 256;
	width = LENGTH(codes);
	_init_chrtrtable(INTEGER(codes), width, code2offset);
	if (with_other) {
		for (i = 0; i < CHRTRTABLE_LENGTH; i++)
			if (code2offset[i] == -1)
				code2offset[i] = width;
		width++;
	}
	return width;
}

static void add_freqs(RoSeq X, SEXP codes, int *freqs)
{
	static int i, offset;

	if (codes == R_NilValue)
		for (i = 0; i < X.nelt; i++, X.elts++)
			freqs[((unsigned char) *X.elts)]++;
	else
		for (i = 0; i < X.nelt; i++, X.elts++) {
			offset = code2offset[(unsigned char) *X.elts];
			if (offset == -1)
				continue;
			freqs[offset]++;
		}
	return;
}

static void add_freqs_by_pos(RoSeq X, SEXP codes, int nrow, int *freqs)
{
	static int i, offset;

	if (codes == R_NilValue)
		for (i = 0; i < X.nelt; i++, X.elts++)
			freqs[(i * nrow) + ((unsigned char) *X.elts)]++;
	else
		for (i = 0; i < X.nelt; i++, X.elts++) {
			offset = code2offset[(unsigned char) *X.elts];
			if (offset == -1)
				continue;
			freqs[(i * nrow) + offset]++;
		}
	return;
}

static void copy_rowbuf_to_row0_in_matrix(int *matrix, int nrow, int ncol)
{
	static int i;

	for (i = 0; i < ncol; i++)
		matrix[i * nrow] = rowbuf[i];
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

static void set_names(SEXP x, SEXP codes, int with_other, int collapse, int which_names)
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
		SET_ELEMENT(dim_names, 1 - which_names, R_NilValue);
		SET_ELEMENT(dim_names, which_names, names);
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
	SEXP ans;
	int ans_length, *freqs;
	RoSeq X;

	ans_length = get_ans_width(codes, LOGICAL(with_other)[0]);
	PROTECT(ans = NEW_INTEGER(ans_length));
	freqs = INTEGER(ans);
	memset(freqs, 0, ans_length * sizeof(int));
	X = _get_XString_asRoSeq(x);
	add_freqs(X, codes, freqs);
	set_names(ans, codes, LOGICAL(with_other)[0], 1, 1);
	UNPROTECT(1);
	return ans;
}

/*
 * --- .Call ENTRY POINT ---
 */
SEXP XStringSet_char_frequency(SEXP x, SEXP codes, SEXP with_other,
		SEXP collapse)
{
	SEXP ans;
	int ans_width, x_length, *freqs, i;
	CachedXStringSet cached_x;
	RoSeq xx;

	ans_width = get_ans_width(codes, LOGICAL(with_other)[0]);
	x_length = _get_XStringSet_length(x);
	cached_x = _new_CachedXStringSet(x);
	if (LOGICAL(collapse)[0]) {
		PROTECT(ans = NEW_INTEGER(ans_width));
		freqs = INTEGER(ans);
		memset(freqs, 0, ans_width * sizeof(int));
		for (i = 0; i < x_length; i++) {
			xx = _get_CachedXStringSet_elt_asRoSeq(&cached_x, i);
			add_freqs(xx, codes, freqs);
		}
	} else {
		PROTECT(ans = allocMatrix(INTSXP, x_length, ans_width));
		for (i = 0, freqs = INTEGER(ans); i < x_length; i++, freqs++) {
			xx = _get_CachedXStringSet_elt_asRoSeq(&cached_x, i);
			memset(rowbuf, 0, ans_width * sizeof(int));
			add_freqs(xx, codes, rowbuf);
			copy_rowbuf_to_row0_in_matrix(freqs, x_length, ans_width);
		}
	}
	set_names(ans, codes, LOGICAL(with_other)[0], LOGICAL(collapse)[0], 1);
	UNPROTECT(1);
	return ans;
}

/*
 * --- .Call ENTRY POINT ---
 */
SEXP XStringSet_char_frequency_by_pos(SEXP x, SEXP codes, SEXP with_other)
{
	SEXP ans;
	int ans_nrow, ans_ncol, ans_width, x_length, *freqs, i;
	CachedXStringSet cached_x;
	RoSeq xx;

	x_length = _get_XStringSet_length(x);
	cached_x = _new_CachedXStringSet(x);
	ans_nrow = get_ans_width(codes, LOGICAL(with_other)[0]);
	ans_ncol = 0;
	for (i = 0; i < x_length; i++) {
		xx = _get_CachedXStringSet_elt_asRoSeq(&cached_x, i);
		if (ans_ncol < xx.nelt)
			ans_ncol = xx.nelt;
	}
	ans_width = ans_nrow * ans_ncol;
	PROTECT(ans = allocMatrix(INTSXP, ans_nrow, ans_ncol));
	freqs = INTEGER(ans);
	memset(freqs, 0, ans_width * sizeof(int));
	for (i = 0; i < x_length; i++) {
		xx = _get_CachedXStringSet_elt_asRoSeq(&cached_x, i);
		add_freqs_by_pos(xx, codes, ans_nrow, freqs);
	}
	set_names(ans, codes, LOGICAL(with_other)[0], 0, 0);
	UNPROTECT(1);
	return ans;
}

/*
 * --- .Call ENTRY POINT ---
 */
SEXP oligonucleotide_frequency(SEXP x, SEXP base_codes, SEXP width,
		SEXP fast_moving_side)
{
	RoSeq X;
	SEXP ans;
	int ans_len, ans_offset_bitmask, ans_offset, width0, nbit_in_mask,
            right_moves_fastest, i, nb_valid_left_char, twobit;
	const char *c;

	static int eightbit2twobit_lkup[256];

	X = _get_XString_asRoSeq(x);
	if (LENGTH(base_codes) != 4)
		error("'base_codes' must be of length 4");
	_init_chrtrtable(INTEGER(base_codes), LENGTH(base_codes),
			 eightbit2twobit_lkup);
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
	for (i = 0, c = X.elts; i < X.nelt; i++, c++) {
		twobit = eightbit2twobit_lkup[(unsigned char) *c];
		if (twobit == -1) {
			nb_valid_left_char = 0;
			continue;
		}
		nb_valid_left_char++;
		ans_offset &= ans_offset_bitmask;
		if (right_moves_fastest) {
			ans_offset <<= 2;
		} else {
			ans_offset >>= 2;
			twobit <<= nbit_in_mask;
		}
		ans_offset += twobit;
		if (nb_valid_left_char < width0)
			continue;
		INTEGER(ans)[ans_offset]++;
	}
	UNPROTECT(1);
	return ans;
}

