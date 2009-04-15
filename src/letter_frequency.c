#include "Biostrings.h"
#include <stdlib.h> /* for abs() */

static ByteTrTable byte2offset;
static int rowbuf[256];

static int get_ans_width(SEXP codes, int with_other)
{
	int width, i;

	if (codes == R_NilValue)
		return 256;
	_init_byte2offset_with_INTEGER(byte2offset, codes, 1);
	width = LENGTH(codes);
	if (with_other) {
		for (i = 0; i < BYTETRTABLE_LENGTH; i++)
			if (byte2offset[i] == NA_INTEGER)
				byte2offset[i] = width;
		width++;
	}
	return width;
}

static void update_freqs(int *freqs, const RoSeq *X, SEXP codes)
{
	const char *x;
	int i, offset;

	x = X->elts;
	if (codes == R_NilValue)
		for (i = 0; i < X->nelt; i++, x++)
			freqs[((unsigned char) *x)]++;
	else
		for (i = 0; i < X->nelt; i++, x++) {
			offset = byte2offset[(unsigned char) *x];
			if (offset == NA_INTEGER)
				continue;
			freqs[offset]++;
		}
	return;
}

/* Note that calling update_freqs2() with shift = 0, nrow = 0 and
   ncol = X->nelt is equivalent to calling update_freqs() */
static void update_freqs2(int *freqs, const RoSeq *X, SEXP codes,
		int shift, int nrow, int ncol)
{
	const char *x;
	int i1, i2, j1, j2, i, offset;

	if (abs(shift) >= ncol)
		return;
	/* i1, i2 are 0-based indices in X->elts
	   (range i1 <= i < i2 must be safe) */
	i1 = 0;
	i2 = X->nelt;
	/* j1, j2 are 0-based column indices in the freqs matrix
	   (range j1 <= j < j2 must be safe) */
	j1 = i1 + shift;
	j2 = i2 + shift;
	if (j1 < 0) {
		i1 -= j1;
		j1 = 0;
	}
	if (j2 > ncol) {
		i2 -= j2 - ncol;
		/* j2 = ncol; not needed */
	}
	x = X->elts + i1;
	freqs += j1 * nrow;
	if (codes == R_NilValue)
		for (i = i1; i < i2; i++, x++, freqs += nrow) {
			freqs[(unsigned char) *x]++;
		}
	else
		for (i = i1; i < i2; i++, x++, freqs += nrow) {
			offset = byte2offset[(unsigned char) *x];
			if (offset == NA_INTEGER)
				continue;
			freqs[offset]++;
		}
	return;
}

static void update_oligo_freqs(int *freqs, int freqs_nrow, int width,
		int invert_twobit_order, const RoSeq *X)
{
	int nbit_in_mask, twobit_mask, nb_valid_left_char,
	    i, twobit_c, twobit_oligo;
	const char *c;

	nbit_in_mask = (width - 1) * 2;
	twobit_mask = (1 << nbit_in_mask) - 1;
	if (invert_twobit_order)
		twobit_mask <<= 2;
	nb_valid_left_char = 0;
	for (i = 0, c = X->elts; i < X->nelt; i++, c++) {
		twobit_c = byte2offset[(unsigned char) *c];
		if (twobit_c == NA_INTEGER) {
			nb_valid_left_char = 0;
			continue;
		}
		nb_valid_left_char++;
		twobit_oligo &= twobit_mask;
		if (invert_twobit_order) {
			twobit_oligo >>= 2;
			twobit_c <<= nbit_in_mask;
		} else {
			twobit_oligo <<= 2;
		}
		twobit_oligo += twobit_c;
		if (nb_valid_left_char < width)
			continue;
		freqs[twobit_oligo * freqs_nrow]++;
	}
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
SEXP XString_letter_frequency(SEXP x, SEXP codes, SEXP with_other)
{
	SEXP ans;
	int ans_length, *freqs;
	RoSeq X;

	ans_length = get_ans_width(codes, LOGICAL(with_other)[0]);
	PROTECT(ans = NEW_INTEGER(ans_length));
	freqs = INTEGER(ans);
	memset(freqs, 0, ans_length * sizeof(int));
	X = _get_XString_asRoSeq(x);
	update_freqs(freqs, &X, codes);
	set_names(ans, codes, LOGICAL(with_other)[0], 1, 1);
	UNPROTECT(1);
	return ans;
}

/*
 * --- .Call ENTRY POINT ---
 */
SEXP XStringSet_letter_frequency(SEXP x, SEXP codes, SEXP with_other,
		SEXP collapse)
{
	SEXP ans;
	int ans_width, x_length, *freqs, i;
	CachedXStringSet cached_x;
	RoSeq X;

	ans_width = get_ans_width(codes, LOGICAL(with_other)[0]);
	x_length = _get_XStringSet_length(x);
	cached_x = _new_CachedXStringSet(x);
	if (LOGICAL(collapse)[0]) {
		PROTECT(ans = NEW_INTEGER(ans_width));
		freqs = INTEGER(ans);
		memset(freqs, 0, ans_width * sizeof(int));
		for (i = 0; i < x_length; i++) {
			X = _get_CachedXStringSet_elt_asRoSeq(&cached_x, i);
			update_freqs(freqs, &X, codes);
		}
	} else {
		PROTECT(ans = allocMatrix(INTSXP, x_length, ans_width));
		for (i = 0, freqs = INTEGER(ans); i < x_length; i++, freqs++) {
			X = _get_CachedXStringSet_elt_asRoSeq(&cached_x, i);
			memset(rowbuf, 0, ans_width * sizeof(int));
			update_freqs(rowbuf, &X, codes);
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
SEXP XStringSet_letter_frequency_by_pos(SEXP x, SEXP codes, SEXP with_other,
		SEXP shift, SEXP width)
{
	SEXP ans;
	int ans_nrow, ans_ncol, ans_length, x_length, *freqs, i, k, s, X_end;
	CachedXStringSet cached_x;
	RoSeq X;

	ans_nrow = get_ans_width(codes, LOGICAL(with_other)[0]);
	x_length = _get_XStringSet_length(x);
	cached_x = _new_CachedXStringSet(x);
	if (width == R_NilValue) {
		if (x_length == 0)
			error("'x' has no element and 'width' is NULL");
		if (LENGTH(shift) == 0)
			error("'shift' has no element");
		ans_ncol = 0;
		for (i = k = 0; i < x_length; i++, k++) {
			if (k >= LENGTH(shift))
				k = 0; /* recycle */
			s = INTEGER(shift)[k];
			if (s == NA_INTEGER)
				error("'shift' contains NAs");
			X = _get_CachedXStringSet_elt_asRoSeq(&cached_x, i);
			X_end = X.nelt + s;
			if (X_end > ans_ncol)
				ans_ncol = X_end;
		}
	} else {
		if (x_length != 0 && LENGTH(shift) == 0)
			error("'shift' has no element");
		ans_ncol = INTEGER(width)[0];
	}
	ans_length = ans_nrow * ans_ncol;
	PROTECT(ans = allocMatrix(INTSXP, ans_nrow, ans_ncol));
	freqs = INTEGER(ans);
	memset(freqs, 0, ans_length * sizeof(int));
	for (i = k = 0; i < x_length; i++, k++) {
		if (k >= LENGTH(shift))
			k = 0; /* recycle */
		s = INTEGER(shift)[k];
		if (s == NA_INTEGER)
			error("'shift' contains NAs");
		X = _get_CachedXStringSet_elt_asRoSeq(&cached_x, i);
		update_freqs2(freqs, &X, codes, s, ans_nrow, ans_ncol);
	}
	set_names(ans, codes, LOGICAL(with_other)[0], 0, 0);
	UNPROTECT(1);
	return ans;
}

/*
 * --- .Call ENTRY POINT ---
 */
SEXP XString_oligonucleotide_frequency(SEXP x, SEXP base_codes, SEXP width,
		SEXP fast_moving_side)
{
	RoSeq X;
	SEXP ans;
	int ans_len, ans_offset_bitmask, ans_offset, width0, nbit_in_mask,
	    right_moves_fastest, i, nb_valid_left_char, twobit;
	const char *c;

	static ByteTrTable eightbit2twobit;

	if (LENGTH(base_codes) != 4)
		error("'base_codes' must be of length 4");
	X = _get_XString_asRoSeq(x);
	_init_byte2offset_with_INTEGER(eightbit2twobit, base_codes, 1);
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
		twobit = eightbit2twobit[(unsigned char) *c];
		if (twobit == NA_INTEGER) {
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

SEXP XStringSet_oligonucleotide_frequency(SEXP x, SEXP base_codes, SEXP width,
		SEXP fast_moving_side, SEXP collapse)
{
	SEXP ans;
	int width0, invert_twobit_order, ans_width, x_length, *freqs, i;
	CachedXStringSet cached_x;
	RoSeq x_elt;

	if (LENGTH(base_codes) != 4)
		error("'base_codes' must be of length 4");
	_init_byte2offset_with_INTEGER(byte2offset, base_codes, 1);
	width0 = INTEGER(width)[0];
	if (width0 < 1 || width0 > 12)
		error("'width' must be >=1 and <= 12");
	invert_twobit_order = strcmp(CHAR(STRING_ELT(fast_moving_side, 0)), "right") != 0;
	ans_width = 1 << (width0 * 2);
	x_length = _get_XStringSet_length(x);
	cached_x = _new_CachedXStringSet(x);
	if (LOGICAL(collapse)[0]) {
		PROTECT(ans = NEW_INTEGER(ans_width));
		freqs = INTEGER(ans);
		memset(freqs, 0, LENGTH(ans) * sizeof(int));
		for (i = 0; i < x_length; i++) {
			x_elt = _get_CachedXStringSet_elt_asRoSeq(&cached_x, i);
			update_oligo_freqs(freqs, 1, width0,
						invert_twobit_order, &x_elt);
		}
	} else {
		PROTECT(ans = allocMatrix(INTSXP, x_length, ans_width));
		freqs = INTEGER(ans);
		memset(freqs, 0, LENGTH(ans) * sizeof(int));
		for (i = 0; i < x_length; i++, freqs++) {
			x_elt = _get_CachedXStringSet_elt_asRoSeq(&cached_x, i);
			update_oligo_freqs(freqs, x_length, width0,
						invert_twobit_order, &x_elt);
		}
	}
	UNPROTECT(1);
	return ans;
}

