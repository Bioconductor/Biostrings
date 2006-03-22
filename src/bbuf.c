#include "Biostrings.h"

/*
 * From R:
 *   .Call("sexp_address", 6:4, PACKAGE="Biostrings")
 *   .Call("sexp_address", new("externalptr"), PACKAGE="Biostrings")
 */
SEXP sexp_address(SEXP s)
{
	SEXP string, ans;
	char addr[40]; /* should be enough, even for 128-bit addresses */

	snprintf(addr, sizeof(addr), "%p", s);
	PROTECT(string = allocString(strlen(addr)));
	strcpy(CHAR(string), addr);
	PROTECT(ans = allocVector(STRSXP, 1));
	SET_STRING_ELT(ans, 0, string);
	UNPROTECT(2);
	return ans;
}

/*
 * Print some obscure info about an "externalptr" object.
 * From R:
 *   .Call("xp_show", new("externalptr"), PACKAGE="Biostrings")
 */
SEXP xp_show(SEXP xp)
{
	SEXP s;
	void *p;

	Rprintf("Object of class 'externalptr':\n");
	Rprintf("  xp adress: %p\n", xp);
	p = R_ExternalPtrAddr(xp);
	Rprintf("  R_ExternalPtrAddr(xp): %p\n", p);
	s = R_ExternalPtrTag(xp);
	Rprintf("  R_ExternalPtrTag(xp): %p", s);
	Rprintf("%s\n", TYPEOF(s) == NILSXP ? " (NILSXP)" : "");
	s = R_ExternalPtrProtected(xp);
	Rprintf("  R_ExternalPtrProtected(xp): %p", s);
	Rprintf("%s\n", TYPEOF(s) == NILSXP ? " (NILSXP)" : "");
	return R_NilValue;
}

/*
 * new("externalptr") will always return the same instance of an external
 * pointer object! If you need a new instance, use this function instead.
 * From R:
 *   xp <- .Call("xp_new", PACKAGE="Biostrings")
 */
SEXP xp_new()
{
	return R_MakeExternalPtr(NULL, R_NilValue, R_NilValue);
}

/*
 * We can't rely on the strsplit() R function to split a string into single
 * characters when the string contains junk. For example:
 *   bb <- bbuf(5)
 *   bb[] <- 255
 *   s <- toString(bb)
 *   strsplit(s, NULL)
 * doesn't work!
 * The function below should be safe, whatever the content of 's' is!
 * The length of the returned string is the number of chars in single
 * string s. Not vectorized.
 */
SEXP safe_explode(SEXP s)
{
	SEXP s0, string, ans;
	int s0_length, i;

	s0 = STRING_ELT(s, 0);
	s0_length = LENGTH(s0);

	PROTECT(ans = allocVector(STRSXP, s0_length));
	for (i = 0; i < s0_length; i++) {
		PROTECT(string = allocString(1));
		CHAR(string)[0] = CHAR(s0)[i];
		SET_STRING_ELT(ans, i, string);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * Memory allocation for a "bbuf" object
 * -------------------------------------
 * A "bbuf" object stores its data in an "external" R string (CHARSXP vector).
 * An R string itself stores its data in a char-array.
 * The "R types" of the argument passed to these functions must be:
 *   'bb_xp': externalptr
 *   'length': single integer
 */

/*
 * Alloc an R string of length 'length' and point 'bb_xp' to it.
 * Does NOT initialize the allocated memory!
 */
SEXP bbuf_alloc(SEXP bb_xp, SEXP length)
{
	SEXP tag;
	int tag_length;

	tag_length = INTEGER(length)[0];

	PROTECT(tag = allocString(tag_length));
	/*
	Rprintf("Memory successfully allocated for %d-byte buffer (starting at address %p)\n",
		tag_length, CHAR(tag));
	 */
	R_SetExternalPtrTag(bb_xp, tag);
	UNPROTECT(1);
	return bb_xp;
}

/*
 * Print some info about the R string pointed by 'bb_xp'.
 * From R:
 *   bb <- bbuf(30)
 *   .Call("bbuf_show", bb@xp, PACKAGE="Biostrings")
 */
SEXP bbuf_show(SEXP bb_xp)
{
	SEXP tag;
	int tag_length;

	tag = R_ExternalPtrTag(bb_xp);
	tag_length = LENGTH(tag);
	Rprintf("%d-byte buffer (starting at address %p)\n",
		tag_length, CHAR(tag));
	return R_NilValue;
}

/*
 * Return length of R string pointed by 'bb_xp'.
 * From R:
 *   bb <- bbuf(30)
 *   .Call("bbuf_length", bb@xp, PACKAGE="Biostrings")
 * Called by method length() for "bbuf" objects.
 */
SEXP bbuf_length(SEXP bb_xp)
{
	SEXP tag, ans;
	int tag_length;

	tag = R_ExternalPtrTag(bb_xp);
	tag_length = LENGTH(tag);

	PROTECT(ans = allocVector(INTSXP, 1));
	INTEGER(ans)[0] = tag_length;
	UNPROTECT(1);
	return ans;
}

/*
 * From R:
 *   bb <- bbuf(30)
 *   .Call("bbuf_memcmp", bb@xp, 1:1, bb@xp, 10:10, 21:21, PACKAGE="Biostrings")
 */
SEXP bbuf_memcmp(SEXP bb1_xp, SEXP first1,
		 SEXP bb2_xp, SEXP first2, SEXP width)
{
	SEXP tag1, tag2, ans;
	int i1, i2, n;

	tag1 = R_ExternalPtrTag(bb1_xp);
	i1 = INTEGER(first1)[0] - 1;
	tag2 = R_ExternalPtrTag(bb2_xp);
	i2 = INTEGER(first2)[0] - 1;
	n = INTEGER(width)[0];

	PROTECT(ans = allocVector(INTSXP, 1));
	INTEGER(ans)[0] = _memcmp(CHAR(tag1), i1,
				  CHAR(tag2), i2,
				  n, sizeof(char));
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * READ/WRITE functions
 * ====================
 * The functions in this section implement the read/write operations to a
 * "bbuf" object. The user can choose between 2 interfaces for each read
 * or write operation:
 *
 *   1. The "imin/imax" interface: the bytes to access are specified by 2
 * integers: 'imin' (the position of the first byte to access, the first
 * byte in the buffer being at position 1) and 'imax' (the position of the
 * last byte to access).
 *
 *   2. The "ii" interface: the bytes to access are specified by an
 * integer vector containing their positions in the buffer.
 *
 * The "ii" interface is intended to be used by the subsetting
 * operator [ defined at the R level for "bbuf" objects.
 * R subsetting operator [ can be used to read values from, or write values
 * to an object that contains a collection of values, like a character
 * vector, an integer vector or a logical vector.
 * If x is a vector and i an integer vector of length n with the following
 * properties:
 *   a) i contains no NA values,
 *   b) i can be used to subset x without being "out of bounds" (i.e all
 *      values in i are >= 1 and <= length(x)),
 * then we have the following properties:
 *   1) READING from x: y <- x[i] produces a vector, of the same type than x,
 *      but of the same length than i (length(y) == n).
 *   2) READING from then WRITING to x: x[i] <- x[i] (short for y <- x[i];
 *      x[i] <- y) doesn't modify the values in x.
 *   3) WRITING to then READING from x: if z is a vector of length n and of
 *      the same type than x, then doing x[i] <- z; y <- x[i] guarantees that
 *      y is identical to z only when i contains no repeated value!
 *
 * Functions in this section that implement the "ii" interface
 * respect the above properties.
 *
 * Here are some arguments to these functions that must always be SEXP of the
 * following types:
 *   bb_xp: externalptr
 *   imin, imax: single integers
 *   ii: integer vector containing the subscripts (with no NAs)
 *   enc_xp: externalptr (hash table for encoding)
 *   dec_xp: externalptr (hash table for decoding)
 */


/* ==========================================================================
 * Copy bytes from a "bbuf" to another "bbuf" object.
 * --------------------------------------------------------------------------
 */

SEXP bbuf_copy(SEXP dest_xp, SEXP imin, SEXP imax, SEXP src_xp)
{
	SEXP dest_tag, src_tag;
	int i1, i2;

	dest_tag = R_ExternalPtrTag(dest_xp);
	i1 = INTEGER(imin)[0] - 1;
	i2 = INTEGER(imax)[0] - 1;
	src_tag = R_ExternalPtrTag(src_xp);

	_memcpy_from_range(i1, i2,
			CHAR(dest_tag), LENGTH(dest_tag),
			CHAR(src_tag), LENGTH(src_tag), sizeof(char));
	return dest_xp;
}

SEXP bbuf_copyii(SEXP dest_xp, SEXP ii, SEXP src_xp)
{
	SEXP dest_tag, src_tag;

	dest_tag = R_ExternalPtrTag(dest_xp);
	src_tag = R_ExternalPtrTag(src_xp);

	_memcpy_from_subset(INTEGER(ii), LENGTH(ii),
			CHAR(dest_tag), LENGTH(dest_tag),
			CHAR(src_tag), LENGTH(src_tag), sizeof(char));
	return dest_xp;
}


/* ==========================================================================
 * Read/write chars to a "bbuf" object
 * --------------------------------------------------------------------------
 */

/*
 * Return a single string (character vector of length 1).
 * From R:
 *   bb <- bbuf(15)
 *   bb[] < "Hello"
 *   .Call("bbuf_read_chars", bb@xp, 2:2, 4:4, PACKAGE="Biostrings")
 */
SEXP bbuf_read_chars(SEXP bb_xp, SEXP imin, SEXP imax)
{
	SEXP tag, string, ans;
	int i1, i2, n;

	tag = R_ExternalPtrTag(bb_xp);
	i1 = INTEGER(imin)[0] - 1;
	i2 = INTEGER(imax)[0] - 1;
	n = i2 - i1 + 1;
	PROTECT(string = allocString(n));
	_memcpy_from_range(i1, i2,
			CHAR(string), n,
			CHAR(tag), LENGTH(tag), sizeof(char));
	PROTECT(ans = allocVector(STRSXP, 1));
	SET_STRING_ELT(ans, 0, string);
	UNPROTECT(2);
	return ans;
}

SEXP bbuf_readii_chars(SEXP bb_xp, SEXP ii)
{
	SEXP tag, string, ans;
	int n;

	tag = R_ExternalPtrTag(bb_xp);
	n = LENGTH(ii);

	PROTECT(string = allocString(n));
	_memcpy_from_subset(INTEGER(ii), n,
			CHAR(string), n,
			CHAR(tag), LENGTH(tag), sizeof(char));
	PROTECT(ans = allocVector(STRSXP, 1));
	SET_STRING_ELT(ans, 0, string);
	UNPROTECT(2);
	return ans;
}

/*
 * 'val' must be a non-empty single string (character vector of length 1).
 */
SEXP bbuf_write_chars(SEXP bb_xp, SEXP imin, SEXP imax, SEXP val)
{
	SEXP tag, string;
	int i1, i2;

	tag = R_ExternalPtrTag(bb_xp);
	string = STRING_ELT(val, 0);
	i1 = INTEGER(imin)[0] - 1;
	i2 = INTEGER(imax)[0] - 1;
	_memcpy_to_range(i1, i2,
			CHAR(tag), LENGTH(tag),
			CHAR(string), LENGTH(string), sizeof(char));
	return bb_xp;
}

SEXP bbuf_writeii_chars(SEXP bb_xp, SEXP ii, SEXP val)
{
	SEXP tag, string;

	tag = R_ExternalPtrTag(bb_xp);
	string = STRING_ELT(val, 0);
	_memcpy_to_subset(INTEGER(ii), LENGTH(ii),
			CHAR(tag), LENGTH(tag),
			CHAR(string), LENGTH(string), sizeof(char));
	return bb_xp;
}


/* ==========================================================================
 * Read/write integers to a "bbuf" object
 * --------------------------------------------------------------------------
 */

/*
 * Return an integer vector of length 'imax' - 'imin' + 1.
 * From R:
 *   bb <- bbuf(30)
 *   .Call("bbuf_read_ints", bb@xp, 20:20, 25:25, PACKAGE="Biostrings")
 */
SEXP bbuf_read_ints(SEXP bb_xp, SEXP imin, SEXP imax)
{
	SEXP tag, ans;
	int i1, i2, n, j;

	tag = R_ExternalPtrTag(bb_xp);
	i1 = INTEGER(imin)[0] - 1;
	i2 = INTEGER(imax)[0] - 1;
	if (i1 < 0 || i2 >= LENGTH(tag))
		error("subscript out of bounds");
	n = i2 - i1 + 1;

	PROTECT(ans = allocVector(INTSXP, n));
	for (j = 0; i1 <= i2; i1++, j++) {
		INTEGER(ans)[j] = (unsigned char) CHAR(tag)[i1];
	}
	UNPROTECT(1);
	return ans;
}

/*
 * Return an integer vector of same length than 'ii'.
 * From R:
 *   bb <- bbuf(30)
 *   .Call("bbuf_readii_ints", bb, 25:20, PACKAGE="Biostrings")
 */
SEXP bbuf_readii_ints(SEXP bb_xp, SEXP ii)
{
	SEXP tag, ans;
	int tag_length;
	int n, i, j;

	tag = R_ExternalPtrTag(bb_xp);
	tag_length = LENGTH(tag);
	n = LENGTH(ii);

	PROTECT(ans = allocVector(INTSXP, n));
	for (j = 0; j < n; j++) {
		i = INTEGER(ii)[j] - 1;
		if (i < 0 || i >= tag_length)
			error("subscript out of bounds");
		INTEGER(ans)[j] = (unsigned char) CHAR(tag)[i];
	}
	UNPROTECT(1);
	return ans;
}

/*
 * 'val' must be an integer vector of length > 0.
 */
SEXP bbuf_write_ints(SEXP bb_xp, SEXP imin, SEXP imax, SEXP val)
{
	SEXP tag;
	int val_length;
	int i1, i2, n, j;
	int v;

	tag = R_ExternalPtrTag(bb_xp);
	i1 = INTEGER(imin)[0] - 1;
	i2 = INTEGER(imax)[0] - 1;
	if (i1 < 0 || i2 >= LENGTH(tag))
		error("subscript out of bounds");
	n = i2 - i1 + 1;
	val_length = LENGTH(val);
	if (val_length == 0 && n != 0)
		error("no value provided");

	for (j = 0; i1 <= i2; i1++, j++) {
		if (j >= val_length)
			j = 0; /* recycle */
		v = INTEGER(val)[j];
		if (v < 0 || v >= 256)
			error("value out of range");
		CHAR(tag)[i1] = (char) v;
	}
	if (j != val_length) {
		warning("number of items to replace is not a multiple "
			"of replacement length");
	}
	return bb_xp;
}

SEXP bbuf_writeii_ints(SEXP bb_xp, SEXP ii, SEXP val)
{
	SEXP tag;
	int tag_length, val_length;
	int n, i, j, z;
	int v;

	val_length = LENGTH(val);
	n = LENGTH(ii);
	if (val_length == 0 && n != 0)
		error("no value provided");
	tag = R_ExternalPtrTag(bb_xp);
	tag_length = LENGTH(tag);

	for (j = z = 0; z < n; j++, z++) {
		i = INTEGER(ii)[z] - 1;
		if (i < 0 || i >= tag_length)
			error("subscript out of bounds");
		if (j >= val_length)
			j = 0; /* recycle */
		v = INTEGER(val)[j];
		if (v < 0 || v >= 256)
			error("value out of range");
		CHAR(tag)[i] = (char) v;
	}
	if (j != val_length) {
		warning("number of items to replace is not a multiple "
			"of replacement length");
	}
	return bb_xp;
}


/* ==========================================================================
 * Read/write encoded chars to a "bbuf" object
 * --------------------------------------------------------------------------
 */

/*
 * Return a single string (character vector of length 1).
 */
SEXP bbuf_read_enc_chars(SEXP bb_xp, SEXP imin, SEXP imax, SEXP dec_xp)
{
	SEXP tag, dectag, string, ans;
	int dectag_length;
	int i1, i2, n, j, h, code;
	char dec_hole, letter;

	tag = R_ExternalPtrTag(bb_xp);
	i1 = INTEGER(imin)[0] - 1;
	i2 = INTEGER(imax)[0] - 1;
	if (i1 < 0 || i2 >= LENGTH(tag))
		error("subscript out of bounds");
	n = i2 - i1 + 1;
	dectag = R_ExternalPtrTag(dec_xp);
	dectag_length = LENGTH(dectag);
	dec_hole = CHAR(dectag)[0];

	PROTECT(string = allocString(n));
	for (j = 0; i1 <= i2; i1++, j++) {
		code = (unsigned char) CHAR(tag)[i1];
		h = code + 1;
		if (h >= dectag_length
		    || (letter = CHAR(dectag)[h]) == dec_hole)
			error("unknown code %d in string to decode", code);
		CHAR(string)[j] = letter;
	}
	PROTECT(ans = allocVector(STRSXP, 1));
	SET_STRING_ELT(ans, 0, string);
	UNPROTECT(2);
	return ans;
}

SEXP bbuf_readii_enc_chars(SEXP bb_xp, SEXP ii, SEXP dec_xp)
{
	SEXP tag, dectag, string, ans;
	int tag_length, dectag_length;
	int n, i, j, h, code;
	char dec_hole, letter;

	n = LENGTH(ii);
	dectag = R_ExternalPtrTag(dec_xp);
	dectag_length = LENGTH(dectag);
	dec_hole = CHAR(dectag)[0];
	tag = R_ExternalPtrTag(bb_xp);
	tag_length = LENGTH(tag);

	PROTECT(string = allocString(n));
	for (j = 0; j < n; j++) {
		i = INTEGER(ii)[j] - 1;
		if (i < 0 || i >= tag_length)
			error("subscript out of bounds");
		code = (unsigned char) CHAR(tag)[i];
		h = code + 1;
		if (h >= dectag_length
		    || (letter = CHAR(dectag)[h]) == dec_hole)
			error("unknown code %d in string to decode", code);
		CHAR(string)[j] = letter;
	}
	PROTECT(ans = allocVector(STRSXP, 1));
	SET_STRING_ELT(ans, 0, string);
	UNPROTECT(2);
	return ans;
}

/*
 * The bbuf_write_enc_chars() function is used when initializing
 * a BString object to encode and store the source string in the @data
 * slot of the object.
 * 'val' must be a non-empty single string (character vector of length 1).
 */
SEXP bbuf_write_enc_chars(SEXP bb_xp, SEXP imin, SEXP imax,
		SEXP val, SEXP enc_xp)
{
	SEXP tag, string, enctag;
	int string_length, enctag_length;
	int i1, i2, n, j, h;
	char enc_hole, letter, code;

	tag = R_ExternalPtrTag(bb_xp);
	i1 = INTEGER(imin)[0] - 1;
	i2 = INTEGER(imax)[0] - 1;
	if (i1 < 0 || i2 >= LENGTH(tag))
		error("subscript out of bounds");
	n = i2 - i1 + 1;
	string = STRING_ELT(val, 0);
	string_length = LENGTH(string);
	if (string_length == 0 && n != 0)
		error("no value provided");
	enctag = R_ExternalPtrTag(enc_xp);
	enctag_length = LENGTH(enctag);
	enc_hole = CHAR(enctag)[0];

	for (j = 0; i1 <= i2; i1++, j++) {
		if (j >= string_length)
			j = 0; /* recycle */
		letter = CHAR(string)[j];
		h = ((unsigned char) letter) + 1;
		if (h >= enctag_length
		    || (code = CHAR(enctag)[h]) == enc_hole)
			error("unknown letter '%c' in string to encode",
			      letter);
		CHAR(tag)[i1] = code;
	}
	if (j != string_length) {
		warning("number of items to replace is not a multiple "
			"of replacement length");
	}
	return bb_xp;
}

SEXP bbuf_writeii_enc_chars(SEXP bb_xp, SEXP ii, SEXP val, SEXP enc_xp)
{
	SEXP tag, string, enctag;
	int tag_length, string_length, enctag_length;
	int n, i, j, z, h;
	char enc_hole, letter, code;

	string = STRING_ELT(val, 0);
	string_length = LENGTH(string);
	n = LENGTH(ii);
	if (string_length == 0 && n != 0)
		error("no value provided");
	enctag = R_ExternalPtrTag(enc_xp);
	enctag_length = LENGTH(enctag);
	enc_hole = CHAR(enctag)[0];
	tag = R_ExternalPtrTag(bb_xp);
	tag_length = LENGTH(tag);

	for (j = z = 0; z < n; j++, z++) {
		i = INTEGER(ii)[z] - 1;
		if (i < 0 || i >= tag_length)
			error("subscript out of bounds");
		if (j >= string_length)
			j = 0; /* recycle */
		letter = CHAR(string)[j];
		h = ((unsigned char) letter) + 1;
		if (h >= enctag_length
		    || (code = CHAR(enctag)[h]) == enc_hole)
			error("unknown letter '%c' in string to encode",
			      letter);
		CHAR(tag)[i] = code;
	}
	if (j != string_length) {
		warning("number of items to replace is not a multiple "
			"of replacement length");
	}
	return bb_xp;
}
