#include "Biostrings.h"


static int debug = 0;

SEXP CharBuffer_debug()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'CharBuffer.c'\n", debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'CharBuffer.c'\n");
#endif
	return R_NilValue;
}

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
 *   cb <- CharBuffer(5)
 *   cb[] <- 255
 *   s <- toString(cb)
 *   strsplit(s, NULL, fixed=TRUE)[[1]]
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
 * Memory allocation for a "CharBuffer" object
 * -------------------------------------------
 * A "CharBuffer" object stores its data in an "external" R string (CHARSXP
 * vector). An R string itself stores its data in a char-array.
 * The "R types" of the argument passed to these functions must be:
 *   'cb_xp': externalptr
 *   'length': single integer
 */

/*
 * Alloc an R string of length 'length' and point 'cb_xp' to it.
 * Does NOT initialize the allocated memory!
 */
SEXP CharBuffer_alloc(SEXP cb_xp, SEXP length)
{
	SEXP tag;
	int tag_length;

	tag_length = INTEGER(length)[0];

	PROTECT(tag = allocString(tag_length));
	/*
	Rprintf("Memory successfully allocated for %d-byte CharBuffer object (starting at address %p)\n",
		tag_length, CHAR(tag));
	 */
	R_SetExternalPtrTag(cb_xp, tag);
	UNPROTECT(1);
	return cb_xp;
}

/*
 * Return the single string printed by the show method for "CharBuffer" objects.
 * 'cb_xp' must be the 'xp' slot of a "CharBuffer" object.
 * From R:
 *   cb <- CharBuffer(30)
 *   .Call("CharBuffer_get_show_string", cb@xp, PACKAGE="Biostrings")
 */
SEXP CharBuffer_get_show_string(SEXP cb_xp)
{
	SEXP tag, string, ans;
	int tag_length;

	tag = R_ExternalPtrTag(cb_xp);
	tag_length = LENGTH(tag);
	PROTECT(string = allocString(100));
	sprintf(CHAR(string), "%d-byte CharBuffer object (starting at address %p)",
		tag_length, CHAR(tag));
	PROTECT(ans = allocVector(STRSXP, 1));
	SET_STRING_ELT(ans, 0, string);
	UNPROTECT(2);
	return ans;
}

/*
 * Return length of R string pointed by 'cb_xp'.
 * From R:
 *   cb <- CharBuffer(30)
 *   .Call("CharBuffer_length", cb@xp, PACKAGE="Biostrings")
 * Called by method length() for "CharBuffer" objects.
 */
SEXP CharBuffer_length(SEXP cb_xp)
{
	SEXP tag, ans;
	int tag_length;

	tag = R_ExternalPtrTag(cb_xp);
	tag_length = LENGTH(tag);

	PROTECT(ans = allocVector(INTSXP, 1));
	INTEGER(ans)[0] = tag_length;
	UNPROTECT(1);
	return ans;
}

/*
 * From R:
 *   cb <- CharBuffer(30)
 *   .Call("CharBuffer_memcmp", cb@xp, 1:1, cb@xp, 10:10, 21:21, PACKAGE="Biostrings")
 */
SEXP CharBuffer_memcmp(SEXP cb1_xp, SEXP start1,
		 SEXP cb2_xp, SEXP start2, SEXP width)
{
	SEXP tag1, tag2, ans;
	int i1, i2, n;

#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] CharBuffer_memcmp(): BEGIN\n");
	}
#endif
	tag1 = R_ExternalPtrTag(cb1_xp);
	i1 = INTEGER(start1)[0] - 1;
	tag2 = R_ExternalPtrTag(cb2_xp);
	i2 = INTEGER(start2)[0] - 1;
	n = INTEGER(width)[0];

#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] CharBuffer_memcmp(): ");
		Rprintf("CHAR(tag1)=%p i1=%d CHAR(tag2)=%p i2=%d n=%d\n",
			CHAR(tag1), i1, CHAR(tag2), i2, n);
	}
#endif
	PROTECT(ans = allocVector(INTSXP, 1));
	INTEGER(ans)[0] = Biostrings_memcmp(CHAR(tag1), i1,
					CHAR(tag2), i2,
					n, sizeof(char));
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] CharBuffer_memcmp(): END\n");
	}
#endif
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * READ/WRITE functions
 * ====================
 * The functions in this section implement the read/write operations to a
 * "CharBuffer" object. The user can choose between 2 interfaces for each
 * read or write operation:
 *
 *   1. The "i1i2" interface: the chars to access are specified by 2
 * integers: 'imin' (the position of the first char to access, the first
 * char in the buffer being at position 1) and 'imax' (the position of the
 * last char to access).
 *
 *   2. The "subset" interface: the chars to access are specified by an
 * integer vector containing their positions in the buffer.
 *
 * The "subset" interface is intended to be used by the subsetting
 * operator [ defined at the R level for "CharBuffer" objects.
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
 * Functions in this section that implement the "subset" interface
 * respect the above properties.
 *
 * Here are some arguments to these functions that must always be SEXP of the
 * following types:
 *   src_xp, dest_xp: externalptr
 *   imin, imax: single integers
 *   subset: integer vector containing the subscripts (with no NAs)
 *   hash_xp: externalptr (hash table for encoding/decoding)
 */


/* ==========================================================================
 * Copy chars from a "CharBuffer" to another "CharBuffer" object.
 * --------------------------------------------------------------------------
 */

SEXP CharBuffer_copy_from_i1i2(SEXP dest_xp, SEXP src_xp, SEXP imin, SEXP imax)
{
	SEXP dest, src;
	int i1, i2;

	dest = R_ExternalPtrTag(dest_xp);
	src = R_ExternalPtrTag(src_xp);
	i1 = INTEGER(imin)[0] - 1;
	i2 = INTEGER(imax)[0] - 1;
	Biostrings_memcpy_from_i1i2(i1, i2,
			CHAR(dest), LENGTH(dest),
			CHAR(src), LENGTH(src), sizeof(char));
	return dest_xp;
}

SEXP CharBuffer_copy_from_subset(SEXP dest_xp, SEXP src_xp, SEXP subset)
{
	SEXP dest, src;

	dest = R_ExternalPtrTag(dest_xp);
	src = R_ExternalPtrTag(src_xp);
	Biostrings_memcpy_from_subset(INTEGER(subset), LENGTH(subset),
			CHAR(dest), LENGTH(dest),
			CHAR(src), LENGTH(src), sizeof(char));
	return dest_xp;
}


/* ==========================================================================
 * Read/write chars from/to a "CharBuffer" object
 * --------------------------------------------------------------------------
 */

/*
 * Return a single string (character vector of length 1).
 * From R:
 *   cb <- CharBuffer(15)
 *   cb[] < "Hello"
 *   .Call("CharBuffer_read_chars_from_i1i2", cb@xp, 2:2, 4:4, PACKAGE="Biostrings")
 */
SEXP CharBuffer_read_chars_from_i1i2(SEXP src_xp, SEXP imin, SEXP imax)
{
	SEXP src, dest, ans;
	int i1, i2, n;

	src = R_ExternalPtrTag(src_xp);
	i1 = INTEGER(imin)[0] - 1;
	i2 = INTEGER(imax)[0] - 1;
	n = i2 - i1 + 1;
	PROTECT(dest = allocString(n));
	Biostrings_memcpy_from_i1i2(i1, i2,
			CHAR(dest), n,
			CHAR(src), LENGTH(src), sizeof(char));
	PROTECT(ans = allocVector(STRSXP, 1));
	SET_STRING_ELT(ans, 0, dest);
	UNPROTECT(2);
	return ans;
}

SEXP CharBuffer_read_chars_from_subset(SEXP src_xp, SEXP subset)
{
	SEXP src, dest, ans;
	int n;

	src = R_ExternalPtrTag(src_xp);
	n = LENGTH(subset);
	PROTECT(dest = allocString(n));
	Biostrings_memcpy_from_subset(INTEGER(subset), n,
			CHAR(dest), n,
			CHAR(src), LENGTH(src), sizeof(char));
	PROTECT(ans = allocVector(STRSXP, 1));
	SET_STRING_ELT(ans, 0, dest);
	UNPROTECT(2);
	return ans;
}

/*
 * 'string' must be a non-empty single string (character vector of length 1).
 */
SEXP CharBuffer_write_chars_to_i1i2(SEXP dest_xp, SEXP imin, SEXP imax, SEXP string)
{
	SEXP dest, src;
	int i1, i2;

	dest = R_ExternalPtrTag(dest_xp);
	i1 = INTEGER(imin)[0] - 1;
	i2 = INTEGER(imax)[0] - 1;
	src = STRING_ELT(string, 0);
	Biostrings_memcpy_to_i1i2(i1, i2,
			CHAR(dest), LENGTH(dest),
			CHAR(src), LENGTH(src), sizeof(char));
	return dest_xp;
}

SEXP CharBuffer_write_chars_to_subset(SEXP dest_xp, SEXP subset, SEXP string)
{
	SEXP dest, src;

	dest = R_ExternalPtrTag(dest_xp);
	src = STRING_ELT(string, 0);
	Biostrings_memcpy_to_subset(INTEGER(subset), LENGTH(subset),
			CHAR(dest), LENGTH(dest),
			CHAR(src), LENGTH(src), sizeof(char));
	return dest_xp;
}


/* ==========================================================================
 * Read/write integers from/to a "CharBuffer" object
 * --------------------------------------------------------------------------
 */

/*
 * Return an integer vector of length 'imax' - 'imin' + 1.
 * From R:
 *   cb <- CharBuffer(30)
 *   .Call("CharBuffer_read_ints_from_i1i2", cb@xp, 20:20, 25:25, PACKAGE="Biostrings")
 */
SEXP CharBuffer_read_ints_from_i1i2(SEXP src_xp, SEXP imin, SEXP imax)
{
	SEXP src, ans;
	int i1, i2, n, j;

	src = R_ExternalPtrTag(src_xp);
	i1 = INTEGER(imin)[0] - 1;
	i2 = INTEGER(imax)[0] - 1;
	if (i1 < 0 || i2 >= LENGTH(src))
		error("subscript out of bounds");
	n = i2 - i1 + 1;

	PROTECT(ans = allocVector(INTSXP, n));
	for (j = 0; i1 <= i2; i1++, j++) {
		INTEGER(ans)[j] = (unsigned char) CHAR(src)[i1];
	}
	UNPROTECT(1);
	return ans;
}

/*
 * Return an integer vector of same length than 'subset'.
 * From R:
 *   cb <- CharBuffer(30)
 *   .Call("CharBuffer_read_ints_from_subset", cb, 25:20, PACKAGE="Biostrings")
 */
SEXP CharBuffer_read_ints_from_subset(SEXP src_xp, SEXP subset)
{
	SEXP src, ans;
	int src_length;
	int n, i, j;

	src = R_ExternalPtrTag(src_xp);
	src_length = LENGTH(src);
	n = LENGTH(subset);

	PROTECT(ans = allocVector(INTSXP, n));
	for (j = 0; j < n; j++) {
		i = INTEGER(subset)[j] - 1;
		if (i < 0 || i >= src_length)
			error("subscript out of bounds");
		INTEGER(ans)[j] = (unsigned char) CHAR(src)[i];
	}
	UNPROTECT(1);
	return ans;
}

/*
 * 'val' must be an integer vector of length > 0.
 */
SEXP CharBuffer_write_ints_to_i1i2(SEXP dest_xp, SEXP imin, SEXP imax, SEXP val)
{
	SEXP dest;
	int val_length;
	int i1, i2, n, j;
	int v;

	dest = R_ExternalPtrTag(dest_xp);
	i1 = INTEGER(imin)[0] - 1;
	i2 = INTEGER(imax)[0] - 1;
	if (i1 < 0 || i2 >= LENGTH(dest))
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
		CHAR(dest)[i1] = (char) v;
	}
	if (j != val_length) {
		warning("number of items to replace is not a multiple "
			"of replacement length");
	}
	return dest_xp;
}

SEXP CharBuffer_write_ints_to_subset(SEXP dest_xp, SEXP subset, SEXP val)
{
	SEXP dest;
	int dest_length, val_length;
	int n, i, j, z;
	int v;

	val_length = LENGTH(val);
	n = LENGTH(subset);
	if (val_length == 0 && n != 0)
		error("no value provided");
	dest = R_ExternalPtrTag(dest_xp);
	dest_length = LENGTH(dest);

	for (j = z = 0; z < n; j++, z++) {
		i = INTEGER(subset)[z] - 1;
		if (i < 0 || i >= dest_length)
			error("subscript out of bounds");
		if (j >= val_length)
			j = 0; /* recycle */
		v = INTEGER(val)[j];
		if (v < 0 || v >= 256)
			error("value out of range");
		CHAR(dest)[i] = (char) v;
	}
	if (j != val_length) {
		warning("number of items to replace is not a multiple "
			"of replacement length");
	}
	return dest_xp;
}


/* ==========================================================================
 * Read/write encoded chars from/to a "CharBuffer" object
 * --------------------------------------------------------------------------
 */

/*
 * Return a single string (character vector of length 1).
 */
SEXP CharBuffer_read_enc_chars_from_i1i2(SEXP src_xp, SEXP imin, SEXP imax, SEXP hash_xp)
{
	SEXP dest, src, hash, ans;
	int i1, i2, n;
	char hash_hole;

	src = R_ExternalPtrTag(src_xp);
	i1 = INTEGER(imin)[0] - 1;
	i2 = INTEGER(imax)[0] - 1;
	n = i2 - i1 + 1;
	hash = R_ExternalPtrTag(hash_xp);
	hash_hole = CHAR(hash)[0];
	PROTECT(dest = allocString(n));
	Biostrings_translate_charcpy_from_i1i2(i1, i2,
		CHAR(dest), n, CHAR(src), LENGTH(src),
		CHAR(hash) + 1, LENGTH(hash) - 1, hash_hole);
	PROTECT(ans = allocVector(STRSXP, 1));
	SET_STRING_ELT(ans, 0, dest);
	UNPROTECT(2);
	return ans;
}

SEXP CharBuffer_read_enc_chars_from_subset(SEXP src_xp, SEXP subset, SEXP hash_xp)
{
	SEXP dest, src, hash, ans;
	int n;
	char hash_hole;

	src = R_ExternalPtrTag(src_xp);
	n = LENGTH(subset);
	hash = R_ExternalPtrTag(hash_xp);
	hash_hole = CHAR(hash)[0];
	PROTECT(dest = allocString(n));
	Biostrings_translate_charcpy_from_subset(INTEGER(subset), n,
		CHAR(dest), LENGTH(dest), CHAR(src), LENGTH(src),
		CHAR(hash) + 1, LENGTH(hash) - 1, hash_hole);
	PROTECT(ans = allocVector(STRSXP, 1));
	SET_STRING_ELT(ans, 0, dest);
	UNPROTECT(2);
	return ans;
}

/*
 * The CharBuffer_write_enc_chars_to_i1i2() function is used when initializing
 * a BString object to encode and store the source string in the @data
 * slot of the object.
 * 'string' must be a non-empty single string (character vector of length 1).
 */
SEXP CharBuffer_write_enc_chars_to_i1i2(SEXP dest_xp, SEXP imin, SEXP imax,
		SEXP string, SEXP hash_xp)
{
	SEXP dest, src, hash;
	int i1, i2, n;
	char hash_hole;

	dest = R_ExternalPtrTag(dest_xp);
	i1 = INTEGER(imin)[0] - 1;
	i2 = INTEGER(imax)[0] - 1;
	n = i2 - i1 + 1;
	src = STRING_ELT(string, 0);
	hash = R_ExternalPtrTag(hash_xp);
	hash_hole = CHAR(hash)[0];
	Biostrings_translate_charcpy_to_i1i2(i1, i2,
		CHAR(dest), n, CHAR(src), LENGTH(src),
		CHAR(hash) + 1, LENGTH(hash) - 1, hash_hole);
	return dest_xp;
}

SEXP CharBuffer_write_enc_chars_to_subset(SEXP dest_xp, SEXP subset,
		SEXP string, SEXP hash_xp)
{
	SEXP dest, src, hash;
	int n;
	char hash_hole;

	dest = R_ExternalPtrTag(dest_xp);
	n = LENGTH(subset);
	src = STRING_ELT(string, 0);
	hash = R_ExternalPtrTag(hash_xp);
	hash_hole = CHAR(hash)[0];
	Biostrings_translate_charcpy_to_subset(INTEGER(subset), n,
		CHAR(dest), LENGTH(dest), CHAR(src), LENGTH(src),
		CHAR(hash) + 1, LENGTH(hash) - 1, hash_hole);
	return dest_xp;
}


/* ==========================================================================
 * Read chars from a "CharBuffer" object and convert them to a vector
 * of complexes.
 * --------------------------------------------------------------------------
 */

SEXP CharBuffer_read_complexes_from_i1i2(SEXP src_xp, SEXP imin, SEXP imax, SEXP lookup_table)
{
	SEXP dest, src;
	int i1, i2, n;

	src = R_ExternalPtrTag(src_xp);
	i1 = INTEGER(imin)[0] - 1;
	i2 = INTEGER(imax)[0] - 1;
	n = i2 - i1 + 1;
	PROTECT(dest = allocVector(CPLXSXP, n));
	Biostrings_coerce_to_complex_from_i1i2(i1, i2,
		COMPLEX(dest), n, CHAR(src), LENGTH(src),
		COMPLEX(lookup_table), LENGTH(lookup_table));
	UNPROTECT(1);
	return dest;
}

SEXP CharBuffer_read_complexes_from_subset(SEXP src_xp, SEXP subset, SEXP lookup_table)
{
	SEXP dest, src;
	int n;

	src = R_ExternalPtrTag(src_xp);
	n = LENGTH(subset);
	PROTECT(dest = allocVector(CPLXSXP, n));
	error("not available yet");
	UNPROTECT(1);
	return dest;
}

