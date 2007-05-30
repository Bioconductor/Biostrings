#include "Biostrings.h"


static int debug = 0;

SEXP XRaw_debug()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'XRaw.c'\n", debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'XRaw.c'\n");
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
	SEXP ans;
	char buf[40]; /* should be enough, even for 128-bit addresses */

	snprintf(buf, sizeof(buf), "%p", s);
	PROTECT(ans = allocVector(STRSXP, 1));
	SET_STRING_ELT(ans, 0, mkChar(buf));
	UNPROTECT(1);
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
 *   > r <- as.raw(c(10, 255))
 *   > s <- rawToChar(r)
 *   > s
 *   [1] "\n\xff"
 *   > strsplit(s, NULL, fixed=TRUE)[[1]]
 *   [1] NA
 * doesn't work!
 * The function below should be safe, whatever the content of 's' is!
 * The length of the returned string is the number of chars in single
 * string s. Not vectorized.
 */
SEXP safe_explode(SEXP s)
{
	SEXP s0, ans;
	int s0_length, i;
	char buf[2] = "X"; /* we only care about having buf[1] == 0 */

	s0 = STRING_ELT(s, 0);
	s0_length = LENGTH(s0);

	PROTECT(ans = allocVector(STRSXP, s0_length));
	for (i = 0; i < s0_length; i++) {
		buf[0] = CHAR(s0)[i];
		SET_STRING_ELT(ans, i, mkChar(buf));
	}
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * Memory allocation for a "XRaw" object
 * -------------------------------------
 * An "XRaw" object stores its data in an "external" raw vector (RAWSXP
 * vector). A RAWSXP vector itself stores its data in a char-array.
 * The "R types" of the argument passed to these functions must be:
 *   'xraw_xp': externalptr
 *   'length': single integer
 */

/*
 * Alloc a RAWSXP vector of length 'length' and point 'xraw_xp' to it.
 * Does NOT initialize the allocated memory!
 */
SEXP XRaw_alloc(SEXP xraw_xp, SEXP length)
{
	SEXP tag;
	int tag_length;

	tag_length = INTEGER(length)[0];

	PROTECT(tag = NEW_RAW(tag_length));
	/*
	Rprintf("Memory successfully allocated for %d-byte XRaw object (starting at address %p)\n",
		tag_length, RAW(tag));
	 */
	R_SetExternalPtrTag(xraw_xp, tag);
	UNPROTECT(1);
	return xraw_xp;
}

/*
 * Return the single string printed by the show method for "XRaw" objects.
 * 'xraw_xp' must be the 'xp' slot of a "XRaw" object.
 * From R:
 *   xr <- XRaw(30)
 *   .Call("XRaw_get_show_string", xr@xp, PACKAGE="Biostrings")
 */
SEXP XRaw_get_show_string(SEXP xraw_xp)
{
	SEXP tag, ans;
	int tag_length;
	char buf[100]; /* should be enough... */

	tag = R_ExternalPtrTag(xraw_xp);
	tag_length = LENGTH(tag);
	snprintf(buf, sizeof(buf), "%d-byte XRaw object (starting at address %p)",
		 tag_length, RAW(tag));
	PROTECT(ans = allocVector(STRSXP, 1));
	SET_STRING_ELT(ans, 0, mkChar(buf));
	UNPROTECT(1);
	return ans;
}

/*
 * Return length of R string pointed by 'xraw_xp'.
 * From R:
 *   xr <- XRaw(30)
 *   .Call("XRaw_length", xr@xp, PACKAGE="Biostrings")
 * Called by method length() for "XRaw" objects.
 */
SEXP XRaw_length(SEXP xraw_xp)
{
	SEXP tag, ans;
	int tag_length;

	tag = R_ExternalPtrTag(xraw_xp);
	tag_length = LENGTH(tag);

	PROTECT(ans = allocVector(INTSXP, 1));
	INTEGER(ans)[0] = tag_length;
	UNPROTECT(1);
	return ans;
}

/*
 * From R:
 *   xr <- XRaw(30)
 *   .Call("XRaw_memcmp", xr@xp, 1:1, xr@xp, 10:10, 21:21, PACKAGE="Biostrings")
 */
SEXP XRaw_memcmp(SEXP xraw1_xp, SEXP start1,
		 SEXP xraw2_xp, SEXP start2, SEXP width)
{
	SEXP tag1, tag2, ans;
	int i1, i2, n;

#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] XRaw_memcmp(): BEGIN\n");
	}
#endif
	tag1 = R_ExternalPtrTag(xraw1_xp);
	i1 = INTEGER(start1)[0] - 1;
	tag2 = R_ExternalPtrTag(xraw2_xp);
	i2 = INTEGER(start2)[0] - 1;
	n = INTEGER(width)[0];

#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] XRaw_memcmp(): ");
		Rprintf("RAW(tag1)=%p i1=%d RAW(tag2)=%p i2=%d n=%d\n",
			RAW(tag1), i1, RAW(tag2), i2, n);
	}
#endif
	PROTECT(ans = allocVector(INTSXP, 1));
	INTEGER(ans)[0] = Biostrings_memcmp((char *) RAW(tag1), i1,
					(char *) RAW(tag2), i2,
					n, sizeof(Rbyte));
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] XRaw_memcmp(): END\n");
	}
#endif
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * READ/WRITE functions
 * ====================
 * The functions in this section implement the read/write operations to a
 * "XRaw" object. The user can choose between 2 interfaces for each
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
 * operator [ defined at the R level for "XRaw" objects.
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
 *   src_xraw_xp, dest_xraw_xp: externalptr
 *   imin, imax: single integers
 *   subset: integer vector containing the subscripts (with no NAs)
 *   lkup: lookup table for encoding/decoding (integer or complex vector)
 */


/* ==========================================================================
 * Copy bytes from a "XRaw" to another "XRaw" object.
 * --------------------------------------------------------------------------
 */

SEXP XRaw_copy_from_i1i2(SEXP dest_xraw_xp, SEXP src_xraw_xp, SEXP imin, SEXP imax)
{
	SEXP dest, src;
	int i1, i2;

	dest = R_ExternalPtrTag(dest_xraw_xp);
	src = R_ExternalPtrTag(src_xraw_xp);
	i1 = INTEGER(imin)[0] - 1;
	i2 = INTEGER(imax)[0] - 1;
	Biostrings_memcpy_from_i1i2(i1, i2,
			(char *) RAW(dest), LENGTH(dest),
			(char *) RAW(src), LENGTH(src), sizeof(Rbyte));
	return dest_xraw_xp;
}

SEXP XRaw_copy_from_subset(SEXP dest_xraw_xp, SEXP src_xraw_xp, SEXP subset)
{
	SEXP dest, src;

	dest = R_ExternalPtrTag(dest_xraw_xp);
	src = R_ExternalPtrTag(src_xraw_xp);
	Biostrings_memcpy_from_subset(INTEGER(subset), LENGTH(subset),
			(char *) RAW(dest), LENGTH(dest),
			(char *) RAW(src), LENGTH(src), sizeof(Rbyte));
	return dest_xraw_xp;
}


/* ==========================================================================
 * Read/write chars from/to a "XRaw" object.
 * All the functions in this group assume that sizeof(Rbyte) == sizeof(char).
 * --------------------------------------------------------------------------
 */

/*
 * Return a single string (character vector of length 1).
 * From R:
 *   xr <- XRaw(15)
 *   xr[] < "Hello"
 *   .Call("XRaw_read_chars_from_i1i2", xr@xp, 2:2, 4:4, PACKAGE="Biostrings")
 */
SEXP XRaw_read_chars_from_i1i2(SEXP src_xraw_xp, SEXP imin, SEXP imax)
{
	SEXP src, ans;
	int i1, i2, n;
	char *dest;

	src = R_ExternalPtrTag(src_xraw_xp);
	i1 = INTEGER(imin)[0] - 1;
	i2 = INTEGER(imax)[0] - 1;
	n = i2 - i1 + 1;
	dest = R_ALLOC_STRING(n);
	/* assumes that sizeof(Rbyte) == sizeof(char) */
	Biostrings_memcpy_from_i1i2(i1, i2,
			dest, n, (char *) RAW(src), LENGTH(src),
			sizeof(char));
	PROTECT(ans = allocVector(STRSXP, 1));
	SET_STRING_ELT(ans, 0, mkChar(dest));
	UNPROTECT(1);
	return ans;
}

SEXP XRaw_read_chars_from_subset(SEXP src_xraw_xp, SEXP subset)
{
	SEXP src, ans;
	int n;
	char *dest;

	src = R_ExternalPtrTag(src_xraw_xp);
	n = LENGTH(subset);
	dest = R_ALLOC_STRING(n);
	/* assumes that sizeof(Rbyte) == sizeof(char) */
	Biostrings_memcpy_from_subset(INTEGER(subset), n,
			dest, n, (char *) RAW(src), LENGTH(src),
			sizeof(char));
	PROTECT(ans = allocVector(STRSXP, 1));
	SET_STRING_ELT(ans, 0, mkChar(dest));
	UNPROTECT(1);
	return ans;
}

/*
 * 'string' must be a non-empty single string (character vector of length 1).
 */
SEXP XRaw_write_chars_to_i1i2(SEXP dest_xraw_xp, SEXP imin, SEXP imax, SEXP string)
{
	SEXP dest, src;
	int i1, i2;

	dest = R_ExternalPtrTag(dest_xraw_xp);
	i1 = INTEGER(imin)[0] - 1;
	i2 = INTEGER(imax)[0] - 1;
	src = STRING_ELT(string, 0);
	/* assumes that sizeof(Rbyte) == sizeof(char) */
	Biostrings_memcpy_to_i1i2(i1, i2,
			(char *) RAW(dest), LENGTH(dest),
			CHAR(src), LENGTH(src), sizeof(char));
	return dest_xraw_xp;
}

SEXP XRaw_write_chars_to_subset(SEXP dest_xraw_xp, SEXP subset, SEXP string)
{
	SEXP dest, src;

	dest = R_ExternalPtrTag(dest_xraw_xp);
	src = STRING_ELT(string, 0);
	/* assumes that sizeof(Rbyte) == sizeof(char) */
	Biostrings_memcpy_to_subset(INTEGER(subset), LENGTH(subset),
			(char *) RAW(dest), LENGTH(dest),
			CHAR(src), LENGTH(src), sizeof(char));
	return dest_xraw_xp;
}


/* ==========================================================================
 * Read/write integers from/to a "XRaw" object
 * --------------------------------------------------------------------------
 */

/*
 * Return an integer vector of length 'imax' - 'imin' + 1.
 * From R:
 *   xr <- XRaw(30)
 *   .Call("XRaw_read_ints_from_i1i2", xr@xp, 20:20, 25:25, PACKAGE="Biostrings")
 */
SEXP XRaw_read_ints_from_i1i2(SEXP src_xraw_xp, SEXP imin, SEXP imax)
{
	SEXP src, ans;
	int i1, i2, n, j;

	src = R_ExternalPtrTag(src_xraw_xp);
	i1 = INTEGER(imin)[0] - 1;
	i2 = INTEGER(imax)[0] - 1;
	if (i1 < 0 || i2 >= LENGTH(src))
		error("subscript out of bounds");
	n = i2 - i1 + 1;

	PROTECT(ans = allocVector(INTSXP, n));
	for (j = 0; i1 <= i2; i1++, j++) {
		INTEGER(ans)[j] = (unsigned char) RAW(src)[i1];
	}
	UNPROTECT(1);
	return ans;
}

/*
 * Return an integer vector of same length than 'subset'.
 * From R:
 *   xr <- XRaw(30)
 *   .Call("XRaw_read_ints_from_subset", xr, 25:20, PACKAGE="Biostrings")
 */
SEXP XRaw_read_ints_from_subset(SEXP src_xraw_xp, SEXP subset)
{
	SEXP src, ans;
	int src_length;
	int n, i, j;

	src = R_ExternalPtrTag(src_xraw_xp);
	src_length = LENGTH(src);
	n = LENGTH(subset);

	PROTECT(ans = allocVector(INTSXP, n));
	for (j = 0; j < n; j++) {
		i = INTEGER(subset)[j] - 1;
		if (i < 0 || i >= src_length)
			error("subscript out of bounds");
		INTEGER(ans)[j] = (unsigned char) RAW(src)[i];
	}
	UNPROTECT(1);
	return ans;
}

/*
 * 'val' must be an integer vector of length > 0.
 */
SEXP XRaw_write_ints_to_i1i2(SEXP dest_xraw_xp, SEXP imin, SEXP imax, SEXP val)
{
	SEXP dest;
	int val_length;
	int i1, i2, n, j;
	int v;

	dest = R_ExternalPtrTag(dest_xraw_xp);
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
		RAW(dest)[i1] = (char) v;
	}
	if (j != val_length) {
		warning("number of items to replace is not a multiple "
			"of replacement length");
	}
	return dest_xraw_xp;
}

SEXP XRaw_write_ints_to_subset(SEXP dest_xraw_xp, SEXP subset, SEXP val)
{
	SEXP dest;
	int dest_length, val_length;
	int n, i, j, z;
	int v;

	val_length = LENGTH(val);
	n = LENGTH(subset);
	if (val_length == 0 && n != 0)
		error("no value provided");
	dest = R_ExternalPtrTag(dest_xraw_xp);
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
		RAW(dest)[i] = (char) v;
	}
	if (j != val_length) {
		warning("number of items to replace is not a multiple "
			"of replacement length");
	}
	return dest_xraw_xp;
}


/* ==========================================================================
 * Read/write encoded chars from/to a "XRaw" object
 * --------------------------------------------------------------------------
 */

/*
 * Return a single string (character vector of length 1).
 */
SEXP XRaw_read_enc_chars_from_i1i2(SEXP src_xraw_xp, SEXP imin, SEXP imax, SEXP lkup)
{
	SEXP src, ans;
	int i1, i2, n;
	char *dest;

	src = R_ExternalPtrTag(src_xraw_xp);
	i1 = INTEGER(imin)[0] - 1;
	i2 = INTEGER(imax)[0] - 1;
	n = i2 - i1 + 1;
	dest = R_ALLOC_STRING(n);
	Biostrings_translate_charcpy_from_i1i2(i1, i2,
			dest, n, (char *) RAW(src), LENGTH(src),
			INTEGER(lkup), LENGTH(lkup));
	PROTECT(ans = allocVector(STRSXP, 1));
	SET_STRING_ELT(ans, 0, mkChar(dest));
	UNPROTECT(1);
	return ans;
}

SEXP XRaw_read_enc_chars_from_subset(SEXP src_xraw_xp, SEXP subset, SEXP lkup)
{
	SEXP src, ans;
	int n;
	char *dest;

	src = R_ExternalPtrTag(src_xraw_xp);
	n = LENGTH(subset);
	dest = R_ALLOC_STRING(n);
	Biostrings_translate_charcpy_from_subset(INTEGER(subset), n,
			dest, n, (char *) RAW(src), LENGTH(src),
			INTEGER(lkup), LENGTH(lkup));
	PROTECT(ans = allocVector(STRSXP, 1));
	SET_STRING_ELT(ans, 0, mkChar(dest));
	UNPROTECT(1);
	return ans;
}

/*
 * The XRaw_write_enc_chars_to_i1i2() function is used when initializing
 * a BString object to encode and store the source string in the @data
 * slot of the object.
 * 'string' must be a non-empty single string (character vector of length 1).
 */
SEXP XRaw_write_enc_chars_to_i1i2(SEXP dest_xraw_xp, SEXP imin, SEXP imax,
		SEXP string, SEXP lkup)
{
	SEXP dest, src;
	int i1, i2, n;

	dest = R_ExternalPtrTag(dest_xraw_xp);
	i1 = INTEGER(imin)[0] - 1;
	i2 = INTEGER(imax)[0] - 1;
	n = i2 - i1 + 1;
	src = STRING_ELT(string, 0);
	Biostrings_translate_charcpy_to_i1i2(i1, i2,
		(char *) RAW(dest), n, CHAR(src), LENGTH(src),
		INTEGER(lkup), LENGTH(lkup));
	return dest_xraw_xp;
}

SEXP XRaw_write_enc_chars_to_subset(SEXP dest_xraw_xp, SEXP subset,
		SEXP string, SEXP lkup)
{
	SEXP dest, src;
	int n;

	dest = R_ExternalPtrTag(dest_xraw_xp);
	n = LENGTH(subset);
	src = STRING_ELT(string, 0);
	Biostrings_translate_charcpy_to_subset(INTEGER(subset), n,
		(char *) RAW(dest), LENGTH(dest), CHAR(src), LENGTH(src),
		INTEGER(lkup), LENGTH(lkup));
	return dest_xraw_xp;
}


/* ==========================================================================
 * Read chars from a "XRaw" object and convert them to a vector
 * of complexes.
 * --------------------------------------------------------------------------
 */

SEXP XRaw_read_complexes_from_i1i2(SEXP src_xraw_xp, SEXP imin, SEXP imax, SEXP lkup)
{
	SEXP dest, src;
	int i1, i2, n;

	src = R_ExternalPtrTag(src_xraw_xp);
	i1 = INTEGER(imin)[0] - 1;
	i2 = INTEGER(imax)[0] - 1;
	n = i2 - i1 + 1;
	PROTECT(dest = allocVector(CPLXSXP, n));
	Biostrings_coerce_to_complex_from_i1i2(i1, i2,
		COMPLEX(dest), n, (char *) RAW(src), LENGTH(src),
		COMPLEX(lkup), LENGTH(lkup));
	UNPROTECT(1);
	return dest;
}

SEXP XRaw_read_complexes_from_subset(SEXP src_xraw_xp, SEXP subset, SEXP lkup)
{
	SEXP dest, src;
	int n;

	src = R_ExternalPtrTag(src_xraw_xp);
	n = LENGTH(subset);
	PROTECT(dest = allocVector(CPLXSXP, n));
	error("not available yet");
	UNPROTECT(1);
	return dest;
}

