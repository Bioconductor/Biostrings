#include "Biostrings.h"

static int debug = 0;

SEXP Biostrings_debug_XRaw()
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
 *   .Call("Biostrings_sexp_address", 6:4, PACKAGE="Biostrings")
 *   .Call("Biostrings_sexp_address", new("externalptr"), PACKAGE="Biostrings")
 */
SEXP Biostrings_sexp_address(SEXP s)
{
	SEXP ans;
	char buf[40]; /* should be enough, even for 128-bit addresses */

	snprintf(buf, sizeof(buf), "%p", s);
	PROTECT(ans = NEW_CHARACTER(1));
	SET_STRING_ELT(ans, 0, mkChar(buf));
	UNPROTECT(1);
	return ans;
}

/*
 * Print some obscure info about an "externalptr" object.
 * From R:
 *   .Call("Biostrings_xp_show", new("externalptr"), PACKAGE="Biostrings")
 */
SEXP Biostrings_xp_show(SEXP xp)
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
 *   xp <- .Call("Biostrings_xp_new", PACKAGE="Biostrings")
 */
SEXP Biostrings_xp_new()
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
SEXP Biostrings_safe_explode(SEXP s)
{
	SEXP s0, ans;
	int s0_length, i;
	char buf[2] = "X"; /* we only care about having buf[1] == 0 */

	s0 = STRING_ELT(s, 0);
	s0_length = LENGTH(s0);

	PROTECT(ans = NEW_CHARACTER(s0_length));
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
SEXP Biostrings_XRaw_alloc(SEXP xraw_xp, SEXP length)
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
 *   .Call("Biostrings_XRaw_get_show_string", xr@xp, PACKAGE="Biostrings")
 */
SEXP Biostrings_XRaw_get_show_string(SEXP xraw_xp)
{
	SEXP tag, ans;
	int tag_length;
	char buf[100]; /* should be enough... */

	tag = R_ExternalPtrTag(xraw_xp);
	tag_length = LENGTH(tag);
	snprintf(buf, sizeof(buf), "%d-byte XRaw object (starting at address %p)",
		 tag_length, RAW(tag));
	PROTECT(ans = NEW_CHARACTER(1));
	SET_STRING_ELT(ans, 0, mkChar(buf));
	UNPROTECT(1);
	return ans;
}

/*
 * Return length of R string pointed by 'xraw_xp'.
 * From R:
 *   xr <- XRaw(30)
 *   .Call("Biostrings_XRaw_length", xr@xp, PACKAGE="Biostrings")
 * Called by method length() for "XRaw" objects.
 */
SEXP Biostrings_XRaw_length(SEXP xraw_xp)
{
	SEXP tag, ans;
	int tag_length;

	tag = R_ExternalPtrTag(xraw_xp);
	tag_length = LENGTH(tag);

	PROTECT(ans = NEW_INTEGER(1));
	INTEGER(ans)[0] = tag_length;
	UNPROTECT(1);
	return ans;
}

/*
 * From R:
 *   xr <- XRaw(30)
 *   .Call("Biostrings_XRaw_memcmp", xr@xp, 1:1, xr@xp, 10:10, 21:21, PACKAGE="Biostrings")
 */
SEXP Biostrings_XRaw_memcmp(SEXP xraw1_xp, SEXP start1,
		 SEXP xraw2_xp, SEXP start2, SEXP width)
{
	SEXP tag1, tag2, ans;
	int i1, i2, n;

#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] Biostrings_XRaw_memcmp(): BEGIN\n");
	}
#endif
	tag1 = R_ExternalPtrTag(xraw1_xp);
	i1 = INTEGER(start1)[0] - 1;
	tag2 = R_ExternalPtrTag(xraw2_xp);
	i2 = INTEGER(start2)[0] - 1;
	n = INTEGER(width)[0];

#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] Biostrings_XRaw_memcmp(): ");
		Rprintf("RAW(tag1)=%p i1=%d RAW(tag2)=%p i2=%d n=%d\n",
			RAW(tag1), i1, RAW(tag2), i2, n);
	}
#endif
	PROTECT(ans = NEW_INTEGER(1));
	INTEGER(ans)[0] = _Biostrings_memcmp((char *) RAW(tag1), i1,
					(char *) RAW(tag2), i2,
					n, sizeof(Rbyte));
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] Biostrings_XRaw_memcmp(): END\n");
	}
#endif
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * READ/WRITE functions
 * ====================
 *
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

SEXP Biostrings_XRaw_copy_from_i1i2(SEXP dest_xraw_xp, SEXP src_xraw_xp, SEXP imin, SEXP imax)
{
	SEXP dest, src;
	int i1, i2;

	dest = R_ExternalPtrTag(dest_xraw_xp);
	src = R_ExternalPtrTag(src_xraw_xp);
	i1 = INTEGER(imin)[0] - 1;
	i2 = INTEGER(imax)[0] - 1;
	_Biostrings_memcpy_from_i1i2(i1, i2,
			(char *) RAW(dest), LENGTH(dest),
			(char *) RAW(src), LENGTH(src), sizeof(Rbyte));
	return dest_xraw_xp;
}

SEXP Biostrings_XRaw_copy_from_subset(SEXP dest_xraw_xp, SEXP src_xraw_xp, SEXP subset)
{
	SEXP dest, src;

	dest = R_ExternalPtrTag(dest_xraw_xp);
	src = R_ExternalPtrTag(src_xraw_xp);
	_Biostrings_memcpy_from_subset(INTEGER(subset), LENGTH(subset),
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
 *   .Call("Biostrings_XRaw_read_chars_from_i1i2", xr@xp, 2:2, 4:4, PACKAGE="Biostrings")
 */
SEXP Biostrings_XRaw_read_chars_from_i1i2(SEXP src_xraw_xp, SEXP imin, SEXP imax)
{
	SEXP src, ans;
	int i1, i2, n;
	char *dest;

	src = R_ExternalPtrTag(src_xraw_xp);
	i1 = INTEGER(imin)[0] - 1;
	i2 = INTEGER(imax)[0] - 1;
	n = i2 - i1 + 1;
	dest = _Biostrings_alloc_string(n);
	/* assumes that sizeof(Rbyte) == sizeof(char) */
	_Biostrings_memcpy_from_i1i2(i1, i2,
			dest, n, (char *) RAW(src), LENGTH(src),
			sizeof(char));
	PROTECT(ans = NEW_CHARACTER(1));
	SET_STRING_ELT(ans, 0, mkChar(dest));
	UNPROTECT(1);
	return ans;
}

SEXP Biostrings_XRaw_read_chars_from_subset(SEXP src_xraw_xp, SEXP subset)
{
	SEXP src, ans;
	int n;
	char *dest;

	src = R_ExternalPtrTag(src_xraw_xp);
	n = LENGTH(subset);
	dest = _Biostrings_alloc_string(n);
	/* assumes that sizeof(Rbyte) == sizeof(char) */
	_Biostrings_memcpy_from_subset(INTEGER(subset), n,
			dest, n, (char *) RAW(src), LENGTH(src),
			sizeof(char));
	PROTECT(ans = NEW_CHARACTER(1));
	SET_STRING_ELT(ans, 0, mkChar(dest));
	UNPROTECT(1);
	return ans;
}

/*
 * 'string' must be a non-empty single string (character vector of length 1).
 */
SEXP Biostrings_XRaw_write_chars_to_i1i2(SEXP dest_xraw_xp, SEXP imin, SEXP imax, SEXP string)
{
	SEXP dest, src;
	int i1, i2;

	dest = R_ExternalPtrTag(dest_xraw_xp);
	i1 = INTEGER(imin)[0] - 1;
	i2 = INTEGER(imax)[0] - 1;
	src = STRING_ELT(string, 0);
	/* assumes that sizeof(Rbyte) == sizeof(char) */
	_Biostrings_memcpy_to_i1i2(i1, i2,
			(char *) RAW(dest), LENGTH(dest),
			CHAR(src), LENGTH(src), sizeof(char));
	return dest_xraw_xp;
}

SEXP Biostrings_XRaw_write_chars_to_subset(SEXP dest_xraw_xp, SEXP subset, SEXP string)
{
	SEXP dest, src;

	dest = R_ExternalPtrTag(dest_xraw_xp);
	src = STRING_ELT(string, 0);
	/* assumes that sizeof(Rbyte) == sizeof(char) */
	_Biostrings_memcpy_to_subset(INTEGER(subset), LENGTH(subset),
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

	PROTECT(ans = NEW_INTEGER(n));
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

	PROTECT(ans = NEW_INTEGER(n));
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
	dest = _Biostrings_alloc_string(n);
	_Biostrings_translate_charcpy_from_i1i2(i1, i2,
			dest, n, (char *) RAW(src), LENGTH(src),
			INTEGER(lkup), LENGTH(lkup));
	PROTECT(ans = NEW_CHARACTER(1));
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
	dest = _Biostrings_alloc_string(n);
	_Biostrings_translate_charcpy_from_subset(INTEGER(subset), n,
			dest, n, (char *) RAW(src), LENGTH(src),
			INTEGER(lkup), LENGTH(lkup));
	PROTECT(ans = NEW_CHARACTER(1));
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
	int i1, i2;

	dest = R_ExternalPtrTag(dest_xraw_xp);
	i1 = INTEGER(imin)[0] - 1;
	i2 = INTEGER(imax)[0] - 1;
	src = STRING_ELT(string, 0);
	_Biostrings_translate_charcpy_to_i1i2(i1, i2,
			(char *) RAW(dest), LENGTH(dest),
			CHAR(src), LENGTH(src),
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
	_Biostrings_translate_charcpy_to_subset(INTEGER(subset), n,
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
	PROTECT(dest = NEW_COMPLEX(n));
	_Biostrings_coerce_to_complex_from_i1i2(i1, i2,
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
	PROTECT(dest = NEW_COMPLEX(n));
	error("not available yet");
	UNPROTECT(1);
	return dest;
}


/****************************************************************************
 * 6 convenience functions
 * =======================
 *
 */

/* NOT a Call() entry point! */
SEXP mkXRaw(SEXP tag)
{
	SEXP ans;

	PROTECT(ans = NEW_OBJECT(MAKE_CLASS("XRaw")));
	SET_SLOT(ans, mkChar("xp"), R_MakeExternalPtr(NULL, tag, R_NilValue));
	UNPROTECT(1);
        return ans;
}

/* NOT a Call() entry point! */
SEXP CHARSXP_to_XRaw(SEXP x, SEXP start, SEXP nchar, SEXP lkup)
{
	SEXP ans, dest;
	int offset, length;
	const char *src;

	offset = _start2offset(INTEGER(start)[0]);
	length = _nchar2length(INTEGER(nchar)[0], offset, LENGTH(x));
	src = CHAR(x) + offset;

	PROTECT(dest = NEW_RAW(length));
	if (lkup == R_NilValue) {
		_Biostrings_memcpy_to_i1i2(0, length - 1,
			(char *) RAW(dest), length,
			src, length, sizeof(char));
	} else {
		_Biostrings_translate_charcpy_to_i1i2(0, length - 1,
			(char *) RAW(dest), length,
			src, length,
			INTEGER(lkup), LENGTH(lkup));
	}
	ans = mkXRaw(dest);
	UNPROTECT(1);
	return ans;
}

SEXP char_to_XRaw(SEXP x, SEXP start, SEXP nchar, SEXP lkup)
{
	return CHARSXP_to_XRaw(STRING_ELT(x, 0), start, nchar, lkup);
}

SEXP copy_subXRaw(SEXP x, SEXP start, SEXP nchar, SEXP lkup)
{
	SEXP ans;

	error("copy_subXRaw() not ready yet");
	return R_NilValue;
}

/*
 * No checking of 'start' and 'nchar' lengths! We assume they are always good...
 */
SEXP STRSXP_to_XRaw(SEXP x, SEXP start, SEXP nchar, SEXP lkup)
{
	SEXP dest, x_elt, ans;
	int x_len, dest_len, i, *start_elt, *nchar_elt, offset, length;
	char *dest_elt;
	const char *src;

	x_len = LENGTH(x);
	dest_len = 0;
	for (i = 0, nchar_elt = INTEGER(nchar); i < x_len; i++, nchar_elt++)
		dest_len += *nchar_elt;

	PROTECT(dest = NEW_RAW(dest_len));
	dest_elt = (char *) RAW(dest);
	for (i = 0, start_elt = INTEGER(start), nchar_elt = INTEGER(nchar);
             i < x_len;
             i++, start_elt++, nchar_elt++) {
		x_elt = STRING_ELT(x, i);
		offset = _start2offset(*start_elt);
		length = _nchar2length(*nchar_elt, offset, LENGTH(x_elt));
		src = CHAR(x_elt) + offset;
		if (lkup == R_NilValue) {
			_Biostrings_memcpy_to_i1i2(0, length - 1,
				dest_elt, length,
				src, length, sizeof(char));
		} else {
			_Biostrings_translate_charcpy_to_i1i2(0, length - 1,
				dest_elt, length,
				src, length,
				INTEGER(lkup), LENGTH(lkup));
		}
		dest_elt += length;
	}
	ans = mkXRaw(dest);
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * I/O functions
 * =============
 *
 */

#define FASTALINE_MAX 20000

/*
 XRaw_loadFASTA().
 Load a FASTA file into an XRaw object.

 Return a named list of 4 elements:

   - "nbyte" (single integer): number of bytes that were written to the XRaw
     object. XRaw_loadFASTA() starts to write at position 1 (the first byte)
     in the XRaw object. If nbyte_max is the length of the XRaw object, we
     should always have 1 <= nbyte <= nbyte_max. An error is raised if the
     FASTA file contains no data or if the size of the data exceeds the
     capacity (i.e. the length) of the XRaw object i.e. if XRaw_loadFASTA()
     tries to write more than nbyte_max bytes to the XRaw object (it does not
     try to resize it). 

   - "start" and "end": 2 integer vectors of same length L (L should always be
     >= 1) containing the start/end positions of the FASTA sequences relatively
     to the first byte of the XRaw object. The start/end positions define a set
     of views such that:
       (a) there is at least one view,
       (b) all views have a width >= 1 ('all(start <= end)' is TRUE),
       (c) the views are not overlapping and are sorted from left to right,
       (d) the leftmost view starts at position 1 (start[1] == 1) and the
           rightmost view ends at position nbyte (end[L] == nbyte).

   - "desc" (character vector): descriptions of the FASTA sequences. The length
     of 'desc' is L too.

 XRaw_loadFASTA() is designed to be very fast:

     library(Biostrings)
     filepath <- "~hpages/BSgenome_srcdata/UCSC/hg18/chr1.fa"
     filesize <- file.info(filepath)$size
     if (is.na(filesize))
         stop(filepath, ": file not found")
     filesize <- as.integer(filesize)
     if (is.na(filesize))
         stop(filepath, ": file is too big")
     xraw <- Biostrings:::XRaw(filesize)

     # Takes < 1 second on lamb1!
     .Call("XRaw_loadFASTA", xraw@xp, filepath, "", NULL, PACKAGE="Biostrings")

     # or use safe wrapper:
     Biostrings:::XRaw.loadFASTA(xraw, filepath)

   Compare to the time it takes to load Hsapiens$chr1 from the
   BSgenome.Hsapiens.UCSC.hg18 package: 1.340 second on lamb1!
   Conclusion: we could put and load directly the FASTA files in the
   BSgenome.* packages. The DNAString and BStringViews instances would be
   created on the fly. No need to store them in .rda files anymore!

*/

SEXP XRaw_loadFASTA(SEXP xraw_xp, SEXP filepath, SEXP collapse, SEXP lkup)
{
	SEXP ans, ans_names, ans_elt, dest;
	const char *path, *coll;
	FILE *infile;
	long int lineno;
	char line[FASTALINE_MAX+1], desc[FASTALINE_MAX+1];
	int nbyte_max, gaplen, line_len, status, view_start, i1, i2;
	char c0;

	dest = R_ExternalPtrTag(xraw_xp);
	nbyte_max = LENGTH(dest);
	path = translateChar(STRING_ELT(filepath, 0));
	coll = CHAR(STRING_ELT(collapse, 0));
	gaplen = strlen(coll);

	if ((infile = fopen(path, "r")) == NULL)
		error("cannot open file");
	lineno = i1 = 0;
	status = 0; /* 0: expecting desc; 1: expecting seq; 2: no expectation */
	_Biostrings_reset_viewsbuf(0);
	while ((line_len = fgets_rtrimmed(line, FASTALINE_MAX+1, infile)) != -1) {
	/* while (fgets(line, FASTALINE_MAX+1, infile) != NULL) { */
		lineno++;
		if (line_len >= FASTALINE_MAX) {
			fclose(infile);
			error("file contains lines that are too long");
		}
		if (line_len == 0)
			continue;
		c0 = line[0];
		if (c0 == ';')
			continue;
		if (c0 != '>') {
			if (status == 0) {
				fclose(infile);
				error("file does not seem to be FASTA");
			}
			i2 = i1 + line_len - 1;
			if (lkup == R_NilValue)
				_Biostrings_memcpy_to_i1i2(i1, i2,
					(char *) RAW(dest), nbyte_max,
					line, line_len, sizeof(char));
			else
				_Biostrings_translate_charcpy_to_i1i2(i1, i2,
					(char *) RAW(dest), nbyte_max,
					line, line_len,
					INTEGER(lkup), LENGTH(lkup));
			i1 = i2 + 1;
			status = 2;
			continue;
		}
		if (status == 1) {
			fclose(infile);
			error("file does not seem to be FASTA");
		}
		if (status == 2) {
			_Biostrings_append_view(view_start, i1, desc);
			if (gaplen != 0) {
				i2 = i1 + gaplen - 1;
				_Biostrings_memcpy_to_i1i2(i1, i2,
					(char *) RAW(dest), nbyte_max,
					coll, gaplen, sizeof(char));
				i1 = i2 + 1;
			}
		}
		view_start = i1 + 1;
		strcpy(desc, line + 1);
		status = 1;
	}
	fclose(infile);
	if (status != 2)
		error("file does not seem to be FASTA");
	_Biostrings_append_view(view_start, i1, desc);

	PROTECT(ans = NEW_LIST(4));
	/* set the names */
	PROTECT(ans_names = NEW_CHARACTER(4));
	SET_STRING_ELT(ans_names, 0, mkChar("nbyte"));
	SET_STRING_ELT(ans_names, 1, mkChar("start"));
	SET_STRING_ELT(ans_names, 2, mkChar("end"));
	SET_STRING_ELT(ans_names, 3, mkChar("desc"));
	SET_NAMES(ans, ans_names);
	UNPROTECT(1);
	/* set the "nbyte" element */
	PROTECT(ans_elt = NEW_INTEGER(1));
	INTEGER(ans_elt)[0] = i1;
	SET_ELEMENT(ans, 0, ans_elt);
	UNPROTECT(1);
	/* set the "start" element */
	PROTECT(ans_elt = _Biostrings_viewsbuf_start_asINTEGER());
	SET_ELEMENT(ans, 1, ans_elt);
	UNPROTECT(1);
	/* set the "end" element */
	PROTECT(ans_elt = _Biostrings_viewsbuf_end_asINTEGER());
	SET_ELEMENT(ans, 2, ans_elt);
	UNPROTECT(1);
	/* set the "desc" element */
	PROTECT(ans_elt = _Biostrings_viewsbuf_desc_asCHARACTER());
	SET_ELEMENT(ans, 3, ans_elt);
	UNPROTECT(1);
	/* ans is ready */
	UNPROTECT(1);
	return ans;
}

