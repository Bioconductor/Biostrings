/****************************************************************************
 *                    Basic manipulation of XRaw objects                    *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"
#include <stdlib.h> /* for realloc() */

static int debug = 0;

SEXP debug_XRaw_class()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'XRaw_class.c'\n", debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'XRaw_class.c'\n");
#endif
	return R_NilValue;
}


/****************************************************************************
 * Things in this section are not necessarily directly related to XRaw
 * objects but they needed to go somewhere so here they are...
 */

const char *_get_class(SEXP x)
{
	return CHAR(STRING_ELT(GET_CLASS(x), 0));
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


/****************************************************************************
 * Allocating memory for an XRaw object.
 *
 * An XRaw object stores its data in an "external" raw vector (RAWSXP
 * vector). A RAWSXP vector itself stores its data in a char-array.
 * The "R types" of the argument passed to these functions must be:
 *   'xraw_xp': externalptr
 *   'length': single integer
 */

/*
 * Alloc an RAWSXP vector of length 'length' and point 'xraw_xp' to it.
 * Does NOT initialize the allocated memory!
 */
SEXP Biostrings_XRaw_alloc(SEXP xraw_xp, SEXP length)
{
	SEXP tag;
	int tag_length;

	tag_length = INTEGER(length)[0];

	PROTECT(tag = NEW_RAW(tag_length));
	/*
	Rprintf("Memory successfully allocated for %d-byte XRaw object (data starting at memory address %p)\n",
		tag_length, RAW(tag));
	*/
	R_SetExternalPtrTag(xraw_xp, tag);
	UNPROTECT(1);
	return xraw_xp;
}


/****************************************************************************
 * Getting information about an XRaw object.
 */

SEXP _get_XRaw_tag(SEXP x)
{
	return R_ExternalPtrTag(GET_SLOT(x, install("xp")));
}

int _get_XRaw_length(SEXP x)
{
	return LENGTH(_get_XRaw_tag(x));
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
	snprintf(buf, sizeof(buf), "%d-byte XRaw object (data starting at memory address %p)",
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


/****************************************************************************
 * Making new XRaw objects.
 */

/*
 * Do NOT make this a .Call() entry point!
 * Its argument is NOT duplicated so it would be a disaster if it was
 * coming from the user space.
 */
SEXP _new_XRaw(SEXP tag)
{
	SEXP ans;

	PROTECT(ans = NEW_OBJECT(MAKE_CLASS("XRaw")));
	SET_SLOT(ans, mkChar("xp"), R_MakeExternalPtr(NULL, tag, R_NilValue));
	UNPROTECT(1);
        return ans;
}

SEXP _new_XRaw_from_RoSeqs(RoSeqs seqs, SEXP lkup)
{
	SEXP tag, ans;
	int tag_length, i;
	const RoSeq *seq;
	char *dest;

	tag_length = 0;
	for (i = 0, seq = seqs.elts; i < seqs.nelt; i++, seq++)
		tag_length += seq->nelt;
	PROTECT(tag = NEW_RAW(tag_length));
	dest = (char *) RAW(tag);
	for (i = 0, seq = seqs.elts; i < seqs.nelt; i++, seq++) {
		if (lkup == R_NilValue) {
			_Biostrings_memcpy_to_i1i2(0, seq->nelt - 1,
				dest, seq->nelt,
				seq->elts, seq->nelt, sizeof(char));
		} else {
			_Biostrings_translate_charcpy_to_i1i2(0, seq->nelt - 1,
				dest, seq->nelt,
				seq->elts, seq->nelt,
				INTEGER(lkup), LENGTH(lkup));
		}
		dest += seq->nelt;
	}
	PROTECT(ans = _new_XRaw(tag));
	UNPROTECT(2);
	return ans;
}


/****************************************************************************
 * Writing an RoSeq object to an XRaw object.
 */

void _write_RoSeq_to_XRaw(SEXP x, int offset, const RoSeq *seq,
		const int *chrtrtable)
{
	char *dest;

	dest = (char *) RAW(_get_XRaw_tag(x)) + offset;
	_copy_seq(dest, seq->elts, seq->nelt, chrtrtable);
	return;
}


/****************************************************************************
 * Converting an RoSeqs struct (array of arrays of const chars) into a set
 * of sequences.
 * FIXME: This has nothing to do with XRaw objects! Find another place for it
 */

SEXP _new_CHARSXP_from_RoSeq(const RoSeq *seq, SEXP lkup)
{
	// IMPORTANT: We use user-controlled memory for this private memory
	// pool so it is persistent between calls to .Call().
	// It will last until the end of the R session and can only grow
	// during the session. It is NOT a memory leak!
	static int bufsize = 0;
	static char *buf = NULL;
	int new_bufsize;
	char *new_buf;

	new_bufsize = seq->nelt + 1;
	if (new_bufsize > bufsize) {
		new_buf = (char *) realloc(buf, new_bufsize);
		if (new_buf == NULL)
			error("_new_CHARSXP_from_RoSeq(): "
			      "call to realloc() failed");
		buf = new_buf;
		bufsize = new_bufsize;
	}
	if (lkup == R_NilValue) {
		_Biostrings_memcpy_to_i1i2(0, seq->nelt - 1,
			buf, seq->nelt,
			seq->elts, seq->nelt, sizeof(char));
	} else {
		_Biostrings_translate_charcpy_to_i1i2(0, seq->nelt - 1,
			buf, seq->nelt,
			seq->elts, seq->nelt,
			INTEGER(lkup), LENGTH(lkup));
	}
	buf[seq->nelt] = 0;
	return mkChar(buf);
}

SEXP _new_STRSXP_from_RoSeqs(RoSeqs seqs, SEXP lkup)
{
	SEXP ans;
	int i;
	const RoSeq *seq;

	PROTECT(ans = NEW_CHARACTER(seqs.nelt));
	for (i = 0, seq = seqs.elts; i < seqs.nelt; i++, seq++)
		SET_STRING_ELT(ans, i, _new_CHARSXP_from_RoSeq(seq, lkup));
	UNPROTECT(1);
	return ans;
}

