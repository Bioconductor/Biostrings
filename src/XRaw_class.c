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

SEXP _new_XRaw_from_RoSeqs(const RoSeqs *seqs, SEXP lkup)
{
	SEXP tag, ans;
	int tag_length, i;
	const RoSeq *seq;
	char *dest;

	tag_length = 0;
	for (i = 0, seq = seqs->elts; i < seqs->nelt; i++, seq++)
		tag_length += seq->nelt;
	PROTECT(tag = NEW_RAW(tag_length));
	dest = (char *) RAW(tag);
	for (i = 0, seq = seqs->elts; i < seqs->nelt; i++, seq++) {
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

/*
 * --- .Call ENTRY POINT ---
 * Arguments:
 *   x: a character vector;
 *   start/width: integer vectors of the same length as 'x' and describing a
 *                valid "narrowing" of 'x';
 *   lkup: lookup table for encoding the letters in 'x';
 *   collapse: not yet supported.
 * TODO: Support the 'collapse' argument
 */
SEXP new_XRaw_from_STRSXP(SEXP x, SEXP start, SEXP width,
		SEXP collapse, SEXP lkup)
{
	int nseq;
	RoSeqs seqs;

	nseq = LENGTH(start);
	if (collapse == R_NilValue) {
		if (nseq != 1)
			error("'collapse' must be specified when the number "
			      "of input sequences is not exactly 1");
	} else {
		if (LENGTH(collapse) != 1
		 || LENGTH(STRING_ELT(collapse, 0)) != 0)
			error("'collapse' can only be NULL "
			      "or the empty string for now");
	}
	seqs = _new_RoSeqs_from_STRSXP(nseq, x);
	_narrow_RoSeqs(&seqs, start, width);
	return _new_XRaw_from_RoSeqs(&seqs, lkup);
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

