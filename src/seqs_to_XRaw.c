/****************************************************************************
 *              Turning a set of sequences into an XRaw object              *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"

static int debug = 0;

SEXP Biostrings_debug_seqs_to_XRaw()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'seqs_to_XRaw.c'\n", debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'seqs_to_XRaw.c'\n");
#endif
	return R_NilValue;
}

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

