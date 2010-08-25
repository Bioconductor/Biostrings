/****************************************************************************
 *                  Basic manipulation of XString objects                   *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"
#include "IRanges_interface.h"
#include <stdlib.h>  /* for realloc() */

static int debug = 0;

SEXP debug_XString_class()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in file %s\n",
		debug ? "on" : "off", __FILE__);
#else
	Rprintf("Debug mode not available in file %s\n", __FILE__);
#endif
	return R_NilValue;
}


/****************************************************************************
 * Encoding/decoding XString data.
 */

static ByteTrTable DNA_enc_byte2code, DNA_dec_byte2code,
		 RNA_enc_byte2code, RNA_dec_byte2code;

const ByteTrTable *get_enc_byte2code(const char *classname)
{
	if (strcmp(classname, "DNAString") == 0)
		return (const ByteTrTable *) &DNA_enc_byte2code;
	else if (strcmp(classname, "RNAString") == 0)
		return (const ByteTrTable *) &RNA_enc_byte2code;
	return NULL;
}

const ByteTrTable *get_dec_byte2code(const char *classname)
{
	if (strcmp(classname, "DNAString") == 0)
		return (const ByteTrTable *) &DNA_dec_byte2code;
	else if (strcmp(classname, "RNAString") == 0)
		return (const ByteTrTable *) &RNA_dec_byte2code;
	return NULL;
}

/* --- .Call ENTRY POINT --- */
SEXP init_DNAlkups(SEXP enc_lkup, SEXP dec_lkup)
{
	_init_ByteTrTable_with_lkup(DNA_enc_byte2code, enc_lkup);
	_init_ByteTrTable_with_lkup(DNA_dec_byte2code, dec_lkup);
	return R_NilValue;
}

char _DNAencode(char c)
{
	int code;

	code = DNA_enc_byte2code[(unsigned char) c];
	if (code == NA_INTEGER)
		error("_DNAencode(): key %d not in lookup table", (int) c);
	return code;
}

char _DNAdecode(char code)
{
	int c;

	c = DNA_dec_byte2code[(unsigned char) code];
	if (c == NA_INTEGER)
		error("_DNAdecode(): key %d not in lookup table", (int) code);
	return c;
}

/* --- .Call ENTRY POINT --- */
SEXP init_RNAlkups(SEXP enc_lkup, SEXP dec_lkup)
{
	_init_ByteTrTable_with_lkup(RNA_enc_byte2code, enc_lkup);
	_init_ByteTrTable_with_lkup(RNA_dec_byte2code, dec_lkup);
	return R_NilValue;
}

char _RNAencode(char c)
{
	int code;

	code = RNA_enc_byte2code[(unsigned char) c];
	if (code == NA_INTEGER)
		error("_RNAencode(): key %d not in lookup table", (int) c);
	return code;
}

char _RNAdecode(char code)
{
	int c;

	c = RNA_dec_byte2code[(unsigned char) code];
	if (c == NA_INTEGER)
		error("_RNAdecode(): key %d not in lookup table", (int) code);
	return (char) c;
}


/****************************************************************************
 * From CHARACTER to XString and vice-versa.
 */

void _copy_CHARSXP_to_cachedCharSeq(cachedCharSeq *dest, SEXP src,
		int start_in_src, const int *lkup, int lkup_length)
{
	int i1, i2;

	i1 = start_in_src - 1;
	i2 = i1 + dest->length - 1;
	/* dest.seq is a const char * so we need to cast it
	   to char * before we can write to it */
	Ocopy_bytes_from_i1i2_with_lkup(i1, i2,
			(char *) dest->seq, dest->length,
			CHAR(src), LENGTH(src),
			lkup, lkup_length);
	return;
}

SEXP _new_CHARSXP_from_cachedCharSeq(const cachedCharSeq *x, SEXP lkup)
{
	// IMPORTANT: We use user-controlled memory for this private memory
	// pool so it is persistent between calls to .Call().
	// It will last until the end of the R session and can only grow
	// during the session. It is NOT a memory leak!
	static int buflength = 0;
	static char *buf = NULL;
	int new_buflength;
	char *new_buf;

	if (lkup == R_NilValue)
		return mkCharLen(x->seq, x->length);
	new_buflength = x->length;
	if (new_buflength > buflength) {
		new_buf = (char *) realloc(buf, new_buflength);
		if (new_buf == NULL)
			error("_new_CHARSXP_from_cachedCharSeq(): "
			      "call to realloc() failed");
		buf = new_buf;
		buflength = new_buflength;
	}
	Ocopy_bytes_to_i1i2_with_lkup(0, x->length - 1,
		buf, buflength,
		x->seq, x->length,
		INTEGER(lkup), LENGTH(lkup));
	return mkCharLen(buf, x->length);
}

/* --- .Call ENTRY POINT --- */
SEXP new_XString_from_CHARACTER(SEXP classname,
		SEXP x, SEXP start, SEXP width, SEXP lkup)
{
	SEXP x_elt, ans;
	cachedCharSeq cached_ans;
	const int *lkup0;
	int lkup_length;

	if (LENGTH(x) != 1)
		error("zero or more than one input sequence");
	x_elt = STRING_ELT(x, 0);
	if (x_elt == NA_STRING)
		error("input sequence is NA");
	PROTECT(ans = alloc_XRaw(CHAR(STRING_ELT(classname, 0)),
				 INTEGER(width)[0]));
	cached_ans = cache_XRaw(ans);
	if (lkup == R_NilValue) {
		lkup0 = NULL;
	} else {
		lkup0 = INTEGER(lkup);
		lkup_length = LENGTH(lkup);
	}
	_copy_CHARSXP_to_cachedCharSeq(&cached_ans, x_elt,
			INTEGER(start)[0], lkup0, lkup_length);
	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP new_CHARACTER_from_XString(SEXP x, SEXP lkup)
{
	cachedCharSeq cached_x;
	SEXP ans, ans_elt;

	cached_x = cache_XRaw(x);
	PROTECT(ans = NEW_CHARACTER(1));
	PROTECT(ans_elt = _new_CHARSXP_from_cachedCharSeq(&cached_x, lkup));
	SET_STRING_ELT(ans, 0, ans_elt);
	UNPROTECT(2);
	return ans;
}

