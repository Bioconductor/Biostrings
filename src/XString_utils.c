/****************************************************************************
 *  Low-level manipulation of XString, XStringSet and XStringList objects   *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"

static int debug = 0;

SEXP Biostrings_debug_XString_utils()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'XString_utils.c'\n", debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'XString_utils.c'\n");
#endif
	return R_NilValue;
}


/****************************************************************************
 * Encoding/decoding XString data.
 */

#define TRTABLE_LENGTH 256

static int DNA_enc_trtable[TRTABLE_LENGTH],
	   DNA_dec_trtable[TRTABLE_LENGTH],
	   RNA_enc_trtable[TRTABLE_LENGTH],
	   RNA_dec_trtable[TRTABLE_LENGTH];

static const int *get_enc_trtable(const char *class)
{
	if (strcmp(class, "DNAString") == 0)
		return DNA_enc_trtable;
	else if (strcmp(class, "RNAString") == 0)
		return RNA_enc_trtable;
	return NULL;
}

static void copy_lkup(const int *lkup1, int len1, int *lkup2, int len2)
{
	int i;

	if (len1 > len2)
		error("Biostrings internal error in copy_lkup(): len1 > len2");
	for (i = 0; i < len1; i++)
		lkup2[i] = lkup1[i];
	for ( ; i < len2; i++)
		lkup2[i] = NA_INTEGER;
	return;
}

SEXP init_DNAlkups(SEXP enc_lkup, SEXP dec_lkup)
{
	copy_lkup(INTEGER(enc_lkup), LENGTH(enc_lkup),
		  DNA_enc_trtable, TRTABLE_LENGTH);
	copy_lkup(INTEGER(dec_lkup), LENGTH(dec_lkup),
		  DNA_dec_trtable, TRTABLE_LENGTH);
	return R_NilValue;
}

char _DNAencode(char c)
{
	int code;

	code = DNA_enc_trtable[(unsigned char) c];
	if (code == NA_INTEGER)
		error("_DNAencode: key %d not in lookup table", (int) c);
	return code;
}

char _DNAdecode(char code)
{
	int c;

	c = DNA_dec_trtable[(unsigned char) code];
	if (c == NA_INTEGER)
		error("_DNAdecode: key %d not in lookup table", (int) code);
	return c;
}

SEXP init_RNAlkups(SEXP enc_lkup, SEXP dec_lkup)
{
	copy_lkup(INTEGER(enc_lkup), LENGTH(enc_lkup),
		  RNA_enc_trtable, TRTABLE_LENGTH);
	copy_lkup(INTEGER(dec_lkup), LENGTH(dec_lkup),
		  RNA_dec_trtable, TRTABLE_LENGTH);
	return R_NilValue;
}

char _RNAencode(char c)
{
	int code;

	code = RNA_enc_trtable[(unsigned char) c];
	if (code == NA_INTEGER)
		error("_RNAencode: key %d not in lookup table", (int) c);
	return code;
}

char _RNAdecode(char code)
{
	int c;

	c = RNA_dec_trtable[(unsigned char) code];
	if (c == NA_INTEGER)
		error("_RNAdecode: key %d not in lookup table", (int) code);
	return (char) c;
}


/****************************************************************************
 * Low-level manipulation of XString objects.
 */

static SEXP get_XString_data(SEXP x)
{
	return GET_SLOT(x, install("data"));
}

const char *_get_XString_charseq(SEXP x, int *length)
{
	SEXP xp;
	int offset;

	xp = GET_SLOT(get_XString_data(x), install("xp"));
	offset = INTEGER(GET_SLOT(x, install("offset")))[0];
	*length = INTEGER(GET_SLOT(x, install("length")))[0];
	return (const char *) (RAW(R_ExternalPtrTag(xp)) + offset);
}

/*
 * NOT a .Call() entry point!
 */
SEXP _new_XString(const char *class, SEXP data, int offset, int length)
{
	SEXP class_def, ans;

	class_def = MAKE_CLASS(class);
	PROTECT(ans = NEW_OBJECT(class_def));
	SET_SLOT(ans, mkChar("data"), data);
	SET_SLOT(ans, mkChar("offset"), ScalarInteger(offset));
	SET_SLOT(ans, mkChar("length"), ScalarInteger(length));
	UNPROTECT(1);
	return ans;
}

SEXP _new_XString_from_RoSeqs(const char *class, RoSeqs seqs)
{
	const int *enc_lkup;
        SEXP lkup, data, ans;

	enc_lkup = get_enc_trtable(class);
	if (enc_lkup == NULL) {
		lkup = R_NilValue;
	} else {
		PROTECT(lkup = NEW_INTEGER(TRTABLE_LENGTH));
		copy_lkup(enc_lkup, TRTABLE_LENGTH,
			  INTEGER(lkup), LENGTH(lkup));
	}
	PROTECT(data = _new_XRaw_from_RoSeqs(seqs, lkup));
	PROTECT(ans = _new_XString(class, data, 0, _get_XRaw_length(data)));
	if (enc_lkup == NULL)
		UNPROTECT(2);
	else
		UNPROTECT(3);
	return ans;
}


/****************************************************************************
 * Utilities for creating an XString instance in 2 steps: first create the
 * skeleton (with junk data in it), then fill it with some character data.
 */

/*
 * Allocate only. The 'data' slot is not initialized (it contains junk,
 * or zeros).
 */
SEXP _alloc_XString(const char *class, int length)
{
	SEXP tag, data, ans;

	PROTECT(tag = NEW_RAW(length));
	PROTECT(data = _new_XRaw(tag));
	PROTECT(ans = _new_XString(class, data, 0, length));
	UNPROTECT(3);
	return ans;
}

void _write_RoSeq_to_XString(SEXP x, int start, RoSeq seq)
{
	int offset;

	offset = INTEGER(GET_SLOT(x, install("offset")))[0];
	_write_RoSeq_to_XRaw(get_XString_data(x), offset + start - 1, seq,
			get_enc_trtable(get_class(x)), TRTABLE_LENGTH);
	return;
}


/****************************************************************************
 * Low-level manipulation of XStringSet objects.
 */

int _get_XStringSet_length(SEXP x)
{
	// Because an XStringSet object IS an .IRanges object
	return _get_IRanges_length(x);
}

static SEXP get_XStringSet_super(SEXP x)
{
        return GET_SLOT(x, install("super"));
}

const char *_get_XStringSet_charseq(SEXP x, int i, int *nchar)
{
	SEXP super;
	int start, super_length;
	const char *super_seq;

	start = _get_IRanges_start0(x)[i];
	*nchar = _get_IRanges_width0(x)[i];
	super = get_XStringSet_super(x);
	super_seq = _get_XString_charseq(super, &super_length);
	return super_seq + start - 1;
}

/*
 * Do NOT try to make this a .Call() entry point!
 * Its arguments are NOT duplicated so it would be a disaster if they were
 * coming from the user space.
 */
static SEXP new_XStringSet_from_IRanges_and_super(SEXP iranges, SEXP super)
{
	char classbuf[80]; // longest string will be "DNAStringSet"
	SEXP class_def, ans, ranges_slot;

	snprintf(classbuf, sizeof(classbuf), "%sSet", get_class(super));
	class_def = MAKE_CLASS(classbuf);
	PROTECT(ans = NEW_OBJECT(class_def));
	PROTECT(ranges_slot = duplicate(GET_SLOT(iranges, install("ranges"))));
	SET_SLOT(ans, mkChar("ranges"), ranges_slot);
	SET_SLOT(ans, mkChar("super"), super);
	UNPROTECT(2);
	return ans;
}

/*
 * Assume that the sequences in 'seqs' are NOT already encoded.
 */
SEXP _new_XStringSet_from_RoSeqs(const char *baseClass, RoSeqs seqs)
{
	SEXP iranges, super, ans;

#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _new_XStringSet_from_RoSeqs(): BEGIN\n");
	}
#endif
	PROTECT(iranges = _new_IRanges_from_RoSeqs(seqs));
	PROTECT(super = _new_XString_from_RoSeqs(baseClass, seqs));
	PROTECT(ans = new_XStringSet_from_IRanges_and_super(iranges, super));
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _new_XStringSet_from_RoSeqs(): END\n");
	}
#endif
	UNPROTECT(3);
	return ans;
}

/*
 * Does NOT duplicate 'x'. The @ranges slot is modified in place!
 */
void _set_XStringSet_names(SEXP x, SEXP names)
{
	SEXP iranges, ranges_slot;

	PROTECT(iranges = _replace_IRanges_names(x, names));
	PROTECT(ranges_slot = duplicate(GET_SLOT(iranges, install("ranges"))));
	SET_SLOT(x, mkChar("ranges"), ranges_slot);
	UNPROTECT(2);
	return;
}


/****************************************************************************
 * Utilities for creating an XStringSet instance in 2 steps: first create the
 * skeleton (with junk data in it), then fill it with some character data.
 */

/*
 * Allocate only. The 'start', 'width' and 'super' slots are not initialized
 * (they contain junk, or zeros).
 */
SEXP _alloc_XStringSet(const char *baseClass, int length, int super_length)
{
	SEXP iranges, super, ans;

#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _alloc_XStringSet(): BEGIN\n");
		Rprintf("[DEBUG] _alloc_XStringSet(): "
			" baseClass=%s length=%d super_length=%d\n",
			baseClass, length, super_length);
	}
#endif
	PROTECT(iranges = _alloc_IRanges(length));
	PROTECT(super = _alloc_XString(baseClass, super_length));
	PROTECT(ans = new_XStringSet_from_IRanges_and_super(iranges, super));
	UNPROTECT(3);
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _alloc_XStringSet(): END\n");
	}
#endif
	return ans;
}

void _write_RoSeq_to_XStringSet_elt(SEXP x, int i, RoSeq seq)
{
	int *start, *width, new_start;

	start = INTEGER(_get_IRanges_start(x)) + i;
	width = INTEGER(_get_IRanges_width(x)) + i;
	if (i == 0)
		new_start = 1;
	else
		new_start = *(start - 1) + *(width - 1);
	_write_RoSeq_to_XString(get_XStringSet_super(x), new_start, seq);
	*start = new_start;
	*width = seq.nelt;
	return;
}


/****************************************************************************
 * Low-level manipulation of XStringList objects.
 */

int _get_XStringList_length(SEXP x)
{
	return LENGTH(GET_SLOT(x, install("seqs")));
}

const char *_get_XStringList_charseq(SEXP x, int i, int *nchar)
{
	SEXP seqs, seq;

	seqs = GET_SLOT(x, install("seqs"));
	seq = VECTOR_ELT(seqs, i);
	return _get_XString_charseq(seq, nchar);
}

/* 'x_seqs' must be the list, NOT the XStringList object!
 * TODO: make this work directly on the XStringList object and use the
 * 2 helper functions above to simplify the code.
 */
SEXP XStrings_to_nchars(SEXP x_seqs)
{
	SEXP ans, x_seq;
	int nseq, i, *seq_length;

	nseq = LENGTH(x_seqs);
	PROTECT(ans = NEW_INTEGER(nseq));
	for (i = 0, seq_length = INTEGER(ans); i < nseq; i++, seq_length++) {
		x_seq = VECTOR_ELT(x_seqs, i);
		_get_XString_charseq(x_seq, seq_length);
	}
	UNPROTECT(1);
	return ans;
}

