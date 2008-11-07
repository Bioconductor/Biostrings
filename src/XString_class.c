/****************************************************************************
 *                  Basic manipulation of XString objects                   *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"
#include "IRanges_interface.h"

static int debug = 0;

SEXP debug_XString_class()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'XString_class.c'\n", debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'XString_class.c'\n");
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

/*
 * --- .Call ENTRY POINT ---
 */
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

/*
 * --- .Call ENTRY POINT ---
 */
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
 * Low-level manipulation of XString objects.
 */

RoSeq _get_XString_asRoSeq(SEXP x)
{
	RoSeq seq;
	SEXP tag;
	int offset;

	tag = get_SequencePtr_tag(get_XSequence_xdata(x));
	offset = INTEGER(get_XSequence_offset(x))[0];
	seq.elts = (const char *) (RAW(tag) + offset);
	seq.nelt = INTEGER(get_XSequence_length(x))[0];
	return seq;
}

/*
 * Creating a set of sequences (RoSeqs struct) from an XString object.
 */
static RoSeqs new_RoSeqs_from_XString(int nelt, SEXP x)
{
	RoSeqs seqs;
	RoSeq *elt1;
	int i;

	seqs = _alloc_RoSeqs(nelt);
	for (i = 0, elt1 = seqs.elts; i < nelt; i++, elt1++)
		*elt1 = _get_XString_asRoSeq(x);
	return seqs;
}

/*
 * --- .Call ENTRY POINT ---
 * Arguments:
 *   x: an XString object;
 *   start/width: integer vectors of the same length as 'x' and describing a
 *                set of valid ranges in 'x';
 *   lkup: lookup table for (re)encoding the letters in 'x'.
 */
SEXP new_RawPtr_from_XString(SEXP x, SEXP start, SEXP width, SEXP lkup)
{
	int nseq;
	RoSeqs seqs;

	nseq = LENGTH(start);
	seqs = new_RoSeqs_from_XString(nseq, x);
	_narrow_RoSeqs(&seqs, start, width);
	return _new_RawPtr_from_RoSeqs(&seqs, lkup);
}

/*
 * Making an XString object from the sequences referenced by a RoSeqs struct.
 * Assume that these sequences are NOT already encoded.
 */
SEXP _new_XString_from_RoSeqs(const char *classname, const RoSeqs *seqs)
{
	const ByteTrTable *byte2code;
	SEXP lkup, xdata, ans;

	byte2code = get_enc_byte2code(classname);
	PROTECT(lkup = _new_lkup_from_ByteTrTable(byte2code));
	PROTECT(xdata = _new_RawPtr_from_RoSeqs(seqs, lkup));
	PROTECT(ans = new_XSequence(classname, xdata, 0, get_SequencePtr_length(xdata)));
	UNPROTECT(3);
	return ans;
}


/****************************************************************************
 * Utilities for creating an XString instance in 2 steps: first create the
 * skeleton (with junk data in it), then fill it with some character data.
 */

/*
 * Allocate only. The 'xdata' slot is not initialized (it contains junk,
 * or zeros).
 */
SEXP _alloc_XString(const char *classname, int length)
{
	SEXP tag, xdata, ans;

	PROTECT(tag = NEW_RAW(length));
	PROTECT(xdata = new_SequencePtr("RawPtr", tag));
	PROTECT(ans = new_XSequence(classname, xdata, 0, length));
	UNPROTECT(3);
	return ans;
}

void _write_RoSeq_to_XString(SEXP x, int start, const RoSeq *seq, int encode)
{
	int offset;
	const ByteTrTable *enc_byte2code;

	offset = INTEGER(get_XSequence_offset(x))[0];
	enc_byte2code = encode ? get_enc_byte2code(get_classname(x)) : NULL;
	_write_RoSeq_to_RawPtr(get_XSequence_xdata(x), offset + start - 1, seq, enc_byte2code);
	return;
}

