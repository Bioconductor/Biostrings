/****************************************************************************
 *                     RoSeq/RoSeqs low-level utilities                     *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"
#include "IRanges_interface.h"
#include <S.h> /* for Salloc() */

static int debug = 0;

SEXP debug_RoSeq_utils()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'RoSeq_utils.c'\n",
		debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'RoSeq_utils.c'\n");
#endif
	return R_NilValue;
}

RoSeqs _alloc_RoSeqs(int nelt)
{
	RoSeqs seqs;

	seqs.elts = Salloc((long) nelt, RoSeq);
	seqs.nelt = nelt;
	return seqs;
}


/*****************************************************************************
 * "Narrowing" a RoSeqs struct.
 */

void _narrow_RoSeqs(RoSeqs *seqs, SEXP start, SEXP width)
{
	int i, s, w;
	const int *s_p, *w_p;
	RoSeq *seq;

	if (LENGTH(start) != seqs->nelt || LENGTH(width) != seqs->nelt)
		error("Biostrings internal error in _narrow_RoSeqs(): "
		      "'start' and 'width' must have the same length as 'seqs'");
	for (i = 0, seq = seqs->elts, s_p = INTEGER(start), w_p = INTEGER(width);
	     i < seqs->nelt;
	     i++, seq++, s_p++, w_p++)
	{
		s = *s_p;
		w = *w_p;
		if (s == NA_INTEGER || w == NA_INTEGER)
			error("Biostrings internal error in _narrow_RoSeqs():"
			      "NAs in 'start' or 'width' are not supported");
		s--; // 0-based start (offset)
		if (s < 0 || w < 0 || s + w > seq->nelt)
			error("Biostrings internal error in _narrow_RoSeqs():"
			      "invalid narrowing");
		seq->elts += s;
		seq->nelt = w;
	}
	return;
}


/*****************************************************************************
 * From a RoSeq struct to a character string.
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
		IRanges_memcpy_to_i1i2(0, seq->nelt - 1,
			buf, seq->nelt,
			seq->elts, seq->nelt, sizeof(char));
	} else {
		IRanges_charcpy_to_i1i2_with_lkup(0, seq->nelt - 1,
			buf, seq->nelt,
			seq->elts, seq->nelt,
			INTEGER(lkup), LENGTH(lkup));
	}
	buf[seq->nelt] = 0;
	return mkChar(buf);
}


/*****************************************************************************
 * From a character vector to a RoSeqs struct and vice versa.
 */

RoSeqs _new_RoSeqs_from_STRSXP(int nelt, SEXP x)
{
	RoSeqs seqs;
	RoSeq *elt1;
	SEXP elt2;
	int i;

	if (nelt > LENGTH(x))
		error("_new_RoSeqs_from_STRSXP(): "
		      "'nelt' must be <= 'LENGTH(x)'");
	seqs = _alloc_RoSeqs(nelt);
	for (i = 0, elt1 = seqs.elts; i < nelt; i++, elt1++) {
		elt2 = STRING_ELT(x, i);
		if (elt2 == NA_STRING)
			error("input sequence %d is NA", i+1);
		elt1->elts = CHAR(elt2);
		elt1->nelt = LENGTH(elt2);
	}
	return seqs;
}

SEXP _new_STRSXP_from_RoSeqs(const RoSeqs *seqs, SEXP lkup)
{
	SEXP ans;
	int i;
	const RoSeq *seq;

	PROTECT(ans = NEW_CHARACTER(seqs->nelt));
	for (i = 0, seq = seqs->elts; i < seqs->nelt; i++, seq++)
		SET_STRING_ELT(ans, i, _new_CHARSXP_from_RoSeq(seq, lkup));
	UNPROTECT(1);
	return ans;
}


/*****************************************************************************
 * From a RoSeqs struct to a RawPtr object.
 */

SEXP _new_RawPtr_from_RoSeqs(const RoSeqs *seqs, SEXP lkup)
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
                        IRanges_memcpy_to_i1i2(0, seq->nelt - 1,
                                dest, seq->nelt,
                                seq->elts, seq->nelt, sizeof(char));
                } else {
                        IRanges_charcpy_to_i1i2_with_lkup(0, seq->nelt - 1,
                                dest, seq->nelt,
                                seq->elts, seq->nelt,
                                INTEGER(lkup), LENGTH(lkup));
                }
                dest += seq->nelt;
        }
        PROTECT(ans = new_SequencePtr("RawPtr", tag));
        UNPROTECT(2);
        return ans;
}


/*****************************************************************************
 * From a character vector to a RawPtr object.
 *
 * --- .Call ENTRY POINT ---
 * Arguments:
 *   x: a character vector;
 *   start/width: integer vectors of the same length as 'x' and describing a
 *                valid "narrowing" of 'x';
 *   lkup: lookup table for encoding the letters in 'x';
 *   collapse: not yet supported.
 * TODO: Support the 'collapse' argument
 */

SEXP new_RawPtr_from_STRSXP(SEXP x, SEXP start, SEXP width,
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
        return _new_RawPtr_from_RoSeqs(&seqs, lkup);
}


/*****************************************************************************
 * From a CharAEAE buffer to a RoSeqs struct.
 */

RoSeqs _new_RoSeqs_from_CharAEAE(const CharAEAE *char_aeae)
{
	RoSeqs seqs;
	RoSeq *elt1;
	CharAE *elt2;
	int i;

	seqs = _alloc_RoSeqs(char_aeae->nelt);
	for (i = 0, elt1 = seqs.elts, elt2 = char_aeae->elts;
	     i < char_aeae->nelt;
	     i++, elt1++, elt2++)
	{
		elt1->elts = elt2->elts;
		elt1->nelt = elt2->nelt;
	}
	return seqs;
}


/*****************************************************************************
 * From a RoSeqs struct to an IRanges object.
 * Only the lengths of the sequences holded by RoSeqs are considered.
 */

SEXP _new_IRanges_from_RoSeqs(const char *classname, const RoSeqs *seqs)
{
	const RoSeq *seq;
	SEXP start, width, ans;
	int *start_elt, *width_elt, *start_prev_elt, i;

#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _new_IRanges_from_RoSeqs(): BEGIN\n");
	}
#endif
	seq = seqs->elts;
	PROTECT(start = NEW_INTEGER(seqs->nelt));
	PROTECT(width = NEW_INTEGER(seqs->nelt));
	start_elt = INTEGER(start);
	width_elt = INTEGER(width);
	if (seqs->nelt >= 1) {
		*(start_elt++) = 1;
		*(width_elt++) = seq->nelt;
	}
	if (seqs->nelt >= 2)
		for (i = 1, start_prev_elt = INTEGER(start); i < seqs->nelt; i++) {
			*(start_elt++) = *(start_prev_elt++) + (seq++)->nelt;
			*(width_elt++) = seq->nelt;
		}
	PROTECT(ans = new_IRanges(classname, start, width, R_NilValue));
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _new_IRanges_from_RoSeqs(): END\n");
	}
#endif
	UNPROTECT(3);
	return ans;
}


/****************************************************************************
 * Writing a RoSeq object to a RawPtr object.
 */

void _write_RoSeq_to_RawPtr(SEXP x, int offset, const RoSeq *seq,
		const ByteTrTable *byte2code)
{
        char *dest;

        dest = (char *) RAW(get_SequencePtr_tag(x)) + offset;
        _copy_seq(dest, seq->elts, seq->nelt, byte2code);
        return;
}


/*****************************************************************************
 * Getting the order of a RoSeqs struct.
 *
 * The implementation below tries to optimize performance and memory footprint
 * by using a zero-copy approach. This is achieved at the (modest) cost of
 * using the 'base_seq' static variable. Would that be a problem in a
 * multithreading context?
 */

static int cmp_RoSeq(const void *p1, const void *p2)
{
	const RoSeq *seq1, *seq2;
	int min_nelt, ret;

	seq1 = (const RoSeq *) p1;
	seq2 = (const RoSeq *) p2;
	min_nelt = seq1->nelt <= seq2->nelt ? seq1->nelt : seq2->nelt;
	ret = memcmp(seq1->elts, seq2->elts, min_nelt);
	return ret != 0 ? ret : seq1->nelt - seq2->nelt;
}

static const RoSeq *base_seq;

static int cmp_RoSeq_indices(const void *p1, const void *p2)
{
	int i1, i2;

	i1 = *((const int *) p1);
	i2 = *((const int *) p2);
	return cmp_RoSeq(base_seq + i1, base_seq + i2);
}

void _get_RoSeqs_order(const RoSeqs *seqs, int *order)
{
	int i;

	base_seq = seqs->elts - 1; // because we will sort 1-based indices
	for (i = 0; i < seqs->nelt; i++)
		order[i] = i + 1; // 1-based indices
	qsort(order, seqs->nelt, sizeof(int), cmp_RoSeq_indices);
	return;
}

