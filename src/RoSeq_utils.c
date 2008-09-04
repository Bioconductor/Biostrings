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

SEXP _new_IRanges_from_RoSeqs(const char *class, const RoSeqs *seqs)
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
	PROTECT(ans = new_IRanges(class, start, width, R_NilValue));
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _new_IRanges_from_RoSeqs(): END\n");
	}
#endif
	UNPROTECT(3);
	return ans;
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
 * Getting the order of a RoSeqs struct.
 *
 * The implementation below uses a zero-copy approach for optimal performance.
 * This is achieved at the (modest) cost of using the 'base_seq' static
 * variable.
 */

static int cmp_RoSeq(const void *p1, const void *p2)
{
	const RoSeq *seq1, *seq2;
	int min_nelt, ret;

	seq1 = (const RoSeq *) p1;
	seq2 = (const RoSeq *) p2;
	if (seq1->nelt <= seq2->nelt)
		min_nelt = seq1->nelt;
	else
		min_nelt = seq2->nelt;
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

