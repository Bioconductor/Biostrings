/****************************************************************************
 *                 Basic manipulation of XStringSet objects                 *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"
#include "IRanges_interface.h"

static int debug = 0;

SEXP debug_XStringSet_class()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'XStringSet_class.c'\n", debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'XStringSet_class.c'\n");
#endif
	return R_NilValue;
}

SEXP _get_XStringSet_super(SEXP x)
{
	return GET_SLOT(x, install("super"));
}

const char *_get_XStringSet_baseClass(SEXP x)
{
	return get_classname(_get_XStringSet_super(x));
}

SEXP _get_XStringSet_ranges(SEXP x)
{
	return GET_SLOT(x, install("ranges"));
}

int _get_XStringSet_length(SEXP x)
{
	return get_IRanges_length(_get_XStringSet_ranges(x));
}

SEXP _get_XStringSet_width(SEXP x)
{
	return get_IRanges_width(_get_XStringSet_ranges(x));
}

CachedXStringSet _new_CachedXStringSet(SEXP x)
{
	CachedXStringSet cached_x;
	SEXP super, ranges;
	RoSeq seq;

	super = _get_XStringSet_super(x);
	seq = _get_XString_asRoSeq(super);
	// We need to discard the const qualifier
	cached_x.super_elts = (char *) seq.elts;
	cached_x.super_nelt = seq.nelt;

	ranges = _get_XStringSet_ranges(x);
	cached_x.start = INTEGER(get_IRanges_start(ranges));
	cached_x.width = INTEGER(get_IRanges_width(ranges));

	cached_x.baseClass = get_classname(super);
	cached_x.enc_byte2code = get_enc_byte2code(cached_x.baseClass);
	cached_x.dec_byte2code = get_dec_byte2code(cached_x.baseClass);

	return cached_x;
}

RoSeq _get_CachedXStringSet_elt_asRoSeq(CachedXStringSet *x, int i)
{
	RoSeq seq;

	seq.elts = x->super_elts + x->start[i] - 1;
	seq.nelt = x->width[i];
	return seq;
}

RoSeq _get_XStringSet_elt_asRoSeq(SEXP x, int i)
{
	CachedXStringSet cached_x;

	cached_x = _new_CachedXStringSet(x);
	return _get_CachedXStringSet_elt_asRoSeq(&cached_x, i);
}

/*
 * Creating a set of sequences (RoSeqs struct) from an XStringSet object.
 */
RoSeqs _new_RoSeqs_from_XStringSet(int nelt, SEXP x)
{
	RoSeqs seqs;
	CachedXStringSet cached_x;
	RoSeq *elt1;
	int i;

	if (nelt > _get_XStringSet_length(x))
		error("_new_RoSeqs_from_XStringSet(): "
		      "'nelt' must be <= '_get_XStringSet_length(x)'");
	seqs = _alloc_RoSeqs(nelt);
	cached_x = _new_CachedXStringSet(x);
	for (i = 0, elt1 = seqs.elts; i < nelt; i++, elt1++)
		*elt1 = _get_CachedXStringSet_elt_asRoSeq(&cached_x, i);
	return seqs;
}

/*
 * Do NOT try to make this a .Call() entry point!
 * Its arguments are NOT duplicated so it would be a disaster if they were
 * coming from the user space.
 */
SEXP _new_XStringSet(const char *classname, SEXP super, SEXP ranges)
{
	char classname_buf[80]; // longest string will be "DNAStringSet"
	SEXP classdef, ans;

	if (classname == NULL) {
		snprintf(classname_buf, sizeof(classname_buf), "%sSet", get_classname(super));
		classname = classname_buf;
	}
	classdef = MAKE_CLASS(classname);
	PROTECT(ans = NEW_OBJECT(classdef));
	SET_SLOT(ans, mkChar("super"), super);
	SET_SLOT(ans, mkChar("ranges"), ranges);
	UNPROTECT(1);
	return ans;
}

/*
 * Making an XStringSet object from the sequences referenced by a RoSeqs struct.
 * Assume that these sequences are NOT already encoded.
 */
SEXP _new_XStringSet_from_RoSeqs(const char *baseClass, const RoSeqs *seqs)
{
	SEXP super, ranges, ans;

#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _new_XStringSet_from_RoSeqs(): BEGIN\n");
	}
#endif
	PROTECT(super = _new_XString_from_RoSeqs(baseClass, seqs));
	PROTECT(ranges = _new_IRanges_from_RoSeqs("IRanges", seqs));
	PROTECT(ans = _new_XStringSet(NULL, super, ranges));
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _new_XStringSet_from_RoSeqs(): END\n");
	}
#endif
	UNPROTECT(3);
	return ans;
}

/*
 * x@ranges@NAMES is modified in place!
 */
void _set_XStringSet_names(SEXP x, SEXP names)
{
	SEXP ranges;

	ranges = GET_SLOT(x, install("ranges"));
	set_IRanges_names(ranges, names);
	return;
}


/****************************************************************************
 * Utilities for creating an XStringSet instance in 2 steps: first create the
 * skeleton (with junk data in it), then fill it with some character data.
 */

/*
 * Allocate only. The 'ranges' and 'super' slots are not initialized (they
 * contain junk, or zeros).
 */
SEXP _alloc_XStringSet(const char *baseClass, int length, int super_length)
{
	SEXP super, ranges, ans;

#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _alloc_XStringSet(): BEGIN\n");
		Rprintf("[DEBUG] _alloc_XStringSet(): "
			" baseClass=%s length=%d super_length=%d\n",
			baseClass, length, super_length);
	}
#endif
	PROTECT(super = _alloc_XString(baseClass, super_length));
	PROTECT(ranges = alloc_IRanges("IRanges", length));
	PROTECT(ans = _new_XStringSet(NULL, super, ranges));
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _alloc_XStringSet(): END\n");
	}
#endif
	UNPROTECT(3);
	return ans;
}

void _write_RoSeq_to_CachedXStringSet_elt(CachedXStringSet *x, int i,
		const RoSeq *seq, int encode)
{
	int new_start;
	const ByteTrTable *byte2code;

	if (i == 0) {
		new_start = 1;
	} else {
		new_start = x->start[i - 1] + x->width[i - 1];
	}
	byte2code = encode ? x->enc_byte2code : NULL;
	_copy_seq(x->super_elts + new_start - 1, seq->elts, seq->nelt, byte2code);
	x->start[i] = new_start;
	x->width[i] = seq->nelt;
	return;
}

void _write_RoSeq_to_XStringSet_elt(SEXP x, int i, const RoSeq *seq, int encode)
{
	CachedXStringSet cached_x;

	cached_x = _new_CachedXStringSet(x);
	_write_RoSeq_to_CachedXStringSet_elt(&cached_x, i, seq, encode);
	return;
}


/****************************************************************************
 * Coercion.
 */

/*
 * --- .Call ENTRY POINT ---
 */
SEXP XStringSet_as_STRSXP(SEXP x, SEXP lkup)
{
	SEXP ans;
	int x_length, i;
	CachedXStringSet cached_x;
	RoSeq xx;

	x_length = _get_XStringSet_length(x);
	cached_x = _new_CachedXStringSet(x);
	PROTECT(ans = NEW_CHARACTER(x_length));
	for (i = 0; i < x_length; i++) {
		xx = _get_CachedXStringSet_elt_asRoSeq(&cached_x, i);
		SET_STRING_ELT(ans, i, _new_CHARSXP_from_RoSeq(&xx, lkup));
	}
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * Getting the order of an XStringSet object.
 */

/*
 * --- .Call ENTRY POINT ---
 */
SEXP XStringSet_order(SEXP x)
{
	SEXP ans;
	int x_length;
	RoSeqs seqs;

	x_length = _get_XStringSet_length(x);
	seqs = _new_RoSeqs_from_XStringSet(x_length, x);
	PROTECT(ans = NEW_INTEGER(x_length));
	_get_RoSeqs_order(&seqs, INTEGER(ans), 1);
	UNPROTECT(1);
	return ans;
}

/*
 * --- .Call ENTRY POINT ---
 */
SEXP XStringSet_rank(SEXP x)
{
	SEXP ans;
	int x_length;
	RoSeqs seqs;

	x_length = _get_XStringSet_length(x);
	seqs = _new_RoSeqs_from_XStringSet(x_length, x);
	PROTECT(ans = NEW_INTEGER(x_length));
	_get_RoSeqs_rank(&seqs, INTEGER(ans));
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * Getting the duplicates of an XStringSet object.
 */

/*
 * --- .Call ENTRY POINT ---
 */
SEXP XStringSet_duplicated(SEXP x)
{
	SEXP ans;
	int x_length;
	RoSeqs seqs;

	x_length = _get_XStringSet_length(x);
	seqs = _new_RoSeqs_from_XStringSet(x_length, x);
	PROTECT(ans = NEW_LOGICAL(x_length));
	_get_RoSeqs_duplicated(&seqs, INTEGER(ans));
	UNPROTECT(1);
	return ans;
}

SEXP XStringSet_not_duplicated(SEXP x)
{
	SEXP ans;
	int x_length;
	RoSeqs seqs;

	x_length = _get_XStringSet_length(x);
	seqs = _new_RoSeqs_from_XStringSet(x_length, x);
	PROTECT(ans = NEW_LOGICAL(x_length));
	_get_RoSeqs_not_duplicated(&seqs, INTEGER(ans));
	UNPROTECT(1);
	return ans;
}
