/****************************************************************************
 *                 Basic manipulation of XStringSet objects                 *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"

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

static SEXP get_XStringSet_super(SEXP x)
{
        return GET_SLOT(x, install("super"));
}

const char *_get_XStringSet_baseClass(SEXP x)
{
	return _get_class(get_XStringSet_super(x));
}

int _get_XStringSet_length(SEXP x)
{
	// Because an XStringSet object IS an IRanges object
	return _get_IRanges_length(x);
}

CachedXStringSet _new_CachedXStringSet(SEXP x)
{
	CachedXStringSet cached_x;
	SEXP super, tag;
	int offset;

	cached_x.start = INTEGER(_get_IRanges_start(x));
	cached_x.width = INTEGER(_get_IRanges_width(x));

	super = get_XStringSet_super(x);
	tag = _get_XRaw_tag(_get_XString_xdata(super));
	offset = INTEGER(GET_SLOT(super, install("offset")))[0];
	cached_x.super_elts = (char *) RAW(tag) + offset;
	cached_x.super_nelt = INTEGER(GET_SLOT(super, install("length")))[0];

	cached_x.baseClass = _get_class(super);
	cached_x.enc_chrtrtable = get_enc_chrtrtable(cached_x.baseClass);
	cached_x.dec_chrtrtable = get_dec_chrtrtable(cached_x.baseClass);

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
 * Do NOT try to make this a .Call() entry point!
 * Its arguments are NOT duplicated so it would be a disaster if they were
 * coming from the user space.
 */
static SEXP new_XStringSet_from_IRanges_and_super(SEXP ranges, SEXP super)
{
	char classbuf[80]; // longest string will be "DNAStringSet"
	SEXP class_def, ans;

	snprintf(classbuf, sizeof(classbuf), "%sSet", _get_class(super));
	class_def = MAKE_CLASS(classbuf);
	PROTECT(ans = NEW_OBJECT(class_def));
	_copy_IRanges_slots(ans, ranges);
	SET_SLOT(ans, mkChar("super"), super);
	UNPROTECT(1);
	return ans;
}

/*
 * Assume that the sequences in 'seqs' are NOT already encoded.
 */
SEXP _new_XStringSet_from_RoSeqs(const char *baseClass, RoSeqs seqs)
{
	SEXP ranges, super, ans;

#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _new_XStringSet_from_RoSeqs(): BEGIN\n");
	}
#endif
	PROTECT(ranges = _new_IRanges_from_RoSeqs("LockedIRanges", seqs));
	PROTECT(super = _new_XString_from_RoSeqs(baseClass, seqs));
	PROTECT(ans = new_XStringSet_from_IRanges_and_super(ranges, super));
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _new_XStringSet_from_RoSeqs(): END\n");
	}
#endif
	UNPROTECT(3);
	return ans;
}

/*
 * Does NOT duplicate 'x'. The @NAMES slot is modified in place!
 */
void _set_XStringSet_names(SEXP x, SEXP names)
{
	_set_IRanges_names(x, names);
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
	SEXP ranges, super, ans;

#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _alloc_XStringSet(): BEGIN\n");
		Rprintf("[DEBUG] _alloc_XStringSet(): "
			" baseClass=%s length=%d super_length=%d\n",
			baseClass, length, super_length);
	}
#endif
	PROTECT(ranges = _alloc_IRanges("LockedIRanges", length));
	PROTECT(super = _alloc_XString(baseClass, super_length));
	PROTECT(ans = new_XStringSet_from_IRanges_and_super(ranges, super));
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
	const int *chrtrtable;

	if (i == 0) {
		new_start = 1;
	} else {
		new_start = x->start[i - 1] + x->width[i - 1];
	}
	chrtrtable = encode ? x->enc_chrtrtable : NULL;
	_copy_seq(x->super_elts + new_start - 1, seq->elts, seq->nelt, chrtrtable);
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
 * Low-level manipulation of XStringList objects.
 */

int _get_XStringList_length(SEXP x)
{
	return LENGTH(GET_SLOT(x, install("seqs")));
}

RoSeq _get_XStringList_elt_asRoSeq(SEXP x, int i)
{
	SEXP seqs;

	seqs = GET_SLOT(x, install("seqs"));
	return _get_XString_asRoSeq(VECTOR_ELT(seqs, i));
}

/*
 * --- .Call ENTRY POINT ---
 * 'x_seqs' must be the list, NOT the XStringList object!
 * TODO: make this work directly on the XStringList object and use the
 * 2 helper functions above to simplify the code.
 */
SEXP XStrings_to_nchars(SEXP x_seqs)
{
	SEXP ans;
	int nseq, i, *ans_elt;
	RoSeq seq;

	nseq = LENGTH(x_seqs);
	PROTECT(ans = NEW_INTEGER(nseq));
	for (i = 0, ans_elt = INTEGER(ans); i < nseq; i++, ans_elt++) {
		seq = _get_XString_asRoSeq(VECTOR_ELT(x_seqs, i));
		*ans_elt = seq.nelt;
	}
	UNPROTECT(1);
	return ans;
}

