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


/****************************************************************************
 * C-level slot getters.
 *
 * Be careful that these functions do NOT duplicate the returned slot.
 * Thus they cannot be made .Call() entry points!
 */

static SEXP
	super_symbol = NULL,
	ranges_symbol = NULL;

SEXP _get_XStringSet_super(SEXP x)
{
	INIT_STATIC_SYMBOL(super)
	return GET_SLOT(x, super_symbol);
}

SEXP _get_XStringSet_ranges(SEXP x)
{
	INIT_STATIC_SYMBOL(ranges)
	return GET_SLOT(x, ranges_symbol);
}

/* Not strict "slot getters" but very much like. */

int _get_XStringSet_length(SEXP x)
{
	return get_IRanges_length(_get_XStringSet_ranges(x));
}

SEXP _get_XStringSet_width(SEXP x)
{
	return get_IRanges_width(_get_XStringSet_ranges(x));
}

// FIXME: '_get_XStringSet_xsbaseclassname()' won't always return the
// xsbaseclassname e.g. if 'super' belongs to a class that extends
// DNAString (there is no such thing yet).
const char *_get_XStringSet_xsbaseclassname(SEXP x)
{
	return get_classname(_get_XStringSet_super(x));
}


/****************************************************************************
 * C-level abstract getters.
 */

cachedXStringSet _cache_XStringSet(SEXP x)
{
	cachedXStringSet cached_x;
	SEXP super, ranges;

	cached_x.classname = get_classname(x);

	super = _get_XStringSet_super(x);
	// FIXME: 'get_classname(super)' won't always return the
	// xsbaseclassname e.g. if 'super' belongs to a class that extends
	// DNAString (there is no such thing yet).
	cached_x.xsbaseclassname = get_classname(super);
	cached_x.super = cache_XRaw(super);

	ranges = _get_XStringSet_ranges(x);
	cached_x.ranges = cache_IRanges(ranges);
	
	cached_x.enc_byte2code = get_enc_byte2code(cached_x.xsbaseclassname);
	cached_x.dec_byte2code = get_dec_byte2code(cached_x.xsbaseclassname);
	return cached_x;
}

int _get_cachedXStringSet_length(const cachedXStringSet *cached_x)
{
	return get_cachedIRanges_length(&(cached_x->ranges));
}

cachedCharSeq _get_cachedXStringSet_elt(const cachedXStringSet *cached_x, int i)
{
	cachedCharSeq charseq;

	charseq.seq = cached_x->super.seq +
		      get_cachedIRanges_elt_start(&(cached_x->ranges), i) - 1;
	charseq.length = get_cachedIRanges_elt_width(&(cached_x->ranges), i);
	return charseq;
}


/****************************************************************************
 * C-level slot setters.
 *
 * Be careful that these functions do NOT duplicate the assigned value!
 */

static void set_XStringSet_super(SEXP x, SEXP value)
{
	INIT_STATIC_SYMBOL(super)
	SET_SLOT(x, super_symbol, value);
	return;
}

static void set_XStringSet_ranges(SEXP x, SEXP value)
{
	INIT_STATIC_SYMBOL(ranges)
	SET_SLOT(x, ranges_symbol, value);
	return;
}

/* x@ranges@NAMES is modified in-place! */
void _set_XStringSet_names(SEXP x, SEXP names)
{
	set_IRanges_names(_get_XStringSet_ranges(x), names);
	return;
}


/****************************************************************************
 * C-level constructors.
 */

/* Be careful that this constructor does NOT duplicate its arguments before
   putting them in the slots of the returned object.
   So don't try to make it a .Call() entry point! */
SEXP _new_XStringSet(const char *classname, SEXP super, SEXP ranges)
{
	char classname_buf[80]; // longest string will be "DNAStringSet"
	SEXP classdef, ans;

	if (classname == NULL) {
		snprintf(classname_buf, sizeof(classname_buf),
			"%sSet", get_classname(super));
		classname = classname_buf;
	}
	PROTECT(classdef = MAKE_CLASS(classname));
	PROTECT(ans = NEW_OBJECT(classdef));
	set_XStringSet_super(ans, super);
	set_XStringSet_ranges(ans, ranges);
	UNPROTECT(2);
	return ans;
}

/* Allocation WITHOUT initialization.
   The 'ranges' and 'super' slots are not initialized (they contain junk). */
SEXP _alloc_XStringSet(const char *xsbaseclassname, int length, int super_length)
{
	SEXP super, ranges, ans;

	PROTECT(super = _alloc_XString(xsbaseclassname, super_length));
	PROTECT(ranges = alloc_IRanges("IRanges", length));
	PROTECT(ans = _new_XStringSet(NULL, super, ranges));
	UNPROTECT(3);
	return ans;
}

/* Making an XStringSet object from the sequences referenced by a RoSeqs struct.
   Assume that these sequences are NOT already encoded. */
SEXP _new_XStringSet_from_RoSeqs(const char *xsbaseclassname, const RoSeqs *seqs)
{
	SEXP super, ranges, ans;

	PROTECT(super = _new_XString_from_RoSeqs(xsbaseclassname, seqs));
	PROTECT(ranges = _new_IRanges_from_RoSeqs("IRanges", seqs));
	PROTECT(ans = _new_XStringSet(NULL, super, ranges));
	UNPROTECT(3);
	return ans;
}


/****************************************************************************
 * Creating a set of sequences (RoSeqs struct) from an XStringSet object.
 */

RoSeqs _new_RoSeqs_from_XStringSet(int nelt, SEXP x)
{
	RoSeqs seqs;
	cachedXStringSet cached_x;
	cachedCharSeq *elt1;
	int i;

	if (nelt > _get_XStringSet_length(x))
		error("_new_RoSeqs_from_XStringSet(): "
		      "'nelt' must be <= '_get_XStringSet_length(x)'");
	seqs = _alloc_RoSeqs(nelt);
	cached_x = _cache_XStringSet(x);
	for (i = 0, elt1 = seqs.elts; i < nelt; i++, elt1++)
		*elt1 = _get_cachedXStringSet_elt(&cached_x, i);
	return seqs;
}


/****************************************************************************
 * Coercion.
 */

/* Note that XStringSet_unlist() is VERY similar to XString_xscat().
   Maybe both could be unified under a fast c() for XRaw objects. */
SEXP XStringSet_unlist(SEXP x)
{
	SEXP ans;
	int x_length, ans_length, write_start, i;
	cachedXStringSet cached_x;
	cachedCharSeq xx;

	cached_x = _cache_XStringSet(x);
	x_length = _get_cachedXStringSet_length(&cached_x);

	/* 1st pass: determine 'ans_length' */
	ans_length = 0;
	for (i = 0; i < x_length; i++) {
		xx = _get_cachedXStringSet_elt(&cached_x, i);
		ans_length += xx.length;
	}
	PROTECT(ans = _alloc_XString(_get_XStringSet_xsbaseclassname(x), ans_length));

	/* 2nd pass: fill 'ans' */
	write_start = 1;
	for (i = 0; i < x_length; i++) {
		xx = _get_cachedXStringSet_elt(&cached_x, i);
		_write_RoSeq_to_XString(ans, write_start, &xx, 0);
		write_start += xx.length;
	}
	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP XStringSet_as_STRSXP(SEXP x, SEXP lkup)
{
	SEXP ans;
	int x_length, i;
	cachedXStringSet cached_x;
	cachedCharSeq xx;

	cached_x = _cache_XStringSet(x);
	x_length = _get_cachedXStringSet_length(&cached_x);
	PROTECT(ans = NEW_CHARACTER(x_length));
	for (i = 0; i < x_length; i++) {
		xx = _get_cachedXStringSet_elt(&cached_x, i);
		SET_STRING_ELT(ans, i, _new_CHARSXP_from_RoSeq(&xx, lkup));
	}
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * Getting the order or rank of an XStringSet object.
 */

/* --- .Call ENTRY POINT --- */
SEXP XStringSet_is_unsorted(SEXP x, SEXP strictly)
{
	SEXP ans;
	int x_length;
	RoSeqs seqs;

	x_length = _get_XStringSet_length(x);
	seqs = _new_RoSeqs_from_XStringSet(x_length, x);
	PROTECT(ans = NEW_LOGICAL(1));
	LOGICAL(ans)[0] = _get_RoSeqs_is_unsorted(&seqs, LOGICAL(strictly)[0]);
	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT --- */
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

/* --- .Call ENTRY POINT --- */
SEXP XStringSet_rank(SEXP x)
{
	SEXP ans;
	int x_length, *order;
	RoSeqs seqs;

	x_length = _get_XStringSet_length(x);
	seqs = _new_RoSeqs_from_XStringSet(x_length, x);
	order = (int *) R_alloc(seqs.nelt, sizeof(int));
	_get_RoSeqs_order(&seqs, order, 0);
	PROTECT(ans = NEW_INTEGER(x_length));
	_get_RoSeqs_rank(&seqs, order, INTEGER(ans));
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * Getting the duplicates of an XStringSet object.
 */

/* --- .Call ENTRY POINT --- */
SEXP XStringSet_duplicated(SEXP x)
{
	SEXP ans;
	int x_length, *order;
	RoSeqs seqs;

	x_length = _get_XStringSet_length(x);
	seqs = _new_RoSeqs_from_XStringSet(x_length, x);
	order = (int *) R_alloc(seqs.nelt, sizeof(int));
	_get_RoSeqs_order(&seqs, order, 0);
	PROTECT(ans = NEW_LOGICAL(x_length));
	_get_RoSeqs_duplicated(&seqs, order, INTEGER(ans));
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * Identical value matching within an XStringSet object.
 */

/* --- .Call ENTRY POINT --- */
SEXP XStringSet_match(SEXP x, SEXP table, SEXP nomatch)
{
	SEXP ans;
	int x_length, table_length, *seqs_order, *set_order, *match_buffer;
	RoSeqs seqs, set;

	x_length = _get_XStringSet_length(x);
	table_length = _get_XStringSet_length(table);
	seqs = _new_RoSeqs_from_XStringSet(x_length, x);
	seqs_order = (int *) R_alloc(seqs.nelt, sizeof(int));
	_get_RoSeqs_order(&seqs, seqs_order, 0);
	set = _new_RoSeqs_from_XStringSet(table_length, table);
	set_order = (int *) R_alloc(set.nelt, sizeof(int));
	_get_RoSeqs_order(&set, set_order, 0);
	match_buffer = (int *) R_alloc(set.nelt, sizeof(int));
	PROTECT(ans = NEW_INTEGER(x_length));
	_get_RoSeqs_match(&seqs, &set, INTEGER(nomatch)[0], seqs_order, set_order,
			          match_buffer, INTEGER(ans));
	UNPROTECT(1);
	return ans;
}

