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
	Rprintf("Debug mode turned %s in file %s\n",
		debug ? "on" : "off", __FILE__);
#else
	Rprintf("Debug mode not available in file %s\n", __FILE__);
#endif
	return R_NilValue;
}


/****************************************************************************
 * C-level slot getters.
 *
 * Be careful that these functions do NOT duplicate the returned slot.
 * Thus they cannot be made .Call() entry points!
 */

/* Not strict "slot getters" but very much like. */

int _get_XStringSet_length(SEXP x)
{
	return get_XVectorList_length(x);
}

SEXP _get_XStringSet_width(SEXP x)
{
	return get_XVectorList_width(x);
}

const char *_get_XStringSet_xsbaseclassname(SEXP x)
{
	return get_Sequence_elementType(x);
}


/****************************************************************************
 * C-level abstract getters.
 */

cachedXStringSet _cache_XStringSet(SEXP x)
{
	return cache_XVectorList(x);
}

int _get_cachedXStringSet_length(const cachedXStringSet *cached_x)
{
	return get_cachedXVectorList_length(cached_x);
}

cachedCharSeq _get_cachedXStringSet_elt(const cachedXStringSet *cached_x, int i)
{
	return get_cachedXRawList_elt(cached_x, i);
}


/****************************************************************************
 * C-level slot setters.
 *
 */

/* WARNING: x@ranges@NAMES is modified in-place! */
void _set_XStringSet_names(SEXP x, SEXP names)
{
	set_XVectorList_names(x, names);
	return;
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
	SEXP ans_tag, ans;
	int x_length, ans_length, tag_offset, i;
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
	PROTECT(ans_tag = NEW_RAW(ans_length));

	/* 2nd pass: fill 'ans' */
	tag_offset = 0;
	for (i = 0; i < x_length; i++) {
		xx = _get_cachedXStringSet_elt(&cached_x, i);
		Ocopy_bytes_to_i1i2_with_lkup(tag_offset,
				tag_offset + xx.length - 1,
                                (char *) RAW(ans_tag), LENGTH(ans_tag),
                                xx.seq, xx.length,
                                NULL, 0);
		tag_offset += xx.length;
	}

	/* Make 'ans' */
	PROTECT(ans = new_XRaw_from_tag(_get_XStringSet_xsbaseclassname(x),
					ans_tag));
	UNPROTECT(2);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP XStringSet_as_STRSXP(SEXP x, SEXP lkup)
{
	SEXP ans, ans_elt;
	int x_length, i;
	cachedXStringSet cached_x;
	cachedCharSeq xx;

	cached_x = _cache_XStringSet(x);
	x_length = _get_cachedXStringSet_length(&cached_x);
	PROTECT(ans = NEW_CHARACTER(x_length));
	for (i = 0; i < x_length; i++) {
		xx = _get_cachedXStringSet_elt(&cached_x, i);
		PROTECT(ans_elt = _new_CHARSXP_from_cachedCharSeq(&xx, lkup));
		SET_STRING_ELT(ans, i, ans_elt);
		UNPROTECT(1);
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

