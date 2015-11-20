#include "Biostrings.h"
#include "XVector_interface.h"
#include "IRanges_interface.h"
#include "S4Vectors_interface.h"

#include <limits.h>  /* for INT_MAX and INT_MIN */
#include <string.h>  /* for memcpy */
#include <stdlib.h>  /* for malloc and free */


static int compute_length_after_replacements(
		const Chars_holder *x_holder,
		const IRanges_holder *at_holder,
		const XStringSet_holder *value_holder,
		int *nb_replacements,
		int *new_length)
{
	int x_len, at_len, value_len, i, at_start, at_width;
	long long delta;
	Chars_holder value_elt_holder;

	x_len = x_holder->length;
	at_len = get_length_from_IRanges_holder(at_holder);
	value_len = _get_length_from_XStringSet_holder(value_holder);
	if (value_len != at_len)
		return -1;
	*nb_replacements = at_len;
	delta = 0;
	for (i = 0; i < at_len; i++) {
		at_start = get_start_elt_from_IRanges_holder(at_holder, i);
		at_width = get_width_elt_from_IRanges_holder(at_holder, i);
		if (at_start < 1 || at_start + at_width - 1 > x_len)
			return -2;
		value_elt_holder =
			_get_elt_from_XStringSet_holder(value_holder, i);
		delta += value_elt_holder.length - at_width;
	}
	if (delta < INT_MIN)
		*new_length = -1;
	else if (delta > INT_MAX)
		*new_length = NA_INTEGER;
	else
		*new_length = safe_int_add(x_len, (int) delta);
	return 0;
}

typedef struct ranges_order_bufs_t {
	int *start;
	int *width;
	int *order;
} RangesOrderBufs;

static void free_RangesOrderBufs(RangesOrderBufs *bufs)
{
	if (bufs->start != NULL)
		free(bufs->start);
	if (bufs->width != NULL)
		free(bufs->width);
	if (bufs->order != NULL)
		free(bufs->order);
	return;
}

static int alloc_RangesOrderBufs(RangesOrderBufs *bufs, int buflength)
{
	bufs->start = (int *) malloc(sizeof(int) * buflength);
	bufs->width = (int *) malloc(sizeof(int) * buflength);
	bufs->order = (int *) malloc(sizeof(int) * buflength);
	if (bufs->start != NULL && bufs->width != NULL && bufs->order != NULL)
		return 0;
	free_RangesOrderBufs(bufs);
	return -1;
}

static int replace_at(const Chars_holder *x_holder,
		      const IRanges_holder *at_holder,
		      const XStringSet_holder *value_holder,
		      RangesOrderBufs *bufs,
		      char *dest)
{
	int at_len, i, dest_offset, x_offset, k, x_chunk_len;
	Chars_holder value_elt_holder;

	at_len = get_length_from_IRanges_holder(at_holder);
	for (i = 0; i < at_len; i++) {
		bufs->start[i] =
			get_start_elt_from_IRanges_holder(at_holder, i);
		bufs->width[i] =
			get_width_elt_from_IRanges_holder(at_holder, i);
	}
	get_order_of_int_pairs(bufs->start, bufs->width, at_len, 0,
			       bufs->order, 0);
	dest_offset = x_offset = 0;
	for (k = 0; k < at_len; k++) {
		i = bufs->order[k];
		x_chunk_len = bufs->start[i] - x_offset - 1;
		if (x_chunk_len < 0)
			return -1;
		if (x_chunk_len != 0) {
			memcpy(dest + dest_offset, x_holder->ptr + x_offset,
			       sizeof(char) * x_chunk_len);
			dest_offset += x_chunk_len;	
			x_offset += x_chunk_len;
		}
		value_elt_holder =
			_get_elt_from_XStringSet_holder(value_holder, i);
		if (value_elt_holder.length != 0) {
			memcpy(dest + dest_offset, value_elt_holder.ptr,
			       sizeof(char) * value_elt_holder.length);
			dest_offset += value_elt_holder.length;
		}
		x_offset += bufs->width[i];
	}
	x_chunk_len = x_holder->length - x_offset;
	if (x_chunk_len != 0)
		memcpy(dest + dest_offset, x_holder->ptr + x_offset,
		       sizeof(char) * x_chunk_len);
	return 0;
}

/* --- .Call ENTRY POINT ---
 * Args:
 *   x:     An XString object.
 *   at:    An IRanges object.
 *   value: An XStringSet object of the same seqtype as 'x'.
 */
SEXP XString_replaceAt(SEXP x, SEXP at, SEXP value)
{
	Chars_holder x_holder;
	IRanges_holder at_holder;
	XStringSet_holder value_holder;
	int ret_code, nb_replacements, ans_len;

	const char *ans_classname;
	SEXP ans;
	RangesOrderBufs bufs;
	Chars_holder ans_holder;

	/* Compute 'ans_len' and 'nb_replacements' */
	x_holder = hold_XRaw(x);
	at_holder = hold_IRanges(at);
	value_holder = _hold_XStringSet(value);
	ret_code = compute_length_after_replacements(
				&x_holder,
				&at_holder,
				&value_holder,
				&nb_replacements,
				&ans_len);
	if (ret_code == -1)
		error("'at' and 'value' must have the same length");
	if (ret_code == -2)
		error("some ranges in 'at' are off-limits "
		      "with respect to sequence 'x'");
	if (ans_len == NA_INTEGER)
		error("replacements in 'x' will produce a "
		      "sequence that is too long\n  (i.e. with more "
		      "than '.Machine$integer.max' letters)");
	if (ans_len < 0)
		error("'at' must contain disjoint ranges (see '?isDisjoint')");

	/* Allocate 'ans' and 'bufs' */
	ans_classname = get_classname(x);
	PROTECT(ans = alloc_XRaw(ans_classname, ans_len));
	ret_code = alloc_RangesOrderBufs(&bufs, nb_replacements);
	if (ret_code == -1) {
		UNPROTECT(1);
		error("Biostrings internal error in "
		      "XString_replaceAt():\n\n  "
		      "    memory allocation failed");
	}

	/* Fill 'ans' */
	ans_holder = hold_XRaw(ans);
	ret_code = replace_at(&x_holder,
			      &at_holder,
			      &value_holder,
			      &bufs,
			      (char *) ans_holder.ptr);
	free_RangesOrderBufs(&bufs);
	UNPROTECT(1);
	if (ret_code == -1)
		error("'at' must contain disjoint ranges (see '?isDisjoint')");
	return ans;
}

/* --- .Call ENTRY POINT ---
 * Args:
 *   x:     An XStringSet object.
 *   at:    A CompressedIRangesList object.
 *   value: An XStringSetList object of the same seqtype as 'x'.
 */
SEXP XStringSet_replaceAt(SEXP x, SEXP at, SEXP value)
{
	XStringSet_holder x_holder;
	CompressedIRangesList_holder at_holder;
	XStringSetList_holder value_holder;
	int x_len, at_len, value_len, max_replacements,
	    i, ret_code = 0, nb_replacements, ans_width_elt;
	SEXP ans_width;
	Chars_holder x_elt_holder;
	IRanges_holder at_elt_holder;
	XStringSet_holder value_elt_holder;

	const char *ans_classname, *ans_elt_type;
	SEXP ans;
	RangesOrderBufs bufs;
	XStringSet_holder ans_holder;
	Chars_holder ans_elt_holder;

	x_holder = _hold_XStringSet(x);
	at_holder = hold_CompressedIRangesList(at);
	value_holder = _hold_XStringSetList(value);
	x_len = _get_length_from_XStringSet_holder(&x_holder);
	at_len = get_length_from_CompressedIRangesList_holder(&at_holder);
	value_len = _get_length_from_XStringSetList_holder(&value_holder);
	if (at_len != x_len || value_len != x_len)
		error("'x', 'at' and 'value' must have the same length");

	/* 1st pass: compute 'ans_width' and 'max_replacements' */
	PROTECT(ans_width = NEW_INTEGER(x_len));
	max_replacements = 0;
	for (i = 0; i < x_len; i++) {
		x_elt_holder = _get_elt_from_XStringSet_holder(
					&x_holder, i);
		at_elt_holder = get_elt_from_CompressedIRangesList_holder(
					&at_holder, i);
		value_elt_holder = _get_elt_from_XStringSetList_holder(
					&value_holder, i);
		ret_code = compute_length_after_replacements(
					&x_elt_holder,
					&at_elt_holder,
					&value_elt_holder,
					&nb_replacements,
					&ans_width_elt);
		if (ret_code == -1) {
			UNPROTECT(1);
			error("'at[[%d]]' and 'value[[%d]]' don't have "
			      "the same length. 'at' and 'value'\n  must "
			      "have the same shape, that is, they must have "
			      "the same length and\n  'at[[i]]' and "
			      "'value[[i]]' must have the same length for "
			      "all 'i'.", i + 1, i + 1);
		}
		if (ret_code == -2) {
			UNPROTECT(1);
			error("some ranges in 'at[[%d]]' are off-limits "
			      "with respect to sequence 'x[[%d]]'",
			      i + 1, i + 1);
		}
		if (ans_width_elt == NA_INTEGER) {
			UNPROTECT(1);
			error("replacements in 'x[[%d]]' will produce a "
			      "sequence that is too long\n  (i.e. with more "
			      "than '.Machine$integer.max' letters)", i + 1);
		}
		if (ans_width_elt < 0) {
			UNPROTECT(1);
			error("'at[[%d]]' must contain disjoint ranges "
			      "(see '?isDisjoint')", i + 1);
		}
		INTEGER(ans_width)[i] = ans_width_elt;
		if (nb_replacements > max_replacements)
			max_replacements = nb_replacements;
	}

	/* Allocate 'ans' and 'bufs' */
	ans_classname = get_classname(x);
	ans_elt_type = _get_XStringSet_xsbaseclassname(x);
	PROTECT(ans = alloc_XRawList(ans_classname, ans_elt_type, ans_width));
	ret_code = alloc_RangesOrderBufs(&bufs, max_replacements);
	if (ret_code == -1) {
		UNPROTECT(2);
		error("Biostrings internal error in "
		      "XStringSet_replaceAt():\n\n  "
		      "    memory allocation failed");
	}

	/* 2nd pass: fill 'ans' */
	ans_holder = _hold_XStringSet(ans);
	for (i = 0; i < x_len; i++) {
		ans_elt_holder = _get_elt_from_XStringSet_holder(
					&ans_holder, i);
		x_elt_holder = _get_elt_from_XStringSet_holder(
					&x_holder, i);
		at_elt_holder = get_elt_from_CompressedIRangesList_holder(
					&at_holder, i);
		value_elt_holder = _get_elt_from_XStringSetList_holder(
					&value_holder, i);
		ret_code = replace_at(&x_elt_holder,
				      &at_elt_holder,
				      &value_elt_holder,
				      &bufs,
				      (char *) ans_elt_holder.ptr);
		if (ret_code == -1)
			break;
	}

	free_RangesOrderBufs(&bufs);
	UNPROTECT(2);
	if (ret_code == -1)
		error("'at[[%d]]' must contain disjoint ranges "
		      "(see '?isDisjoint')", i + 1);
	return ans;
}

