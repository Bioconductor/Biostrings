/****************************************************************************
 *             Low-level manipulation of XStringSetList objects             *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"
#include "IRanges_interface.h"
#include "S4Vectors_interface.h"


/****************************************************************************
 * C-level abstract getters.
 */

XStringSetList_holder _hold_XStringSetList(SEXP x)
{
	XStringSetList_holder x_holder;
	SEXP x_end;

	x_holder.classname = get_classname(x);
	x_end = get_PartitioningByEnd_end(get_CompressedList_partitioning(x));
	x_holder.length = LENGTH(x_end);
	x_holder.end = INTEGER(x_end);
	x_holder.unlistData_holder = _hold_XStringSet(
			get_CompressedList_unlistData(x));
	return x_holder;
}

int _get_length_from_XStringSetList_holder(
		const XStringSetList_holder *x_holder)
{
	return x_holder->length;
}

XStringSet_holder _get_elt_from_XStringSetList_holder(
		const XStringSetList_holder *x_holder, int i)
{
	int offset, length;

	offset = i == 0 ? 0 : x_holder->end[i - 1];
	length = x_holder->end[i] - offset;
	return _get_linear_subset_from_XStringSet_holder(
			&(x_holder->unlistData_holder),
			offset, length);
}

