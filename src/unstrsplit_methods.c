#include "Biostrings.h"
#include "IRanges_interface.h"
#include "XVector_interface.h"


static int compute_joined_strings_length(const XStringSet_holder *x_holder,
					 int sep_len)
{
	int x_len, joined_length, i;
	Chars_holder x_elt_holder;

	x_len = _get_length_from_XStringSet_holder(x_holder);
	joined_length = 0;
	if (x_len != 0) {
		for (i = 0; i < x_len; i++) {
			x_elt_holder = _get_elt_from_XStringSet_holder(
							x_holder, i);
			joined_length += x_elt_holder.length;
		}
		joined_length += (x_len - 1) * sep_len;
	}
	return joined_length;
}

static void join_strings_in_buf(char *dest, const XStringSet_holder *x_holder,
				const char *sep, int sep_len)
{
	int x_len, i;
	Chars_holder x_elt_holder;

	x_len = _get_length_from_XStringSet_holder(x_holder);
	for (i = 0; i < x_len; i++) {
		if (i != 0) {
			memcpy(dest, sep, sep_len);
			dest += sep_len;
		}
		x_elt_holder = _get_elt_from_XStringSet_holder(x_holder, i);
		memcpy(dest, x_elt_holder.seq, x_elt_holder.length);
		dest += x_elt_holder.length;
	}
	return;
}

/* --- .Call ENTRY POINT ---
 * Args:
 *   x:   An XStringSetList object.
 *   sep: An XString object of the same seqtype as 'x'.
 * Returns an XStringSet object parallel to and of the same seqtype as 'x'.
 */
SEXP XStringSetList_unstrsplit(SEXP x, SEXP sep, SEXP seqtype)
{
	XStringSetList_holder x_holder;
	XStringSet_holder x_elt_holder, ans_holder;
	Chars_holder sep_holder, ans_elt_holder;
	int x_len, sep_len, i;
	const char *seqtype0;
	char ans_elt_type[37];  /* longest string should be "DNAString" */
	char ans_classname[40];  /* longest string should be "DNAStringSet" */
	SEXP ans, ans_width, ans_names;

	x_holder = _hold_XStringSetList(x);
	x_len = _get_length_from_XStringSetList_holder(&x_holder);
	sep_holder = hold_XRaw(sep);
	sep_len = sep_holder.length;
	seqtype0 = CHAR(STRING_ELT(seqtype, 0));
	if (snprintf(ans_elt_type, sizeof(ans_elt_type),
		     "%sString", seqtype0) >= sizeof(ans_elt_type))
		error("Biostrings internal error in "
		      "XStringSetList_unstrsplit(): "
		      "'ans_elt_type' buffer too small");
	if (snprintf(ans_classname, sizeof(ans_classname),
		     "%sSet", ans_elt_type) >= sizeof(ans_classname))
		error("Biostrings internal error in "
		      "XStringSetList_unstrsplit(): "
		      "'ans_classname' buffer too small");

	/* 1st pass: compute 'ans_width' */
	PROTECT(ans_width = NEW_INTEGER(x_len));
	for (i = 0; i < x_len; i++) {
		x_elt_holder = _get_elt_from_XStringSetList_holder(
					&x_holder, i);
		INTEGER(ans_width)[i] = compute_joined_strings_length(
					&x_elt_holder, sep_len);
	}

	/* Allocate 'ans' */
	PROTECT(ans = alloc_XRawList(ans_classname, ans_elt_type, ans_width));

	/* 2nd pass: fill 'ans' */
	ans_holder = _hold_XStringSet(ans);
	for (i = 0; i < x_len; i++) {
		x_elt_holder = _get_elt_from_XStringSetList_holder(
					&x_holder, i);
		ans_elt_holder = _get_elt_from_XStringSet_holder(
					&ans_holder, i);
		join_strings_in_buf((char *) ans_elt_holder.seq, &x_elt_holder,
					sep_holder.seq, sep_holder.length);
	}

	PROTECT(ans_names = duplicate(get_CompressedList_names(x)));
	_set_XStringSet_names(ans, ans_names);
	UNPROTECT(3);
	return ans;
}

