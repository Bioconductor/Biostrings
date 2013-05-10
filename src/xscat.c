#include "Biostrings.h"
#include "IRanges_interface.h"
#include "XVector_interface.h"
#include <S.h> /* for Salloc() */


/*
 * --- .Call ENTRY POINT ---
 * Arguments:
 *   args: a non-empty list of XString objects of the same XString base type
 *         (see R/seqtype.R).
 * Note that this function is VERY similar to XStringSet_unlist().
 * Maybe both could be unified under a fast c() for XRaw objects.
 */
SEXP XString_xscat(SEXP args)
{
	int nargs, ans_length, tag_offset, j;
	SEXP arg, ans_tag, ans;
	const char *ans_classname;
	cachedCharSeq cached_arg;

	nargs = LENGTH(args);
	if (nargs == 0)
		error("XString_xscat(): no input");

	/* 1st pass: determine 'ans_classname' and 'ans_length' */
	for (j = 0; j < nargs; j++) {
		arg = VECTOR_ELT(args, j);
		cached_arg = cache_XRaw(arg);
		if (j == 0) {
			ans_length = cached_arg.length;
			ans_classname = get_classname(arg);
		} else {
			ans_length += cached_arg.length;
		}
	}
	PROTECT(ans_tag = NEW_RAW(ans_length));

	/* 2nd pass: fill 'ans_tag' */
	tag_offset = 0;
	for (j = 0; j < nargs; j++) {
		arg = VECTOR_ELT(args, j);
		cached_arg = cache_XRaw(arg);
		Ocopy_bytes_to_i1i2_with_lkup(tag_offset,
				tag_offset + cached_arg.length - 1,
				(char *) RAW(ans_tag), LENGTH(ans_tag),
				cached_arg.seq, cached_arg.length,
				NULL, 0);
		tag_offset += cached_arg.length;
	}

	/* Make 'ans' */
	PROTECT(ans = new_XRaw_from_tag(ans_classname, ans_tag));
	UNPROTECT(2);
	return ans;
}

/*
 * --- .Call ENTRY POINT ---
 * Arguments:
 *   args: a non-empty list of XStringSet objects of the same XString base
 *         type (see R/seqtype.R).
 */
SEXP XStringSet_xscat(SEXP args)
{
	cachedXStringSet *cached_args;
	int nargs, *arg_lengths, *ii, ans_length, tag_offset, i, j,
	    *start, *width;
	unsigned int ans_tag_length;
	SEXP arg, ans_ranges_start, ans_width, ans_tag, ans_ranges, ans;
	const char *ans_element_type;
	cachedCharSeq cached_arg_elt;
	char ans_classname[40];  /* longest string should be "DNAStringSet" */

	nargs = LENGTH(args);
	if (nargs == 0)
		error("XStringSet_xscat(): no input");
	cached_args = Salloc((long) nargs, cachedXStringSet);
	arg_lengths = Salloc((long) nargs, int);
	ii = Salloc((long) nargs, int);

	/* 1st pass: determine 'ans_length' and 'ans_element_type' */
	for (j = 0; j < nargs; j++) {
		arg = VECTOR_ELT(args, j);
		cached_args[j] = _cache_XStringSet(arg);
		arg_lengths[j] = _get_XStringSet_length(arg);
		if (j == 0) {
			ans_length = arg_lengths[j];
			ans_element_type = _get_XStringSet_xsbaseclassname(arg);
		} else {
			if (arg_lengths[j] > ans_length)
				ans_length = arg_lengths[j];
		}
	}
	PROTECT(ans_ranges_start = NEW_INTEGER(ans_length));
	PROTECT(ans_width = NEW_INTEGER(ans_length));

	/* 2nd pass: fill 'ans_ranges_start' and 'ans_width'
	             and determine 'ans_tag_length' */
	ans_tag_length = 0U;
	for (j = 0; j < nargs; j++)
		ii[j] = 0;
	for (i = 0, start = INTEGER(ans_ranges_start), width = INTEGER(ans_width);
	     i < ans_length;
	     i++, start++, width++)
	{
		*start = ans_tag_length + 1U;
		*width = 0;
		for (j = 0; j < nargs; j++) {
			if (ii[j] >= arg_lengths[j])
				ii[j] = 0; /* recycle */
			cached_arg_elt = _get_cachedXStringSet_elt(
							cached_args + j, ii[j]);
			*width += cached_arg_elt.length;
			ii[j]++;
		}
		ans_tag_length += *width;
		if (ans_tag_length > INT_MAX)
			error("XStringSet_xscat(): reached the maximum number "
			      "of letters an XStringSet\n  object can hold (%d), "
			      "sorry!", INT_MAX);
	}
	PROTECT(ans_tag = NEW_RAW(ans_tag_length));

	/* 3rd pass: fill 'ans_tag' */
	tag_offset = 0;
	for (j = 0; j < nargs; j++)
		ii[j] = 0;
	for (i = 0; i < ans_length;  i++) {
		for (j = 0; j < nargs; j++) {
			if (ii[j] >= arg_lengths[j])
				ii[j] = 0; /* recycle */
			cached_arg_elt = _get_cachedXStringSet_elt(
							cached_args + j, ii[j]);
			Ocopy_bytes_to_i1i2_with_lkup(tag_offset,
				tag_offset + cached_arg_elt.length - 1,
				(char *) RAW(ans_tag), LENGTH(ans_tag),
				cached_arg_elt.seq, cached_arg_elt.length,
				NULL, 0);
			tag_offset += cached_arg_elt.length;
			ii[j]++;
		}
	}

	/* Put 'ans' pieces together */
	if (snprintf(ans_classname, sizeof(ans_classname),
			"%sSet", ans_element_type) >= sizeof(ans_classname)) {
		UNPROTECT(3);
		error("Biostrings internal error in XStringSet_xscat(): "
		      "'ans_classname' buffer too small");
	}
	PROTECT(ans_ranges = new_IRanges("IRanges",
				ans_ranges_start,
				ans_width,
				R_NilValue));
	PROTECT(ans = new_XRawList_from_tag(ans_classname, ans_element_type,
				ans_tag, ans_ranges));
	UNPROTECT(5);
	return ans;
}

