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
			ans_classname = get_classname(arg);
			ans_length = cached_arg.length;
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
		memcpy((char *) RAW(ans_tag) + tag_offset,
		       cached_arg.seq,
		       cached_arg.length * sizeof(char));
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
	cachedXStringSet *cached_args, cached_ans;
	int nargs, *arg_lengths, *ii, ans_length, i, j, *width;
	SEXP arg, ans_width, ans;
	const char *ans_element_type;
	cachedCharSeq cached_arg_elt, cached_ans_elt;
	char ans_classname[40];  /* longest string should be "DNAStringSet" */

	nargs = LENGTH(args);
	if (nargs == 0)
		error("XStringSet_xscat(): no input");
	cached_args = Salloc((long) nargs, cachedXStringSet);
	arg_lengths = Salloc((long) nargs, int);
	ii = Salloc((long) nargs, int);

	/* 1st pass: determine 'ans_element_type' and 'ans_length' */
	for (j = 0; j < nargs; j++) {
		arg = VECTOR_ELT(args, j);
		cached_args[j] = _cache_XStringSet(arg);
		arg_lengths[j] = _get_XStringSet_length(arg);
		if (j == 0) {
			ans_element_type = _get_XStringSet_xsbaseclassname(arg);
			ans_length = arg_lengths[j];
		} else {
			if (arg_lengths[j] > ans_length)
				ans_length = arg_lengths[j];
		}
	}
	PROTECT(ans_width = NEW_INTEGER(ans_length));

	/* 2nd pass: fill 'ans_width' */
	for (j = 0; j < nargs; j++)
		ii[j] = 0;
	for (i = 0, width = INTEGER(ans_width); i < ans_length; i++, width++)
	{
		*width = 0;
		for (j = 0; j < nargs; j++) {
			if (ii[j] >= arg_lengths[j])
				ii[j] = 0; /* recycle */
			cached_arg_elt = _get_cachedXStringSet_elt(
						cached_args + j, ii[j]);
			*width += cached_arg_elt.length;
			ii[j]++;
		}
	}

	if (snprintf(ans_classname, sizeof(ans_classname),
			"%sSet", ans_element_type) >= sizeof(ans_classname)) {
		UNPROTECT(1);
		error("Biostrings internal error in XStringSet_xscat(): "
		      "'ans_classname' buffer too small");
	}
	PROTECT(ans = alloc_XRawList(ans_classname, ans_element_type,
				     ans_width));

	/* 3rd pass: fill 'ans' */
	cached_ans = cache_XVectorList(ans);
	for (j = 0; j < nargs; j++)
		ii[j] = 0;
	for (i = 0; i < ans_length;  i++) {
		cached_ans_elt = _get_cachedXStringSet_elt(&cached_ans, i);
		cached_ans_elt.length = 0;
		for (j = 0; j < nargs; j++) {
			if (ii[j] >= arg_lengths[j])
				ii[j] = 0; /* recycle */
			cached_arg_elt = _get_cachedXStringSet_elt(
							cached_args + j, ii[j]);
			/* cached_ans_elt->seq is a const char * so we need to
			   cast it to char * in order to write to it */
			memcpy((char *) cached_ans_elt.seq +
			                cached_ans_elt.length,
			       cached_arg_elt.seq,
			       cached_arg_elt.length * sizeof(char));
			cached_ans_elt.length += cached_arg_elt.length;
			ii[j]++;
		}
	}

	UNPROTECT(2);
	return ans;
}

