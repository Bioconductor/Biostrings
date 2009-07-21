#include "Biostrings.h"
#include "IRanges_interface.h"
#include <S.h> /* for Salloc() */


/*
 * --- .Call ENTRY POINT ---
 * Arguments:
 *   args: a non-empty list of XString objects of the same XString base type
 *         (see R/xsbasetype.R).
 * Note that this function is VERY similar to XStringSet_unlist().
 * Maybe both could be unified under a fast c() for XRaw objects.
 */
SEXP XString_xscat(SEXP args)
{
	int nargs, ans_length, write_start, j;
	SEXP arg, ans;
	const char *ans_classname;
	RoSeq seq;

	nargs = LENGTH(args);
	if (nargs == 0)
		error("XString_xscat(): no input");

	/* 1st pass: determine 'ans_classname' and 'ans_length' */
	for (j = 0; j < nargs; j++) {
		arg = VECTOR_ELT(args, j);
		seq = _get_XString_asRoSeq(arg);
		if (j == 0) {
			ans_length = seq.nelt;
			ans_classname = get_classname(arg);
		} else {
			ans_length += seq.nelt;
		}
	}
	PROTECT(ans = _alloc_XString(ans_classname, ans_length));

	/* 2nd pass: fill 'ans' */
	write_start = 1;
	for (j = 0; j < nargs; j++) {
		arg = VECTOR_ELT(args, j);
		seq = _get_XString_asRoSeq(arg);
		_write_RoSeq_to_XString(ans, write_start, &seq, 0);
		write_start += seq.nelt;
	}
	UNPROTECT(1);
	return ans;
}

/*
 * --- .Call ENTRY POINT ---
 * Arguments:
 *   args: a non-empty list of XStringSet objects of the same XString base
 *         type (see R/xsbasetype.R).
 */
SEXP XStringSet_xscat(SEXP args)
{
	cachedXStringSet *cached_args;
	int nargs, *arg_lengths, *ii, ans_length, write_start, i, j, *start, *width;
	unsigned int ans_super_length;
	SEXP arg, ans_ranges_start, ans_width, ans_super, ans_ranges, ans;
	const char *ans_classname, *ans_xsbaseclassname;
	RoSeq seq;

	nargs = LENGTH(args);
	if (nargs == 0)
		error("XStringSet_xscat(): no input");
	cached_args = Salloc((long) nargs, cachedXStringSet);
	arg_lengths = Salloc((long) nargs, int);
	ii = Salloc((long) nargs, int);

	/* 1st pass: determine 'ans_length', 'ans_classname' and 'ans_xsbaseclassname' */
	for (j = 0; j < nargs; j++) {
		arg = VECTOR_ELT(args, j);
		cached_args[j] = _cache_XStringSet(arg);
		arg_lengths[j] = _get_XStringSet_length(arg);
		if (j == 0) {
			ans_length = arg_lengths[j];
			ans_classname = get_classname(arg);
			ans_xsbaseclassname = _get_XStringSet_xsbaseclassname(arg);
		} else {
			if (arg_lengths[j] > ans_length)
				ans_length = arg_lengths[j];
		}
	}
	PROTECT(ans_ranges_start = NEW_INTEGER(ans_length));
	PROTECT(ans_width = NEW_INTEGER(ans_length));

	/* 2nd pass: fill 'ans_ranges_start' and 'ans_width'
	             and determine 'ans_super_length' */
	ans_super_length = 0U;
	for (j = 0; j < nargs; j++)
		ii[j] = 0;
	for (i = 0, start = INTEGER(ans_ranges_start), width = INTEGER(ans_width);
	     i < ans_length;
	     i++, start++, width++)
	{
		*start = ans_super_length + 1U;
		*width = 0;
		for (j = 0; j < nargs; j++) {
			if (ii[j] >= arg_lengths[j])
				ii[j] = 0; /* recycle */
			seq = _get_cachedXStringSet_elt(cached_args + j, ii[j]);
			*width += seq.nelt;
			ii[j]++;
		}
		ans_super_length += *width;
		if (ans_super_length > INT_MAX)
			error("XStringSet_xscat(): reached the maximum number "
			      "of letters an XStringSet\n  object can hold (%d), "
			      "sorry!", INT_MAX);
	}
	PROTECT(ans_super = _alloc_XString(ans_xsbaseclassname, ans_super_length));

	/* 3rd pass: fill 'ans_super' */
	write_start = 1;
	for (j = 0; j < nargs; j++)
		ii[j] = 0;
	for (i = 0; i < ans_length;  i++) {
		for (j = 0; j < nargs; j++) {
			if (ii[j] >= arg_lengths[j])
				ii[j] = 0; /* recycle */
			seq = _get_cachedXStringSet_elt(cached_args + j, ii[j]);
			_write_RoSeq_to_XString(ans_super, write_start, &seq, 0);
			write_start += seq.nelt;
			ii[j]++;
		}
	}

	/* Put 'ans' pieces together */
	PROTECT(ans_ranges = new_IRanges("IRanges",
				ans_ranges_start,
				ans_width,
				R_NilValue));
	PROTECT(ans = _new_XStringSet(ans_classname, ans_super, ans_ranges));
	UNPROTECT(5);
	return ans;
}

