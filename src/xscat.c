#include "Biostrings.h"
#include "IRanges_interface.h"
#include <S.h> /* for Salloc() */


/*
 * --- .Call ENTRY POINT ---
 */
SEXP XString_xscat(SEXP args)
{
	error("IMPLEMENT ME!");
	return R_NilValue;
}

/*
 * --- .Call ENTRY POINT ---
 */
SEXP XStringSet_xscat(SEXP args)
{
	CachedXStringSet *cached_args;
	int nargs, *arg_lengths, *ii, ans_length, ans_super_length, write_start,
	    i, j, *start, *width;
	SEXP arg, ans_ranges_start, ans_width, ans_super, ans_ranges, ans;
	const char *ans_classname, *ans_baseClass;
	RoSeq seq;

	nargs = LENGTH(args);
	cached_args = Salloc((long) nargs, CachedXStringSet);
	arg_lengths = Salloc((long) nargs, int);
	ii = Salloc((long) nargs, int);
	for (j = 0; j < nargs; j++) {
		arg = VECTOR_ELT(args, j);
		cached_args[j] = _new_CachedXStringSet(arg);
		arg_lengths[j] = _get_XStringSet_length(arg);
		if (j == 0) {
			ans_length = arg_lengths[j];
			ans_classname = get_classname(arg);
			ans_baseClass = _get_XStringSet_baseClass(arg);
		} else {
			if (arg_lengths[j] > ans_length)
				ans_length = arg_lengths[j];
		}
	}
	PROTECT(ans_ranges_start = NEW_INTEGER(ans_length));
	PROTECT(ans_width = NEW_INTEGER(ans_length));
	ans_super_length = 0;
	for (j = 0; j < nargs; j++)
		ii[j] = 0;
	for (i = 0, start = INTEGER(ans_ranges_start), width = INTEGER(ans_width);
	     i < ans_length;
	     i++, start++, width++)
	{
		*start = ans_super_length + 1;
		*width = 0;
		for (j = 0; j < nargs; j++) {
			if (ii[j] >= arg_lengths[j])
				ii[j] = 0; /* recycle */
			seq = _get_CachedXStringSet_elt_asRoSeq(cached_args + j, ii[j]);
			*width += seq.nelt;
			ii[j]++;
		}
		ans_super_length += *width;
	}
	PROTECT(ans_super = _alloc_XString(ans_baseClass, ans_super_length));
	write_start = 1;
	for (j = 0; j < nargs; j++)
		ii[j] = 0;
	for (i = 0; i < ans_length;  i++) {
		for (j = 0; j < nargs; j++) {
			if (ii[j] >= arg_lengths[j])
				ii[j] = 0; /* recycle */
			seq = _get_CachedXStringSet_elt_asRoSeq(cached_args + j, ii[j]);
			_write_RoSeq_to_XString(ans_super, write_start, &seq, 0);
			write_start += seq.nelt;
			ii[j]++;
		}
	}
	PROTECT(ans_ranges = new_IRanges("IRanges",
				ans_ranges_start,
				ans_width,
				R_NilValue));
	PROTECT(ans = _new_XStringSet(ans_classname, ans_super, ans_ranges));
	UNPROTECT(5);
	return ans;
}

