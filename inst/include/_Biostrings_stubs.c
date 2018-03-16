#include "Biostrings_interface.h"

#define DEFINE_CCALLABLE_STUB(retT, stubname, Targs, args) \
typedef retT(*__ ## stubname ## _funtype__)Targs; \
retT stubname Targs \
{ \
	static __ ## stubname ## _funtype__ fun = NULL; \
	if (fun == NULL) \
		fun = (__ ## stubname ## _funtype__) R_GetCCallable("Biostrings", "_" #stubname); \
	return fun args; \
}

/*
 * Using the above macro when retT (the returned type) is void will make Sun
 * Studio 12 C compiler unhappy. So we need to use the following macro to
 * handle that case.
 */
#define DEFINE_NOVALUE_CCALLABLE_STUB(stubname, Targs, args) \
typedef void(*__ ## stubname ## _funtype__)Targs; \
void stubname Targs \
{ \
	static __ ## stubname ## _funtype__ fun = NULL; \
	if (fun == NULL) \
		fun = (__ ## stubname ## _funtype__) R_GetCCallable("Biostrings", "_" #stubname); \
	fun args; \
	return; \
}


/*
 * Stubs for callables defined in XString_class.c
 */

DEFINE_CCALLABLE_STUB(char, DNAencode,
	(char c),
	(     c)
)

DEFINE_CCALLABLE_STUB(char, DNAdecode,
	(char code),
	(     code)
)

DEFINE_CCALLABLE_STUB(char, RNAencode,
	(char c),
	(     c)
)

DEFINE_CCALLABLE_STUB(char, RNAdecode,
	(char code),
	(     code)
)

/*
 * Stubs for callables defined in XStringSet_class.c
 */

DEFINE_CCALLABLE_STUB(int, get_XStringSet_length,
	(SEXP x),
	(     x)
)

DEFINE_CCALLABLE_STUB(XStringSet_holder, hold_XStringSet,
	(SEXP x),
	(     x)
)

DEFINE_CCALLABLE_STUB(int, get_length_from_XStringSet_holder,
	(const XStringSet_holder *x_holder),
	(                         x_holder)
)

DEFINE_CCALLABLE_STUB(Chars_holder, get_elt_from_XStringSet_holder,
	(const XStringSet_holder *x_holder, int i),
	(                         x_holder,     i)
)

DEFINE_CCALLABLE_STUB(XStringSet_holder, get_linear_subset_from_XStringSet_holder,
	(const XStringSet_holder *x_holder, int offset, int length),
	(                         x_holder,     offset,     length)
)

DEFINE_NOVALUE_CCALLABLE_STUB(set_XStringSet_names,
	(SEXP x, SEXP names),
	(     x,      names)
)

/*
 * Stubs for callables defined in XStringSetList_class.c
 */

DEFINE_CCALLABLE_STUB(XStringSetList_holder, hold_XStringSetList,
	(SEXP x),
	(     x)
)

DEFINE_CCALLABLE_STUB(int, get_length_from_XStringSetList_holder,
        (const XStringSetList_holder *x_holder),
        (                             x_holder)
)

DEFINE_CCALLABLE_STUB(XStringSet_holder, get_elt_from_XStringSetList_holder,
        (const XStringSetList_holder *x_holder, int i),
        (                             x_holder,     i)
)

/*
 * Stubs for callables defined in match_reporting.c
 */

DEFINE_NOVALUE_CCALLABLE_STUB(init_match_reporting,
	(const char *ms_mode, int nPSpair),
	(            ms_mode,     nPSpair)
)

DEFINE_NOVALUE_CCALLABLE_STUB(set_active_PSpair,
	(int PSpair_id),
	(    PSpair_id)
)

DEFINE_NOVALUE_CCALLABLE_STUB(set_match_shift,
	(int shift),
	(    shift)
)

DEFINE_NOVALUE_CCALLABLE_STUB(report_match,
	(int start, int width),
	(    start,     width)
)

DEFINE_NOVALUE_CCALLABLE_STUB(drop_reported_matches,
	(),
	()
)

DEFINE_CCALLABLE_STUB(int, get_match_count,
	(),
	()
)

DEFINE_CCALLABLE_STUB(SEXP, reported_matches_asSEXP,
	(),
	()
)

/*
 * Stubs for callables defined in MIndex_class.c
 */

DEFINE_CCALLABLE_STUB(MIndex_holder, hold_MIndex,
	(SEXP x),
	(     x)
)

DEFINE_CCALLABLE_STUB(int, get_length_from_MIndex_holder,
	(const MIndex_holder *x_holder),
	(                     x_holder)
)

DEFINE_CCALLABLE_STUB(int, get_width0_elt_from_MIndex_holder,
	(const MIndex_holder *x_holder, int i),
	(                     x_holder,     i)
)

DEFINE_CCALLABLE_STUB(IRanges_holder, get_elt_from_MIndex_holder,
	(const MIndex_holder *x_holder, int i),
	(                     x_holder,     i)
)

/*
 * Stubs for callables defined in match_pattern_boyermoore.c
 */

DEFINE_CCALLABLE_STUB(int, match_pattern_boyermoore,
	(const Chars_holder *P, const Chars_holder *S, int nfirstmatches, int walk_backward),
	(                    P,                     S,     nfirstmatches,     walk_backward)
)

