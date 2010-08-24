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

DEFINE_CCALLABLE_STUB(const char *, get_XStringSet_xsbaseclassname,
	(SEXP x),
	(     x)
)

DEFINE_CCALLABLE_STUB(cachedXStringSet, cache_XStringSet,
	(SEXP x),
	(     x)
)

DEFINE_CCALLABLE_STUB(int, get_cachedXStringSet_length,
	(const cachedXStringSet *cached_x),
	(                        cached_x)
)

DEFINE_CCALLABLE_STUB(cachedCharSeq, get_cachedXStringSet_elt,
	(const cachedXStringSet *cached_x, int i),
	(                        cached_x,     i)
)

DEFINE_NOVALUE_CCALLABLE_STUB(set_XStringSet_names,
	(SEXP x, SEXP names),
	(     x,      names)
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

DEFINE_CCALLABLE_STUB(cachedMIndex, cache_MIndex,
	(SEXP x),
	(     x)
)

DEFINE_CCALLABLE_STUB(int, get_cachedMIndex_length,
	(const cachedMIndex *cached_x),
	(                    cached_x)
)

DEFINE_CCALLABLE_STUB(int, get_cachedMIndex_elt_width0,
	(const cachedMIndex *cached_x, int i),
	(                    cached_x,     i)
)

DEFINE_CCALLABLE_STUB(cachedIRanges, get_cachedMIndex_elt,
	(const cachedMIndex *cached_x, int i),
	(                    cached_x,     i)
)

/*
 * Stubs for callables defined in match_pattern_boyermoore.c
 */

DEFINE_CCALLABLE_STUB(int, match_pattern_boyermoore,
	(const cachedCharSeq *P, const cachedCharSeq *S, int nfirstmatches, int walk_backward),
	(                     P,                      S,     nfirstmatches,     walk_backward)
)

