#include "Biostrings_interface.h"

#define DEFINE_CCALLABLE_STUB(retT, stubname, Targs, args) \
retT stubname Targs \
{ \
	static retT (*fun)Targs = NULL; \
	if (fun == NULL) \
		fun = (retT (*)Targs) R_GetCCallable("Biostrings", "_" #stubname); \
	return fun args; \
}


/*
 * Stubs for callables defined in RoSeq_utils.c
 */

DEFINE_CCALLABLE_STUB(SEXP, new_STRSXP_from_RoSeqs,
	(const RoSeqs *seqs, SEXP lkup),
	(              seqs,      lkup)
);

DEFINE_CCALLABLE_STUB(RoSeqs, new_RoSeqs_from_CharAEAE,
	(const CharAEAE *char_aeae),
	(                char_aeae)
);

DEFINE_CCALLABLE_STUB(SEXP, new_IRanges_from_RoSeqs,
	(const char *class, const RoSeqs *seqs),
	(            class,               seqs)
);


/*
 * Stubs for callables defined in XString_class.c
 */

DEFINE_CCALLABLE_STUB(char, DNAencode,
	(char c),
	(     c)
);

DEFINE_CCALLABLE_STUB(char, DNAdecode,
	(char code),
	(     code)
);

DEFINE_CCALLABLE_STUB(char, RNAencode,
	(char c),
	(     c)
);

DEFINE_CCALLABLE_STUB(char, RNAdecode,
	(char code),
	(     code)
);

DEFINE_CCALLABLE_STUB(RoSeq, get_XString_asRoSeq,
	(SEXP x),
	(     x)
);


/*
 * Stubs for callables defined in XStringSet_class.c
 */

DEFINE_CCALLABLE_STUB(const char *, get_XStringSet_baseClass,
	(SEXP x),
	(     x)
);

DEFINE_CCALLABLE_STUB(int, get_XStringSet_length,
	(SEXP x),
	(     x)
);

DEFINE_CCALLABLE_STUB(CachedXStringSet, new_CachedXStringSet,
	(SEXP x),
	(     x)
);

DEFINE_CCALLABLE_STUB(RoSeq, get_CachedXStringSet_elt_asRoSeq,
	(CachedXStringSet *x, int i),
	(                  x,     i)
);

DEFINE_CCALLABLE_STUB(RoSeq, get_XStringSet_elt_asRoSeq,
	(SEXP x, int i),
	(     x,     i)
);

DEFINE_CCALLABLE_STUB(SEXP, new_XStringSet_from_RoSeqs,
	(const char *baseClass, const RoSeqs *seqs),
	(            baseClass,               seqs)
);

DEFINE_CCALLABLE_STUB(void, set_XStringSet_names,
	(SEXP x, SEXP names),
	(     x,      names)
);

DEFINE_CCALLABLE_STUB(SEXP, alloc_XStringSet,
	(const char *baseClass, int length, int super_length),
	(            baseClass,     length,     super_length)
);

DEFINE_CCALLABLE_STUB(void, write_RoSeq_to_CachedXStringSet_elt,
	(CachedXStringSet *x, int i, const RoSeq *seq, int encode),
	(                  x,     i,              seq,     encode)
);

DEFINE_CCALLABLE_STUB(void, write_RoSeq_to_XStringSet_elt,
	(SEXP x, int i, const RoSeq *seq, int encode),
	(     x,     i,              seq,     encode)
);


/*
 * Stubs for callables defined in match_reporting.c
 */

DEFINE_CCALLABLE_STUB(void, init_match_reporting,
	(int mrmode),
	(    mrmode)
);

DEFINE_CCALLABLE_STUB(int, report_match,
	(int start, int end),
	(    start,     end)
);

DEFINE_CCALLABLE_STUB(SEXP, reported_matches_asSEXP,
	(),
	()
);

