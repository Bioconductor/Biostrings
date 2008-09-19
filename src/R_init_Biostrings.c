#include "Biostrings.h"

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

#define REGISTER_CCALLABLE(fun) \
	R_RegisterCCallable("Biostrings", #fun, (DL_FUNC) &fun)

static const R_CallMethodDef callMethods[] = {

/* copy_seq.c */
	CALLMETHOD_DEF(debug_copy_seq, 0),

/* utils.c */
	CALLMETHOD_DEF(debug_utils, 0),

/* RoSeq_utils.c */
	CALLMETHOD_DEF(debug_RoSeq_utils, 0),
	CALLMETHOD_DEF(new_RawPtr_from_STRSXP, 5),

/* Dups_utils.c */
	CALLMETHOD_DEF(debug_Dups_utils, 0),
	CALLMETHOD_DEF(Dups_diff, 2),

/* SparseList_utils.c */
	CALLMETHOD_DEF(debug_SparseList_utils, 0),

/* fasta_io.c */
	CALLMETHOD_DEF(debug_fasta_io, 0),
	CALLMETHOD_DEF(RawPtr_loadFASTA, 4),

/* XString_class.c */
	CALLMETHOD_DEF(debug_XString_class, 0),
	CALLMETHOD_DEF(init_DNAlkups, 2),
	CALLMETHOD_DEF(init_RNAlkups, 2),
	CALLMETHOD_DEF(new_RawPtr_from_XString, 4),

/* XStringSet_class.c */
	CALLMETHOD_DEF(debug_XStringSet_class, 0),
	CALLMETHOD_DEF(XStringSet_as_STRSXP, 2),
	CALLMETHOD_DEF(XStringSet_order, 1),

/* char_frequency.c */
	CALLMETHOD_DEF(XString_char_frequency, 3),
	CALLMETHOD_DEF(XStringSet_char_frequency, 4),
	CALLMETHOD_DEF(oligonucleotide_frequency, 4),

/* char_translate.c */
	CALLMETHOD_DEF(XStringSet_char_translate, 3),

/* replace_letter_at.c */
	CALLMETHOD_DEF(XString_replace_letter_at, 6),
	CALLMETHOD_DEF(XString_inplace_replace_letter_at, 4),

/* inject_code.c */
	CALLMETHOD_DEF(inject_code, 4),

/* match_utils.c */
	CALLMETHOD_DEF(debug_match_utils, 0),
	CALLMETHOD_DEF(nmismatch_at, 5),
	CALLMETHOD_DEF(is_matching, 5),
	CALLMETHOD_DEF(nmatch_PairwiseAlignment, 4),

/* match_reporting.c */
	CALLMETHOD_DEF(debug_match_reporting, 0),

/* match_pattern_boyermoore.c */
	CALLMETHOD_DEF(debug_match_pattern_boyermoore, 0),

/* match_pattern_shiftor.c */
	CALLMETHOD_DEF(debug_match_pattern_shiftor, 0),
	CALLMETHOD_DEF(bits_per_long, 0),

/* match_pattern.c */
	CALLMETHOD_DEF(debug_match_pattern, 0),
	CALLMETHOD_DEF(XString_match_pattern, 6),
	CALLMETHOD_DEF(XStringViews_match_pattern, 8),
	CALLMETHOD_DEF(XStringSet_vmatch_pattern, 6),

/* match_BOC.c */
	CALLMETHOD_DEF(debug_match_BOC, 0),
	CALLMETHOD_DEF(match_BOC_preprocess, 12),
	CALLMETHOD_DEF(match_BOC_exact, 16),

/* match_BOC2.c */
	CALLMETHOD_DEF(debug_match_BOC2, 0),
	CALLMETHOD_DEF(match_BOC2_preprocess, 9),
	CALLMETHOD_DEF(match_BOC2_exact, 13),

/* match_PWM.c */
	CALLMETHOD_DEF(PWM_score, 3),
	CALLMETHOD_DEF(match_PWM, 4),

/* find_palindromes.c */
	CALLMETHOD_DEF(debug_find_palindromes, 0),
	CALLMETHOD_DEF(find_palindromes, 6),

/* MIndex_utils.c */
	CALLMETHOD_DEF(debug_MIndex_utils, 0),
	CALLMETHOD_DEF(ByPos_MIndex_endIndex, 3),
	CALLMETHOD_DEF(ByName_MIndex_endIndex, 4),
	CALLMETHOD_DEF(ByPos_MIndex_coverage, 4),
	CALLMETHOD_DEF(ByName_MIndex_coverage, 4),
	CALLMETHOD_DEF(ByPos_MIndex_combine, 1),

/* match_pdict_Twobit.c */
	CALLMETHOD_DEF(debug_match_pdict_Twobit, 0),
	CALLMETHOD_DEF(build_Twobit, 3),

/* match_pdict_ACtree.c */
	CALLMETHOD_DEF(debug_match_pdict_ACtree, 0),
	CALLMETHOD_DEF(free_actree_nodes_buf, 0),
	CALLMETHOD_DEF(build_ACtree, 3),

/* match_pdict.c */
	CALLMETHOD_DEF(debug_match_pdict, 0),
	CALLMETHOD_DEF(XString_match_pdict, 12),
	CALLMETHOD_DEF(XStringViews_match_pdict, 14),

/* pmatchPattern.c */
	CALLMETHOD_DEF(lcprefix, 6),
	CALLMETHOD_DEF(lcsuffix, 6),

/* align_needwunsQS.c */
	CALLMETHOD_DEF(align_needwunsQS, 7),

/* align_pairwiseAlignment.c */
	CALLMETHOD_DEF(XStringSet_align_pairwiseAlignment, 15),
	CALLMETHOD_DEF(XStringSet_align_distance, 13),

/* align_utils.c */
	CALLMETHOD_DEF(AlignedXStringSet_nchar, 1),
	CALLMETHOD_DEF(AlignedXStringSet_align_aligned, 2),
	CALLMETHOD_DEF(align_compareStrings, 6),

	{NULL, NULL, 0}
};

void R_init_Biostrings(DllInfo *info)
{
	/* Lots of code around assumes that sizeof(Rbyte) == sizeof(char) */
	if (sizeof(Rbyte) != sizeof(char))
		error("sizeof(Rbyte) != sizeof(char)");
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);

/* RoSeq_utils.c */
	REGISTER_CCALLABLE(_new_STRSXP_from_RoSeqs);
	REGISTER_CCALLABLE(_new_RoSeqs_from_CharAEAE);
	REGISTER_CCALLABLE(_new_IRanges_from_RoSeqs);

/* XString_class.c */
	REGISTER_CCALLABLE(_DNAencode);
	REGISTER_CCALLABLE(_DNAdecode);
	REGISTER_CCALLABLE(_RNAencode);
	REGISTER_CCALLABLE(_RNAdecode);
	REGISTER_CCALLABLE(_get_XString_asRoSeq);

/* XStringSet_class.c */
	REGISTER_CCALLABLE(_get_XStringSet_baseClass);
	REGISTER_CCALLABLE(_get_XStringSet_length);
	REGISTER_CCALLABLE(_new_CachedXStringSet);
	REGISTER_CCALLABLE(_get_CachedXStringSet_elt_asRoSeq);
	REGISTER_CCALLABLE(_get_XStringSet_elt_asRoSeq);
	REGISTER_CCALLABLE(_new_XStringSet_from_RoSeqs);
	REGISTER_CCALLABLE(_set_XStringSet_names);
	REGISTER_CCALLABLE(_alloc_XStringSet);
	REGISTER_CCALLABLE(_write_RoSeq_to_CachedXStringSet_elt);
	REGISTER_CCALLABLE(_write_RoSeq_to_XStringSet_elt);

/* match_reporting.c */
	REGISTER_CCALLABLE(_init_match_reporting);
	REGISTER_CCALLABLE(_report_match);
	REGISTER_CCALLABLE(_reported_matches_asSEXP);

	return;
}

