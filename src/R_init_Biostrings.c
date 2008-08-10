#include "Biostrings.h"

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

#define REGISTER_CCALLABLE(fun) \
	R_RegisterCCallable("Biostrings", #fun, (DL_FUNC) &fun)

static const R_CallMethodDef callMethods[] = {

/* copy_seq.c */
	CALLMETHOD_DEF(debug_copy_seq, 0),

/* utils.c */
	CALLMETHOD_DEF(debug_utils, 0),
	CALLMETHOD_DEF(Biostrings_length_vectors_in_list, 1),

/* RoSeq_utils.c */
	CALLMETHOD_DEF(debug_RoSeq_utils, 0),

/* SparseList_utils.c */
	CALLMETHOD_DEF(debug_SparseList_utils, 0),

/* XRaw_class.c */
	CALLMETHOD_DEF(debug_XRaw_class, 0),

	CALLMETHOD_DEF(Biostrings_sexp_address, 1),
	CALLMETHOD_DEF(Biostrings_safe_explode, 1),
	CALLMETHOD_DEF(Biostrings_xp_show, 1),
	CALLMETHOD_DEF(Biostrings_xp_new, 0),
	CALLMETHOD_DEF(Biostrings_XRaw_alloc, 2),
	CALLMETHOD_DEF(Biostrings_XRaw_get_show_string, 1),
	CALLMETHOD_DEF(Biostrings_XRaw_length, 1),

/* XRaw_utils.c */
	CALLMETHOD_DEF(debug_XRaw_utils, 0),

	CALLMETHOD_DEF(Biostrings_XRaw_memcmp, 5),

	CALLMETHOD_DEF(Biostrings_XRaw_copy_from_i1i2, 4),
	CALLMETHOD_DEF(Biostrings_XRaw_copy_from_subset, 3),

	CALLMETHOD_DEF(Biostrings_XRaw_read_chars_from_i1i2, 3),
	CALLMETHOD_DEF(Biostrings_XRaw_read_chars_from_subset, 2),
	CALLMETHOD_DEF(Biostrings_XRaw_write_chars_to_i1i2, 4),
	CALLMETHOD_DEF(Biostrings_XRaw_write_chars_to_subset, 3),

	CALLMETHOD_DEF(XRaw_read_ints_from_i1i2, 3),
	CALLMETHOD_DEF(XRaw_read_ints_from_subset, 2),
	CALLMETHOD_DEF(XRaw_write_ints_to_i1i2, 4),
	CALLMETHOD_DEF(XRaw_write_ints_to_subset, 3),

	CALLMETHOD_DEF(XRaw_read_enc_chars_from_i1i2, 4),
	CALLMETHOD_DEF(XRaw_read_enc_chars_from_subset, 3),
	CALLMETHOD_DEF(XRaw_write_enc_chars_to_i1i2, 5),
	CALLMETHOD_DEF(XRaw_write_enc_chars_to_subset, 4),

	CALLMETHOD_DEF(XRaw_read_complexes_from_i1i2, 4),
	CALLMETHOD_DEF(XRaw_read_complexes_from_subset, 3),

	CALLMETHOD_DEF(XRaw_loadFASTA, 4),

/* XInteger.c */
	CALLMETHOD_DEF(debug_XInteger, 0),

	CALLMETHOD_DEF(XInteger_alloc, 2),
	CALLMETHOD_DEF(XInteger_get_show_string, 1),
	CALLMETHOD_DEF(XInteger_length, 1),
	CALLMETHOD_DEF(XInteger_memcmp, 5),

	CALLMETHOD_DEF(XInteger_read_ints_from_i1i2, 3),
	CALLMETHOD_DEF(XInteger_read_ints_from_subset, 2),
	CALLMETHOD_DEF(XInteger_write_ints_to_i1i2, 4),
	CALLMETHOD_DEF(XInteger_write_ints_to_subset, 3),

/* XNumeric.c */
	CALLMETHOD_DEF(debug_XNumeric, 0),

	CALLMETHOD_DEF(XNumeric_alloc, 2),
	CALLMETHOD_DEF(XNumeric_get_show_string, 1),
	CALLMETHOD_DEF(XNumeric_length, 1),
	CALLMETHOD_DEF(XNumeric_memcmp, 5),

	CALLMETHOD_DEF(XNumeric_read_nums_from_i1i2, 3),
	CALLMETHOD_DEF(XNumeric_read_nums_from_subset, 2),
	CALLMETHOD_DEF(XNumeric_write_nums_to_i1i2, 4),
	CALLMETHOD_DEF(XNumeric_write_nums_to_subset, 3),

/* XString_class.c */
	CALLMETHOD_DEF(debug_XString_class, 0),
	CALLMETHOD_DEF(init_DNAlkups, 2),
	CALLMETHOD_DEF(init_RNAlkups, 2),

/* XStringSet_class.c */
	CALLMETHOD_DEF(debug_XStringSet_class, 0),
	CALLMETHOD_DEF(XStringSet_as_STRSXP, 2),

/* seqs_to_seqs.c */
	CALLMETHOD_DEF(debug_seqs_to_seqs,0),

	CALLMETHOD_DEF(copy_subXRaw, 4),
	CALLMETHOD_DEF(new_XRaw_from_STRSXP, 5),
	CALLMETHOD_DEF(new_XRaw_from_XString, 4),

/* char_frequency.c */
	CALLMETHOD_DEF(XString_char_frequency, 3),
	CALLMETHOD_DEF(XStringSet_char_frequency, 4),
	CALLMETHOD_DEF(oligonucleotide_frequency, 4),

/* char_translate.c */
	CALLMETHOD_DEF(XRaw_translate_copy_from_i1i2, 5),
	CALLMETHOD_DEF(XRaw_translate_copy_from_subset, 4),
	CALLMETHOD_DEF(XRaw_reverse_copy_from_i1i2, 4),
	CALLMETHOD_DEF(XRaw_reverse_translate_copy_from_i1i2, 5),
	CALLMETHOD_DEF(XStringSet_char_translate, 3),

/* replace_locs.c */
	CALLMETHOD_DEF(XString_replace_locs_bySTRSXP, 6),
	CALLMETHOD_DEF(XString_inplace_replace_locs_bySTRSXP, 4),

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

/* PDict_utils.c */
	CALLMETHOD_DEF(debug_PDict_utils, 0),
	CALLMETHOD_DEF(Dups_diff, 2),

/* MIndex_utils.c */
	CALLMETHOD_DEF(debug_MIndex_utils, 0),
	CALLMETHOD_DEF(shiftListOfInts, 2),
	CALLMETHOD_DEF(extract_endIndex, 4),
	CALLMETHOD_DEF(ByPos_MIndex_coverage, 4),
	CALLMETHOD_DEF(ByName_MIndex_coverage, 4),
	CALLMETHOD_DEF(ByPos_MIndex_combine, 1),

/* match_pdict_Twobit.c */
	CALLMETHOD_DEF(debug_match_pdict_Twobit, 0),
	CALLMETHOD_DEF(build_Twobit, 2),

/* match_pdict_ACtree.c */
	CALLMETHOD_DEF(debug_match_pdict_ACtree, 0),
	CALLMETHOD_DEF(free_actree_nodes_buf, 0),
	CALLMETHOD_DEF(build_ACtree, 2),

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
	CALLMETHOD_DEF(XStringSet_align_pairwiseAlignment, 18),
	CALLMETHOD_DEF(XStringSet_align_distance, 14),

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
	REGISTER_CCALLABLE(_new_IRanges_from_RoSeqs);

/* XRaw_class.c */
	REGISTER_CCALLABLE(_new_STRSXP_from_RoSeqs);

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

/* seqs_to_seqs.c */
	REGISTER_CCALLABLE(_new_RoSeqs_from_CharAEAE);

/* match_reporting.c */
	REGISTER_CCALLABLE(_init_match_reporting);
	REGISTER_CCALLABLE(_report_match);
	REGISTER_CCALLABLE(_reported_matches_asSEXP);

	return;
}

