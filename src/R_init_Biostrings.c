#include "Biostrings.h"

static const R_CallMethodDef callMethods[] = {

/* copy_seq.c */
	{"debug_copy_seq", (DL_FUNC) &debug_copy_seq, 0},

/* utils.c */
	{"debug_utils", (DL_FUNC) &debug_utils, 0},
	{"Biostrings_length_vectors_in_list", (DL_FUNC) &Biostrings_length_vectors_in_list, 1},

/* bufutils.c */
	{"debug_bufutils", (DL_FUNC) &debug_bufutils, 0},

/* IRanges_class.c */
	{"debug_IRanges_class", (DL_FUNC) &debug_IRanges_class, 0},

/* IRanges_utils.c */
	{"debug_IRanges_utils", (DL_FUNC) &debug_IRanges_utils, 0},

	{"narrow_IRanges", (DL_FUNC) &narrow_IRanges, 4},
	{"int_to_adjacent_ranges", (DL_FUNC) &int_to_adjacent_ranges, 1},
	{"reduce_IRanges", (DL_FUNC) &reduce_IRanges, 2},
	{"IRanges_coverage", (DL_FUNC) &IRanges_coverage, 3},

/* SparseList_utils.c */
	{"debug_SparseList_utils", (DL_FUNC) &debug_SparseList_utils, 0},

/* XRaw_class.c */
	{"debug_XRaw_class", (DL_FUNC) &debug_XRaw_class, 0},

	{"Biostrings_sexp_address", (DL_FUNC) &Biostrings_sexp_address, 1},
	{"Biostrings_safe_explode", (DL_FUNC) &Biostrings_safe_explode, 1},
	{"Biostrings_xp_show", (DL_FUNC) &Biostrings_xp_show, 1},
	{"Biostrings_xp_new", (DL_FUNC) &Biostrings_xp_new, 0},
	{"Biostrings_XRaw_alloc", (DL_FUNC) &Biostrings_XRaw_alloc, 2},
	{"Biostrings_XRaw_get_show_string", (DL_FUNC) &Biostrings_XRaw_get_show_string, 1},
	{"Biostrings_XRaw_length", (DL_FUNC) &Biostrings_XRaw_length, 1},

/* XRaw_utils.c */
	{"debug_XRaw_utils", (DL_FUNC) &debug_XRaw_utils, 0},

	{"Biostrings_XRaw_memcmp", (DL_FUNC) &Biostrings_XRaw_memcmp, 5},

	{"Biostrings_XRaw_copy_from_i1i2", (DL_FUNC) &Biostrings_XRaw_copy_from_i1i2, 4},
	{"Biostrings_XRaw_copy_from_subset", (DL_FUNC) &Biostrings_XRaw_copy_from_subset, 3},

	{"Biostrings_XRaw_read_chars_from_i1i2", (DL_FUNC) &Biostrings_XRaw_read_chars_from_i1i2, 3},
	{"Biostrings_XRaw_read_chars_from_subset", (DL_FUNC) &Biostrings_XRaw_read_chars_from_subset, 2},
	{"Biostrings_XRaw_write_chars_to_i1i2", (DL_FUNC) &Biostrings_XRaw_write_chars_to_i1i2, 4},
	{"Biostrings_XRaw_write_chars_to_subset", (DL_FUNC) &Biostrings_XRaw_write_chars_to_subset, 3},

	{"XRaw_read_ints_from_i1i2", (DL_FUNC) &XRaw_read_ints_from_i1i2, 3},
	{"XRaw_read_ints_from_subset", (DL_FUNC) &XRaw_read_ints_from_subset, 2},
	{"XRaw_write_ints_to_i1i2", (DL_FUNC) &XRaw_write_ints_to_i1i2, 4},
	{"XRaw_write_ints_to_subset", (DL_FUNC) &XRaw_write_ints_to_subset, 3},

	{"XRaw_read_enc_chars_from_i1i2", (DL_FUNC) &XRaw_read_enc_chars_from_i1i2, 4},
	{"XRaw_read_enc_chars_from_subset", (DL_FUNC) &XRaw_read_enc_chars_from_subset, 3},
	{"XRaw_write_enc_chars_to_i1i2", (DL_FUNC) &XRaw_write_enc_chars_to_i1i2, 5},
	{"XRaw_write_enc_chars_to_subset", (DL_FUNC) &XRaw_write_enc_chars_to_subset, 4},

	{"XRaw_read_complexes_from_i1i2", (DL_FUNC) &XRaw_read_complexes_from_i1i2, 4},
	{"XRaw_read_complexes_from_subset", (DL_FUNC) &XRaw_read_complexes_from_subset, 3},

	{"XRaw_loadFASTA", (DL_FUNC) &XRaw_loadFASTA, 4},

/* XInteger.c */
	{"debug_XInteger", (DL_FUNC) &debug_XInteger, 0},

	{"XInteger_alloc", (DL_FUNC) &XInteger_alloc, 2},
	{"XInteger_get_show_string", (DL_FUNC) &XInteger_get_show_string, 1},
	{"XInteger_length", (DL_FUNC) &XInteger_length, 1},
	{"XInteger_memcmp", (DL_FUNC) &XInteger_memcmp, 5},

	{"XInteger_read_ints_from_i1i2", (DL_FUNC) &XInteger_read_ints_from_i1i2, 3},
	{"XInteger_read_ints_from_subset", (DL_FUNC) &XInteger_read_ints_from_subset, 2},
	{"XInteger_write_ints_to_i1i2", (DL_FUNC) &XInteger_write_ints_to_i1i2, 4},
	{"XInteger_write_ints_to_subset", (DL_FUNC) &XInteger_write_ints_to_subset, 3},

/* XNumeric.c */
	{"debug_XNumeric", (DL_FUNC) &debug_XNumeric, 0},

	{"XNumeric_alloc", (DL_FUNC) &XNumeric_alloc, 2},
	{"XNumeric_get_show_string", (DL_FUNC) &XNumeric_get_show_string, 1},
	{"XNumeric_length", (DL_FUNC) &XNumeric_length, 1},
	{"XNumeric_memcmp", (DL_FUNC) &XNumeric_memcmp, 5},

	{"XNumeric_read_nums_from_i1i2", (DL_FUNC) &XNumeric_read_nums_from_i1i2, 3},
	{"XNumeric_read_nums_from_subset", (DL_FUNC) &XNumeric_read_nums_from_subset, 2},
	{"XNumeric_write_nums_to_i1i2", (DL_FUNC) &XNumeric_write_nums_to_i1i2, 4},
	{"XNumeric_write_nums_to_subset", (DL_FUNC) &XNumeric_write_nums_to_subset, 3},

/* XString_class.c */
	{"debug_XString_class", (DL_FUNC) &debug_XString_class, 0},
	{"init_DNAlkups", (DL_FUNC) &init_DNAlkups, 2},
	{"init_RNAlkups", (DL_FUNC) &init_RNAlkups, 2},

/* XStringSet_class.c */
	{"debug_XStringSet_class", (DL_FUNC) &debug_XStringSet_class, 0},
	{"XStringSet_as_STRSXP", (DL_FUNC) &XStringSet_as_STRSXP, 2},
	{"XStrings_to_nchars", (DL_FUNC) &XStrings_to_nchars, 1},

/* seqs_to_seqs.c */
	{"debug_seqs_to_seqs", (DL_FUNC) &debug_seqs_to_seqs, 0},

	{"copy_subXRaw", (DL_FUNC) &copy_subXRaw, 4},
	{"new_XRaw_from_STRSXP", (DL_FUNC) &new_XRaw_from_STRSXP, 5},
	{"new_XRaw_from_XString", (DL_FUNC) &new_XRaw_from_XString, 4},
	{"new_XStringList_from_XRaw", (DL_FUNC) &new_XStringList_from_XRaw, 4},
	{"narrow_XStringList", (DL_FUNC) &narrow_XStringList, 4},
        
/* char_frequency.c */
	{"XString_char_frequency", (DL_FUNC) &XString_char_frequency, 3},
	{"XStringSet_char_frequency", (DL_FUNC) &XStringSet_char_frequency, 4},
	{"oligonucleotide_frequency", (DL_FUNC) &oligonucleotide_frequency, 4},

/* char_translate.c */
        {"XRaw_translate_copy_from_i1i2", (DL_FUNC) &XRaw_translate_copy_from_i1i2, 5},
        {"XRaw_translate_copy_from_subset", (DL_FUNC) &XRaw_translate_copy_from_subset, 4},
        {"XRaw_reverse_copy_from_i1i2", (DL_FUNC) &XRaw_reverse_copy_from_i1i2, 4},
        {"XRaw_reverse_translate_copy_from_i1i2", (DL_FUNC) &XRaw_reverse_translate_copy_from_i1i2, 5},
        {"XStringSet_char_translate", (DL_FUNC) &XStringSet_char_translate, 3},

/* replace_locs.c */
	{"XString_replace_locs_bySTRSXP", (DL_FUNC) &XString_replace_locs_bySTRSXP, 6},
	{"XString_inplace_replace_locs_bySTRSXP", (DL_FUNC) &XString_inplace_replace_locs_bySTRSXP, 4},

/* inject_code.c */
	{"inject_code", (DL_FUNC) &inject_code, 4},

/* match_utils.c */
	{"debug_match_utils", (DL_FUNC) &debug_match_utils, 0},
	{"nmismatch_at", (DL_FUNC) &nmismatch_at, 5},
	{"is_matching", (DL_FUNC) &is_matching, 5},

/* match_reporting.c */
	{"debug_match_reporting", (DL_FUNC) &debug_match_reporting, 0},

/* match_pattern_boyermoore.c */
	{"debug_match_pattern_boyermoore", (DL_FUNC) &debug_match_pattern_boyermoore, 0},

/* match_pattern_shiftor.c */
	{"debug_match_pattern_shiftor", (DL_FUNC) &debug_match_pattern_shiftor, 0},
	{"bits_per_long", (DL_FUNC) &bits_per_long, 0},

/* match_pattern.c */
	{"debug_match_pattern", (DL_FUNC) &debug_match_pattern, 0},
	{"XString_match_pattern", (DL_FUNC) &XString_match_pattern, 6},
	{"XStringViews_match_pattern", (DL_FUNC) &XStringViews_match_pattern, 8},
	{"XStringSet_vmatch_pattern", (DL_FUNC) &XStringSet_vmatch_pattern, 6},

/* match_BOC.c */
	{"debug_match_BOC", (DL_FUNC) &debug_match_BOC, 0},
	{"match_BOC_preprocess", (DL_FUNC) &match_BOC_preprocess, 12},
	{"match_BOC_exact", (DL_FUNC) &match_BOC_exact, 16},

/* match_BOC2.c */
	{"debug_match_BOC2", (DL_FUNC) &debug_match_BOC2, 0},
	{"match_BOC2_preprocess", (DL_FUNC) &match_BOC2_preprocess, 9},
	{"match_BOC2_exact", (DL_FUNC) &match_BOC2_exact, 13},

/* match_PWM.c */
	{"PWM_score", (DL_FUNC) &PWM_score, 3},
	{"match_PWM", (DL_FUNC) &match_PWM, 4},

/* find_palindromes.c */
	{"debug_find_palindromes", (DL_FUNC) &debug_find_palindromes, 0},
	{"find_palindromes", (DL_FUNC) &find_palindromes, 6},

/* PDict_utils.c */
	{"debug_PDict_utils", (DL_FUNC) &debug_PDict_utils, 0},

/* MIndex_utils.c */
	{"debug_MIndex_utils", (DL_FUNC) &debug_MIndex_utils, 0},
	{"shiftListOfInts", (DL_FUNC) &shiftListOfInts, 2},
	{"extract_endIndex", (DL_FUNC) &extract_endIndex, 4},
	{"ByPos_MIndex_coverage", (DL_FUNC) &ByPos_MIndex_coverage, 4},
	{"ByName_MIndex_coverage", (DL_FUNC) &ByName_MIndex_coverage, 4},

/* match_pdict_Twobit.c */
	{"debug_match_pdict_Twobit", (DL_FUNC) &debug_match_pdict_Twobit, 0},
	{"build_Twobit", (DL_FUNC) &build_Twobit, 2},

/* match_pdict_ACtree.c */
	{"debug_match_pdict_ACtree", (DL_FUNC) &debug_match_pdict_ACtree, 0},
	{"free_actree_nodes_buf", (DL_FUNC) &free_actree_nodes_buf, 0},
	{"build_ACtree", (DL_FUNC) &build_ACtree, 2},

/* match_pdict.c */
	{"debug_match_pdict", (DL_FUNC) &debug_match_pdict, 0},
	{"XString_match_pdict", (DL_FUNC) &XString_match_pdict, 12},
	{"XStringViews_match_pdict", (DL_FUNC) &XStringViews_match_pdict, 14},

/* pmatchPattern.c */
	{"lcprefix", (DL_FUNC) &lcprefix, 6},
	{"lcsuffix", (DL_FUNC) &lcsuffix, 6},

/* align_needwunsQS.c */
	{"align_needwunsQS", (DL_FUNC) &align_needwunsQS, 7},

/* align_pairwiseAlignment.c */
	{"XStringSet_align_pairwiseAlignment", (DL_FUNC) &XStringSet_align_pairwiseAlignment, 17},

/* align_utils.c */
	{"AlignedXStringSet_nchar", (DL_FUNC) &AlignedXStringSet_nchar, 1},
	{"AlignedXStringSet_align_aligned", (DL_FUNC) &AlignedXStringSet_align_aligned, 2},
	{"align_compareStrings", (DL_FUNC) &align_compareStrings, 6},

	{NULL, NULL, 0}
};

void R_init_Biostrings(DllInfo *info)
{
	/* Lots of code around assumes that sizeof(Rbyte) == sizeof(char) */
	if (sizeof(Rbyte) != sizeof(char))
		error("sizeof(Rbyte) != sizeof(char)");
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	R_RegisterCCallable("Biostrings", "_new_CharBuf_from_string", (DL_FUNC) &_new_CharBuf_from_string);
	R_RegisterCCallable("Biostrings", "_new_CharBBuf", (DL_FUNC) &_new_CharBBuf);
	R_RegisterCCallable("Biostrings", "_append_string_to_CharBBuf", (DL_FUNC) &_append_string_to_CharBBuf);
	R_RegisterCCallable("Biostrings", "_DNAencode", (DL_FUNC) &_DNAencode);
	R_RegisterCCallable("Biostrings", "_DNAdecode", (DL_FUNC) &_DNAdecode);
	R_RegisterCCallable("Biostrings", "_RNAencode", (DL_FUNC) &_RNAencode);
	R_RegisterCCallable("Biostrings", "_RNAdecode", (DL_FUNC) &_RNAdecode);
	R_RegisterCCallable("Biostrings", "_get_XString_asRoSeq", (DL_FUNC) &_get_XString_asRoSeq);
	R_RegisterCCallable("Biostrings", "_get_XStringSet_baseClass", (DL_FUNC) &_get_XStringSet_baseClass);
	R_RegisterCCallable("Biostrings", "_get_XStringSet_length", (DL_FUNC) &_get_XStringSet_length);
	R_RegisterCCallable("Biostrings", "_new_CachedXStringSet", (DL_FUNC) &_new_CachedXStringSet);
	R_RegisterCCallable("Biostrings", "_get_CachedXStringSet_elt_asRoSeq", (DL_FUNC) &_get_CachedXStringSet_elt_asRoSeq);
	R_RegisterCCallable("Biostrings", "_get_XStringSet_elt_asRoSeq", (DL_FUNC) &_get_XStringSet_elt_asRoSeq);
	R_RegisterCCallable("Biostrings", "_new_XStringSet_from_RoSeqs", (DL_FUNC) &_new_XStringSet_from_RoSeqs);
	R_RegisterCCallable("Biostrings", "_set_XStringSet_names", (DL_FUNC) &_set_XStringSet_names);
	R_RegisterCCallable("Biostrings", "_alloc_XStringSet", (DL_FUNC) &_alloc_XStringSet);
	R_RegisterCCallable("Biostrings", "_write_RoSeq_to_CachedXStringSet_elt", (DL_FUNC) &_write_RoSeq_to_CachedXStringSet_elt);
	R_RegisterCCallable("Biostrings", "_write_RoSeq_to_XStringSet_elt", (DL_FUNC) &_write_RoSeq_to_XStringSet_elt);
	R_RegisterCCallable("Biostrings", "_new_RoSeqs_from_BBuf", (DL_FUNC) &_new_RoSeqs_from_BBuf);
	R_RegisterCCallable("Biostrings", "_new_STRSXP_from_RoSeqs", (DL_FUNC) &_new_STRSXP_from_RoSeqs);
	R_RegisterCCallable("Biostrings", "_init_match_reporting", (DL_FUNC) &_init_match_reporting);
	R_RegisterCCallable("Biostrings", "_report_match", (DL_FUNC) &_report_match);
	R_RegisterCCallable("Biostrings", "_reported_matches_asSEXP", (DL_FUNC) &_reported_matches_asSEXP);
}

