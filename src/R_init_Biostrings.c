#include "Biostrings.h"

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

#define REGISTER_CCALLABLE(fun) \
	R_RegisterCCallable("Biostrings", #fun, (DL_FUNC) &fun)

static R_NativePrimitiveArgType gtestsim_t[11] = {
	INTSXP,  /* int *nrow */
	INTSXP,  /* int *ncol */
	INTSXP,	 /* int *nrowt */
	INTSXP,  /* int *ncolt */
	INTSXP,  /* int *n */
	INTSXP,  /* int *b */
	REALSXP, /* double *expected */
	INTSXP,  /* int *observed */
	REALSXP, /* double *fact */
	INTSXP,  /* int *jwork */
	REALSXP  /* double *results */
};

static const R_CMethodDef cMethods[] = {
	{"gtestsim", (DL_FUNC) &gtestsim, 11, gtestsim_t},
	{NULL, NULL, 0}
};

static const R_CallMethodDef callMethods[] = {

/* utils.c */
	CALLMETHOD_DEF(debug_utils, 0),

/* io_utils.c */
	CALLMETHOD_DEF(new_input_ExternalFilePtr, 1),
	CALLMETHOD_DEF(new_output_ExternalFilePtr, 2),
	CALLMETHOD_DEF(ExternalFilePtr_close, 1),

/* RoSeqs_utils.c */
	CALLMETHOD_DEF(debug_RoSeqs_utils, 0),

/* XString_class.c */
	CALLMETHOD_DEF(debug_XString_class, 0),
	CALLMETHOD_DEF(init_DNAlkups, 2),
	CALLMETHOD_DEF(init_RNAlkups, 2),
	CALLMETHOD_DEF(new_XString_from_CHARACTER, 5),

/* XStringSet_class.c */
	CALLMETHOD_DEF(debug_XStringSet_class, 0),
	CALLMETHOD_DEF(new_XStringSet_from_CHARACTER, 6),
	CALLMETHOD_DEF(XStringSet_unlist, 1),
	CALLMETHOD_DEF(XStringSet_as_STRSXP, 2),
	CALLMETHOD_DEF(XStringSet_is_unsorted, 2),
	CALLMETHOD_DEF(XStringSet_order, 1),
	CALLMETHOD_DEF(XStringSet_duplicated, 1),
	CALLMETHOD_DEF(XStringSet_match, 3),

/* xscat.c */
	CALLMETHOD_DEF(XString_xscat, 1),
	CALLMETHOD_DEF(XStringSet_xscat, 1),

/* XStringSet_io.c */
	CALLMETHOD_DEF(debug_XStringSet_io, 0),
	CALLMETHOD_DEF(fasta_info, 4),
	CALLMETHOD_DEF(read_fasta_in_XStringSet, 6),
	CALLMETHOD_DEF(write_XStringSet_to_fasta, 4),
	CALLMETHOD_DEF(fastq_geometry, 3),
	CALLMETHOD_DEF(read_fastq_in_XStringSet, 6),

/* letter_frequency.c */
	CALLMETHOD_DEF(XString_letter_frequency, 3),
	CALLMETHOD_DEF(XStringSet_letter_frequency, 4),
	CALLMETHOD_DEF(XString_letterFrequencyInSlidingView, 4),
	CALLMETHOD_DEF(XString_oligo_frequency, 7),
	CALLMETHOD_DEF(XStringSet_oligo_frequency, 8),
	CALLMETHOD_DEF(XStringSet_nucleotide_frequency_at, 7),
	CALLMETHOD_DEF(XStringSet_consensus_matrix, 5),

/* extract_transcripts.c */
	CALLMETHOD_DEF(transcript_widths, 2),
	CALLMETHOD_DEF(extract_transcripts, 6),
	CALLMETHOD_DEF(tlocs2rlocs, 5),

/* translate.c */
	CALLMETHOD_DEF(DNAStringSet_translate, 4),

/* replace_letter_at.c */
	CALLMETHOD_DEF(XString_replace_letter_at, 6),
	CALLMETHOD_DEF(XString_inplace_replace_letter_at, 4),

/* inject_code.c */
	CALLMETHOD_DEF(XString_inject_code, 4),

/* SparseList_utils.c */
	CALLMETHOD_DEF(debug_SparseList_utils, 0),

/* match_reporting.c */
	CALLMETHOD_DEF(debug_match_reporting, 0),

/* MIndex_class.c */
	CALLMETHOD_DEF(debug_MIndex_class, 0),
	CALLMETHOD_DEF(ByPos_MIndex_endIndex, 3),
	CALLMETHOD_DEF(SparseMIndex_endIndex, 4),
	CALLMETHOD_DEF(ByPos_MIndex_combine, 1),

/* lowlevel_matching.c */
	CALLMETHOD_DEF(debug_lowlevel_matching, 0),
	CALLMETHOD_DEF(XString_match_pattern_at, 10),
	CALLMETHOD_DEF(XStringSet_vmatch_pattern_at, 10),
	CALLMETHOD_DEF(XStringSet_dist_hamming, 1),

/* match_pattern_boyermoore.c */
	CALLMETHOD_DEF(debug_match_pattern_boyermoore, 0),

/* match_pattern_shiftor.c */
	CALLMETHOD_DEF(debug_match_pattern_shiftor, 0),
	CALLMETHOD_DEF(bits_per_long, 0),

/* match_pattern_indels.c */
	CALLMETHOD_DEF(debug_match_pattern_indels, 0),

/* match_pattern.c */
	CALLMETHOD_DEF(debug_match_pattern, 0),
	CALLMETHOD_DEF(XString_match_pattern, 8),
	CALLMETHOD_DEF(XStringViews_match_pattern, 10),
	CALLMETHOD_DEF(XStringSet_vmatch_pattern, 8),

/* match_BOC.c */
	CALLMETHOD_DEF(debug_match_BOC, 0),
	CALLMETHOD_DEF(match_BOC_preprocess, 12),
	CALLMETHOD_DEF(match_BOC_exact, 16),

/* match_BOC2.c */
	CALLMETHOD_DEF(debug_match_BOC2, 0),
	CALLMETHOD_DEF(match_BOC2_preprocess, 9),
	CALLMETHOD_DEF(match_BOC2_exact, 13),

/* match_PWM.c */
	CALLMETHOD_DEF(PWM_score_starting_at, 4),
	CALLMETHOD_DEF(XString_match_PWM, 5),
	CALLMETHOD_DEF(XStringViews_match_PWM, 7),

/* match_WCP.c */
	CALLMETHOD_DEF(WCP_score_starting_at, 3),
	CALLMETHOD_DEF(XString_match_WCP, 4),
	CALLMETHOD_DEF(XStringViews_match_PWM, 6),

/* find_palindromes.c */
	CALLMETHOD_DEF(debug_find_palindromes, 0),
	CALLMETHOD_DEF(find_palindromes, 6),

/* BitMatrix.c */
	CALLMETHOD_DEF(debug_BitMatrix, 0),

/* PreprocessedTB_class.c */
	CALLMETHOD_DEF(debug_PreprocessedTB_class, 0),

/* match_pdict_utils.c */
	CALLMETHOD_DEF(debug_match_pdict_utils, 0),

/* match_pdict_Twobit.c */
	CALLMETHOD_DEF(debug_match_pdict_Twobit, 0),
	CALLMETHOD_DEF(build_Twobit, 3),

/* BAB_class.c */
	CALLMETHOD_DEF(debug_BAB_class, 0),
	CALLMETHOD_DEF(IntegerBAB_new, 1),

/* match_pdict_ACtree2.c */
	CALLMETHOD_DEF(debug_match_pdict_ACtree2, 0),
	CALLMETHOD_DEF(ACtree2_nodebuf_max_nblock, 0),
	CALLMETHOD_DEF(ACtree2_nodeextbuf_max_nblock, 0),
	CALLMETHOD_DEF(ACtree2_nnodes, 1),
	CALLMETHOD_DEF(ACtree2_print_nodes, 1),
	CALLMETHOD_DEF(ACtree2_summary, 1),
	CALLMETHOD_DEF(ACtree2_build, 5),
	CALLMETHOD_DEF(ACtree2_has_all_flinks, 1),
	CALLMETHOD_DEF(ACtree2_compute_all_flinks, 1),

/* match_pdict.c */
	CALLMETHOD_DEF(debug_match_pdict, 0),
	CALLMETHOD_DEF(match_PDict3Parts_XString, 9),
	CALLMETHOD_DEF(match_XStringSet_XString, 8),
	CALLMETHOD_DEF(match_PDict3Parts_XStringViews, 11),
	CALLMETHOD_DEF(match_XStringSet_XStringViews, 10),
	CALLMETHOD_DEF(vmatch_PDict3Parts_XStringSet, 11),
	CALLMETHOD_DEF(vmatch_XStringSet_XStringSet, 10),

/* align_utils.c */
	CALLMETHOD_DEF(PairwiseAlignedXStringSet_nmatch, 4),
	CALLMETHOD_DEF(AlignedXStringSet_nchar, 1),
	CALLMETHOD_DEF(AlignedXStringSet_align_aligned, 2),
	CALLMETHOD_DEF(PairwiseAlignedFixedSubject_align_aligned, 3),
	CALLMETHOD_DEF(align_compareStrings, 6),

/* pmatchPattern.c */
	CALLMETHOD_DEF(lcprefix, 6),
	CALLMETHOD_DEF(lcsuffix, 6),

/* align_pairwiseAlignment.c */
	CALLMETHOD_DEF(XStringSet_align_pairwiseAlignment, 14),
	CALLMETHOD_DEF(XStringSet_align_distance, 12),

/* align_needwunsQS.c */
	CALLMETHOD_DEF(align_needwunsQS, 7),

/* strutils.c (belonged originally to old matchprobes package) */
	CALLMETHOD_DEF(MP_rna_revcomp, 1),
	CALLMETHOD_DEF(MP_dna_revcomp, 1),
	CALLMETHOD_DEF(MP_revstring, 1),
	CALLMETHOD_DEF(MP_complementSeq, 3),
	CALLMETHOD_DEF(MP_basecontent, 2),
	CALLMETHOD_DEF(MP_longestConsecutive, 2),

/* matchprobes.c (belonged originally to old matchprobes package) */
	CALLMETHOD_DEF(MP_matchprobes, 3),

	{NULL, NULL, 0}
};

void R_init_Biostrings(DllInfo *info)
{
	/* Lots of code around assumes that sizeof(Rbyte) == sizeof(char) */
	if (sizeof(Rbyte) != sizeof(char))
		error("sizeof(Rbyte) != sizeof(char)");
	R_registerRoutines(info, cMethods, NULL, NULL, NULL);
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);

/* XString_class.c */
	REGISTER_CCALLABLE(_DNAencode);
	REGISTER_CCALLABLE(_DNAdecode);
	REGISTER_CCALLABLE(_RNAencode);
	REGISTER_CCALLABLE(_RNAdecode);

/* XStringSet_class.c */
	REGISTER_CCALLABLE(_get_XStringSet_length);
	REGISTER_CCALLABLE(_get_XStringSet_xsbaseclassname);
	REGISTER_CCALLABLE(_cache_XStringSet);
	REGISTER_CCALLABLE(_get_cachedXStringSet_length);
	REGISTER_CCALLABLE(_get_cachedXStringSet_elt);
	REGISTER_CCALLABLE(_set_XStringSet_names);

/* match_reporting.c */
	REGISTER_CCALLABLE(_init_match_reporting);
	REGISTER_CCALLABLE(_set_active_PSpair);
	REGISTER_CCALLABLE(_set_match_shift);
	REGISTER_CCALLABLE(_report_match);
	REGISTER_CCALLABLE(_drop_reported_matches);
	REGISTER_CCALLABLE(_get_match_count);
	REGISTER_CCALLABLE(_reported_matches_asSEXP);

/* MIndex_class.c */
	REGISTER_CCALLABLE(_cache_MIndex);
	REGISTER_CCALLABLE(_get_cachedMIndex_length);
	REGISTER_CCALLABLE(_get_cachedMIndex_elt_width0);
	REGISTER_CCALLABLE(_get_cachedMIndex_elt);

/* match_pattern_boyermoore.c */
	REGISTER_CCALLABLE(_match_pattern_boyermoore);

	return;
}

