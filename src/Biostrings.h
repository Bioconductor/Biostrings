#include "../inst/include/Biostrings_defines.h"
#include <string.h>

#define INIT_STATIC_SYMBOL(NAME) \
{ \
	if (NAME ## _symbol == NULL) \
		NAME ## _symbol = install(# NAME); \
}


/* utils.c */

void _init_ByteTrTable_with_lkup(
	ByteTrTable *byte_tr_table,
	SEXP lkup
);

SEXP _new_lkup_from_ByteTrTable(const ByteTrTable *byte_tr_table);

void _init_byte2offset_with_INTEGER(
	ByteTrTable *byte2offset,
	SEXP bytes,
	int error_on_dup
);

void _init_byte2offset_with_Chars_holder(
	ByteTrTable *byte2offset,
	const Chars_holder *seq,
	const BytewiseOpTable *bytewise_match_table
);

TwobitEncodingBuffer _new_TwobitEncodingBuffer(
	SEXP base_codes,
	int buflength,
	int endianness
);

void _reset_twobit_signature(TwobitEncodingBuffer *teb);

int _shift_twobit_signature(
	TwobitEncodingBuffer *teb,
	char c
);

int _get_twobit_signature(
	TwobitEncodingBuffer *teb,
	const Chars_holder *seq
);

int _get_twobit_signature_at(
	TwobitEncodingBuffer *teb,
	const Chars_holder *seq,
	const int *at,
	int at_length
);


/* RoSeqs_utils.c */

RoSeqs _alloc_RoSeqs(int nelt);


/* XString_class.c */

const ByteTrTable *get_enc_byte2code(const char *classname);

const ByteTrTable *get_dec_byte2code(const char *classname);

SEXP init_DNAlkups(SEXP enc_lkup, SEXP dec_lkup);

char _DNAencode(char c);

char _DNAdecode(char code);

SEXP init_RNAlkups(SEXP enc_lkup, SEXP dec_lkup);

char _RNAencode(char c);

char _RNAdecode(char code);

void _copy_CHARSXP_to_Chars_holder(
	Chars_holder *dest,
	SEXP src,
	int start_in_src,
	const int *lkup,
	int lkup_length
);

SEXP _new_CHARSXP_from_Chars_holder(
	const Chars_holder *x,
	SEXP lkup
);

SEXP new_XString_from_CHARACTER(
	SEXP classname,
	SEXP x,
	SEXP start,
	SEXP width,
	SEXP lkup
);

SEXP new_CHARACTER_from_XString(
	SEXP x,
	SEXP lkup
);


/* XStringSet_class.c */

int _get_XStringSet_length(SEXP x);

SEXP _get_XStringSet_width(SEXP x);

const char *_get_XStringSet_xsbaseclassname(SEXP x);

XStringSet_holder _hold_XStringSet(SEXP x);

int _get_length_from_XStringSet_holder(const XStringSet_holder *x_holder);

Chars_holder _get_elt_from_XStringSet_holder(
	const XStringSet_holder *x_holder,
	int i
);

XStringSet_holder _get_linear_subset_from_XStringSet_holder(
	const XStringSet_holder *x_holder,
	int offset,
	int length
);

void _set_XStringSet_names(
	SEXP x,
	SEXP names
);

SEXP new_XStringSet_from_CHARACTER(
	SEXP classname,
	SEXP element_type,
	SEXP x,
	SEXP start,
	SEXP width,
	SEXP lkup
);

SEXP new_CHARACTER_from_XStringSet(
	SEXP x,
	SEXP lkup
);

RoSeqs _new_RoSeqs_from_XStringSet(
	int nelt,
	SEXP x
);

SEXP XStringSet_unlist(SEXP x);


/* XStringSetList_class.c */

XStringSetList_holder _hold_XStringSetList(SEXP x);

int _get_length_from_XStringSetList_holder(
	const XStringSetList_holder *x_holder
);

XStringSet_holder _get_elt_from_XStringSetList_holder(
	const XStringSetList_holder *x_holder,
	int i
);


/* xscat.c */

SEXP XString_xscat(SEXP args);

SEXP XStringSet_xscat(SEXP args);


/* XStringSet_io.c */

SEXP fasta_index(
	SEXP filexp_list,
	SEXP nrec,
	SEXP skip,
	SEXP seek_first_rec,
	SEXP lkup
);

SEXP read_XStringSet_from_fasta_blocks(
	SEXP seqlength,
	SEXP filexp_list,
	SEXP nrec_list,
	SEXP offset_list,
	SEXP elementType,
	SEXP lkup
);

SEXP write_XStringSet_to_fasta(
	SEXP x,
	SEXP filexp_list,
	SEXP width,
	SEXP lkup
);

SEXP fastq_geometry(
	SEXP filexp_list,
	SEXP nrec,
	SEXP skip,
	SEXP seek_first_rec
);

SEXP read_XStringSet_from_fastq(
	SEXP filexp_list,
	SEXP nrec,
	SEXP skip,
	SEXP seek_first_rec,
	SEXP use_names,
	SEXP elementType,
	SEXP lkup
);

SEXP write_XStringSet_to_fastq(
	SEXP x,
	SEXP filexp_list,
	SEXP qualities,
	SEXP lkup
);


/* letter_frequency.c */

SEXP XString_letter_frequency(
	SEXP x,
	SEXP codes,
	SEXP with_other
);

SEXP XStringSet_letter_frequency(
	SEXP x,
	SEXP collapse,
	SEXP codes,
	SEXP with_other
);

SEXP XString_letterFrequencyInSlidingView(
	SEXP x,
	SEXP view_width,
	SEXP single_codes,
	SEXP colmap,
	SEXP colnames
);

SEXP XStringSet_letterFrequency(
	SEXP x,
	SEXP single_codes,
	SEXP colmap,
	SEXP colnames,
	SEXP collapse
);

SEXP XString_oligo_frequency(
	SEXP x,
	SEXP width,
	SEXP step,
	SEXP as_prob,
	SEXP as_array,
	SEXP fast_moving_side,
	SEXP with_labels,
	SEXP base_codes
);

SEXP XStringSet_oligo_frequency(
	SEXP x,
	SEXP width,
	SEXP step,
	SEXP as_prob,
	SEXP as_array,
	SEXP fast_moving_side,
	SEXP with_labels,
	SEXP simplify_as,
	SEXP base_codes
);

SEXP XStringSet_nucleotide_frequency_at(
	SEXP x,
	SEXP at,
	SEXP as_prob,
	SEXP as_array,
	SEXP fast_moving_side,
	SEXP with_labels,
	SEXP base_codes
);

SEXP XStringSet_consensus_matrix(
	SEXP x,
	SEXP shift,
	SEXP width,
	SEXP with_other,
	SEXP codes
);

SEXP XString_two_way_letter_frequency(
        SEXP x,
        SEXP y,
        SEXP x_codes,
        SEXP y_codes,
        SEXP with_other
);

SEXP XStringSet_two_way_letter_frequency(
        SEXP x,
        SEXP y,
        SEXP collapse,
        SEXP x_codes,
        SEXP y_codes,
        SEXP with_other
);

SEXP XStringSet_two_way_letter_frequency_by_quality(
        SEXP x,
        SEXP y,
        SEXP x_quality,
        SEXP y_quality,
        SEXP codes,
        SEXP quality_codes,
        SEXP with_other
);

/* gtestsim.c */

void gtestsim(
	int *nrow,
	int *ncol,
	int *nrowt,
	int *ncolt,
	int *n,
	int *b,
	double *expected,
	int *observed,
	double *fact,
	int *jwork,
	double *results
);


/* translate.c */

SEXP DNAStringSet_translate(
	SEXP x,
	SEXP skip_code,
	SEXP dna_codes,
	SEXP lkup,
	SEXP init_lkup,
	SEXP if_non_ambig,
	SEXP if_ambig
);

/* replaceAt.c */

SEXP XString_replaceAt(
	SEXP x,
	SEXP at,
	SEXP value
);

SEXP XStringSet_replaceAt(
	SEXP x,
	SEXP at,
	SEXP value
);

/* replace_letter_at.c */

SEXP XString_replace_letter_at(
	SEXP x,
	SEXP at,
	SEXP letter,
	SEXP lkup,
	SEXP if_not_extending,
	SEXP verbose
);

SEXP XString_inplace_replace_letter_at(
	SEXP x,
	SEXP at,
	SEXP letter,
	SEXP lkup
);


/* inject_code.c */

SEXP XString_inject_code(
	SEXP x,
	SEXP start,
	SEXP width,
	SEXP code
);


/* unstrsplit_methods.c */

SEXP XStringSetList_unstrsplit(
	SEXP x,
	SEXP sep,
	SEXP seqtype
);


/* SparseList_utils.c */

SEXP _SparseList_int2symb(int symb_as_int);

int _SparseList_symb2int(SEXP symbol);

SEXP _get_val_from_env(
	SEXP symbol,
	SEXP env,
	int error_on_unbound_value
);

SEXP _get_val_from_SparseList(
	int symb_as_int,
	SEXP env,
	int error_on_unbound_value
);

int _get_int_from_SparseList(
	int symb_as_int,
	SEXP env
);

void _set_env_from_IntAE(
	SEXP env,
	const IntAE *int_ae
);

void _set_env_from_IntAEAE(
	SEXP env,
	const IntAEAE *int_aeae
);


/* match_reporting.c */

int _get_match_storing_code(const char *ms_mode);

MatchBuf _new_MatchBuf(
	int ms_code,
	int nPSpair
);

void _MatchBuf_report_match(
	MatchBuf *match_buf,
	int PSpair_id,
	int start,
	int width
);

void _MatchBuf_flush(MatchBuf *match_buf);

void _MatchBuf_append_and_flush(
	MatchBuf *match_buf1,
	MatchBuf *match_buf2,
	int view_offset
);

SEXP _MatchBuf_which_asINTEGER(const MatchBuf *match_buf);

SEXP _MatchBuf_counts_asINTEGER(const MatchBuf *match_buf);

SEXP _MatchBuf_starts_asLIST(const MatchBuf *match_buf);

SEXP _MatchBuf_ends_asLIST(const MatchBuf *match_buf);

SEXP _MatchBuf_as_Ranges(const MatchBuf *match_buf);

SEXP _MatchBuf_as_SEXP(
	const MatchBuf *match_buf,
	SEXP env
);

void _init_match_reporting(const char *ms_mode, int nPSpair);

void _set_active_PSpair(int PSpair_id);

void _set_match_shift(int shift);

void _report_match(int start, int width);

void _drop_reported_matches();

int _get_match_count();

SEXP _reported_matches_asSEXP();

MatchBuf *_get_internal_match_buf();


/* MIndex_class.c */

MIndex_holder _hold_MIndex(SEXP x);

int _get_length_from_MIndex_holder(const MIndex_holder *x_holder);

int _get_width0_elt_from_MIndex_holder(
	const MIndex_holder *x_holder,
	int i
);

IRanges_holder _get_elt_from_MIndex_holder(
	const MIndex_holder *x_holder,
	int i
);

SEXP ByPos_MIndex_endIndex(
	SEXP x_high2low,
	SEXP x_ends,
	SEXP x_width0
);

SEXP SparseMIndex_endIndex(
	SEXP x_ends_envir,
	SEXP x_width0,
	SEXP x_names,
	SEXP all_names
);

SEXP ByPos_MIndex_combine(SEXP ends_listlist);


/* lowlevel_matching.c */

void _init_bytewise_match_tables();

const BytewiseOpTable *_select_bytewise_match_table(int fixedP, int fixedS);

int _nmismatch_at_Pshift(
	const Chars_holder *P,
	const Chars_holder *S,
	int Pshift,
	int max_nmis,
	const BytewiseOpTable *bytewise_match_table
);

int _nedit_for_Ploffset(
	const Chars_holder *P,
	const Chars_holder *S,
	int Ploffset,
	int max_nedit,
	int loose_Ploffset,
	int *min_width,
	const BytewiseOpTable *bytewise_match_table
);

int _nedit_for_Proffset(
	const Chars_holder *P,
	const Chars_holder *S,
	int Proffset,
	int max_nedit,
	int loose_Proffset,
	int *min_width,
	const BytewiseOpTable *bytewise_match_table
);

SEXP XString_match_pattern_at(
	SEXP pattern,
	SEXP subject,
	SEXP at,
	SEXP at_type,
	SEXP max_mismatch,
	SEXP min_mismatch,
	SEXP with_indels,
	SEXP fixed,
	SEXP ans_type,
	SEXP auto_reduce_pattern
);

SEXP XStringSet_vmatch_pattern_at(
	SEXP pattern,
	SEXP subject,
	SEXP at,
	SEXP at_type,
	SEXP max_mismatch,
	SEXP min_mismatch,
	SEXP with_indels,
	SEXP fixed,
	SEXP ans_type,
	SEXP auto_reduce_pattern
);

SEXP XStringSet_dist_hamming(SEXP x);


/* match_pattern_boyermoore.c */

int _match_pattern_boyermoore(
	const Chars_holder *P,
	const Chars_holder *S,
	int nfirstmatches,
	int walk_backward
);


/* match_pattern_shiftor.c */

SEXP bits_per_long();

void _match_pattern_shiftor(
	const Chars_holder *P,
	const Chars_holder *S,
	int max_nmis,
	int fixedP,
	int fixedS
);


/* match_pattern_indels.c */

void _match_pattern_indels(
	const Chars_holder *P,
	const Chars_holder *S,
	int max_nmis,
	int fixedP,
	int fixedS
);


/* match_pattern.c */

void _match_pattern_XString(
	const Chars_holder *P,
	const Chars_holder *S,
	SEXP max_mismatch,
	SEXP min_mismatch,
	SEXP with_indels,
	SEXP fixed,
	const char *algo
);

void _match_pattern_XStringViews(
	const Chars_holder *P,
	const Chars_holder *S,
	SEXP views_start,
	SEXP views_width,
	SEXP max_mismatch,
	SEXP min_mismatch,
	SEXP with_indels,
	SEXP fixed,
	const char *algo
);

SEXP XString_match_pattern(
	SEXP pattern,
	SEXP subject,
	SEXP max_mismatch,
	SEXP min_mismatch,
	SEXP with_indels,
	SEXP fixed,
	SEXP algorithm,
	SEXP count_only
);

SEXP XStringViews_match_pattern(
	SEXP pattern,
	SEXP subject,
	SEXP views_start,
	SEXP views_width,
	SEXP max_mismatch,
	SEXP min_mismatch,
	SEXP with_indels,
	SEXP fixed,
	SEXP algorithm,
	SEXP count_only
);

SEXP XStringSet_vmatch_pattern(
	SEXP pattern,
	SEXP subject,
	SEXP max_mismatch,
	SEXP min_mismatch,
	SEXP with_indels,
	SEXP fixed,
	SEXP algorithm,
	SEXP ms_mode
);


/* match_PWM.c */

SEXP PWM_score_starting_at(
	SEXP pwm,
	SEXP subject,
	SEXP starting_at,
	SEXP base_codes
);

SEXP XString_match_PWM(
	SEXP pwm,
	SEXP subject,
	SEXP min_score,
	SEXP count_only,
	SEXP base_codes
);

SEXP XStringViews_match_PWM(
	SEXP pwm,
	SEXP subject,
	SEXP views_start,
	SEXP views_width,
	SEXP min_score,
	SEXP count_only,
	SEXP base_codes
);


/* find_palindromes.c */

SEXP find_palindromes(
	SEXP x,
	SEXP min_armlength,
	SEXP max_looplength,
	SEXP max_mismatch,
	SEXP L2R_lkup
);

SEXP palindrome_arm_length(
	SEXP x,
	SEXP max_mismatch,
	SEXP L2R_lkup
);


/* BitMatrix.c */

void _BitCol_set_val(
	BitCol *bitcol,
	BitWord val
);

BitCol _new_BitCol(
	int nbit,
	BitWord val
);

int _BitCol_get_bit(
	const BitCol *bitcol,
	int i
);

void _BitCol_set_bit(
	BitCol *bitcol,
	int i,
	int bit
);

void _BitCol_A_gets_BimpliesA(
	BitCol *A,
	const BitCol *B
);

BitCol _BitMatrix_get_col(
	const BitMatrix *bitmat,
	int j
);

void _BitMatrix_set_col(
	BitMatrix *bitmat,
	int j,
	const BitCol *bitcol
);

void _BitMatrix_set_val(
	BitMatrix *bitmat,
	BitWord val
);

BitMatrix _new_BitMatrix(
	int nrow,
	int ncol,
	BitWord val
);

int _BitMatrix_get_bit(
	const BitMatrix *bitmat,
	int i,
	int j
);

void _BitMatrix_set_bit(
	BitMatrix *bitmat,
	int i,
	int j,
	int bit
);

void _BitMatrix_Rrot1(BitMatrix *bitmat);

void _BitMatrix_grow1rows(
	BitMatrix *bitmat,
	const BitCol *bitcol
);


/* PreprocessedTB_class.c */

SEXP _get_PreprocessedTB_tb(SEXP x);

SEXP _get_PreprocessedTB_dups(SEXP x);

SEXP _get_PreprocessedTB_base_codes(SEXP x);

int _get_PreprocessedTB_length(SEXP x);

int _get_PreprocessedTB_width(SEXP x);

SEXP _get_PreprocessedTB_low2high(SEXP x);

SEXP _get_Twobit_sign2pos_tag(SEXP x);

SEXP _get_ACtree2_nodebuf_ptr(SEXP x);

SEXP _get_ACtree2_nodeextbuf_ptr(SEXP x);

void _init_ppdups_buf(int length);

void _report_ppdup(
	int poffset,
	int P_id
);

SEXP _get_ppdups_buf_asINTEGER();


/* match_pdict_utils.c */

TBMatchBuf _new_TBMatchBuf(
	int tb_length,
	int tb_width,
	const int *head_widths,
	const int *tail_widths
);

void _TBMatchBuf_report_match(
	TBMatchBuf *buf,
	int PSpair_id,
	int end
);

void _TBMatchBuf_flush(TBMatchBuf *buf);

MatchPDictBuf _new_MatchPDictBuf(
	SEXP matches_as,
	int tb_length,
	int tb_width,
	const int *head_widths,
	const int *tail_widths
);

void _MatchPDictBuf_report_match(
	MatchPDictBuf *buf,
	int PSpair_id,
	int tb_end
);

void _MatchPDictBuf_flush(MatchPDictBuf *buf);

void _MatchPDictBuf_append_and_flush(
	MatchBuf *buf1,
	MatchPDictBuf *buf2,
	int view_offset
);

HeadTail _new_HeadTail(
	SEXP pdict_head,
	SEXP pdict_tail,
	SEXP pptb,
	SEXP max_mismatch,
	SEXP fixed,
	int with_ppheadtail
);

void _match_pdict_flanks_at(
	int key0,
	SEXP low2high,
	HeadTail *headtail,
	const Chars_holder *S,
	int tb_end,
	int max_nmis,
	int min_nmis,
	int fixedP,
	int fixedS,
	MatchPDictBuf *matchpdict_buf
);

void _match_pdict_all_flanks(
	SEXP low2high,
	HeadTail *headtail,
	const Chars_holder *S,
	int max_nmis,
	int min_nmis,
	int fixedP,
	int fixedS,
	MatchPDictBuf *matchpdict_buf
);


/* match_pdict_Twobit.c */

SEXP build_Twobit(
	SEXP tb,
	SEXP pp_exclude,
	SEXP base_codes
);

void _match_Twobit(
	SEXP pptb,
	const Chars_holder *S,
	int fixedS,
	TBMatchBuf *tb_matches
);


/* BAB_class.c */

SEXP IntegerBAB_new(SEXP max_nblock);

int *_get_BAB_nblock_ptr(SEXP x);

int *_get_BAB_lastblock_nelt_ptr(SEXP x);

SEXP _get_BAB_blocks(SEXP x);

SEXP _IntegerBAB_addblock(
	SEXP x,
	int block_length
);


/* match_pdict_ACtree2.c */

SEXP ACtree2_nodebuf_max_nblock();

SEXP ACtree2_nodeextbuf_max_nblock();

SEXP ACtree2_nnodes(SEXP pptb);

SEXP ACtree2_print_nodes(SEXP pptb);

SEXP ACtree2_summary(SEXP pptb);

SEXP ACtree2_build(
	SEXP tb,
	SEXP pp_exclude,
	SEXP base_codes,
	SEXP nodebuf_ptr,
	SEXP nodeextbuf_ptr
);

SEXP ACtree2_has_all_flinks(SEXP pptb);

SEXP ACtree2_compute_all_flinks(SEXP pptb);

void _match_tbACtree2(
	SEXP pptb,
	const Chars_holder *S,
	int fixedS,
	TBMatchBuf *tb_matches
);

void _match_pdictACtree2(
	SEXP pptb,
	HeadTail *headtail,
	const Chars_holder *S,
	int max_nmis,
	int min_nmis,
	int fixedP,
	int fixedS,
	MatchPDictBuf *matchpdict_buf
);


/* match_pdict.c */

SEXP match_PDict3Parts_XString(
	SEXP pptb,
	SEXP pdict_head,
	SEXP pdict_tail,
	SEXP subject,
	SEXP max_mismatch,
	SEXP min_mismatch,
	SEXP fixed,
	SEXP matches_as,
	SEXP envir
);

SEXP match_XStringSet_XString(
	SEXP pattern,
	SEXP subject,
	SEXP max_mismatch,
	SEXP min_mismatch,
	SEXP with_indels,
	SEXP fixed,
	SEXP algorithm,
	SEXP matches_as,
	SEXP envir
);

SEXP match_PDict3Parts_XStringViews(
	SEXP pptb,
	SEXP pdict_head,
	SEXP pdict_tail,
	SEXP subject,
	SEXP views_start,
	SEXP views_width,
	SEXP max_mismatch,
	SEXP min_mismatch,
	SEXP fixed,
	SEXP matches_as,
	SEXP envir
);

SEXP match_XStringSet_XStringViews(
	SEXP pattern,
	SEXP subject,
	SEXP views_start,
	SEXP views_width,
	SEXP max_mismatch,
	SEXP min_mismatch,
	SEXP with_indels,
	SEXP fixed,
	SEXP algorithm,
	SEXP matches_as,
	SEXP envir
);

SEXP vmatch_PDict3Parts_XStringSet(
	SEXP pptb,
	SEXP pdict_head,
	SEXP pdict_tail,
	SEXP subject,
	SEXP max_mismatch,
	SEXP min_mismatch,
	SEXP fixed,
	SEXP collapse,
	SEXP weight,
	SEXP matches_as,
	SEXP envir
);

SEXP vmatch_XStringSet_XStringSet(
	SEXP pattern,
	SEXP subject,
	SEXP max_mismatch,
	SEXP min_mismatch,
	SEXP with_indels,
	SEXP fixed,
	SEXP algorithm,
	SEXP collapse,
	SEXP weight,
	SEXP matches_as,
	SEXP envir
);


/* align_utils.c */

SEXP PairwiseAlignments_nmatch(
	SEXP nchar,
	SEXP nmismatch,
	SEXP ninsertion,
	SEXP ndeletion
);

SEXP AlignedXStringSet_nchar(SEXP alignedXStringSet);

SEXP AlignedXStringSet_align_aligned(
	SEXP alignedXStringSet,
	SEXP gapCode
);

SEXP PairwiseAlignmentsSingleSubject_align_aligned(
	SEXP alignment,
	SEXP gapCode,
	SEXP endgapCode
);

SEXP align_compareStrings(
	SEXP patternStrings,
	SEXP subjectStrings,
	SEXP maxNChar,
	SEXP insertionCode,
	SEXP deletionCode,
	SEXP mismatchCode
);


/* pmatchPattern.c */

SEXP lcprefix(
	SEXP s1_xp,
	SEXP s1_offset,
	SEXP s1_length,
	SEXP s2_xp,
	SEXP s2_offset,
	SEXP s2_length
);

SEXP lcsuffix(
	SEXP s1_xp,
	SEXP s1_offset,
	SEXP s1_length,
	SEXP s2_xp,
	SEXP s2_offset,
	SEXP s2_length
);


/* align_pairwiseAlignment.c */

SEXP XStringSet_align_pairwiseAlignment(
	SEXP pattern,
	SEXP subject,
	SEXP type,
	SEXP typeCode,
	SEXP scoreOnly,
	SEXP gapOpening,
	SEXP gapExtension,
	SEXP useQuality,
	SEXP substitutionArray,
	SEXP substitutionArrayDim,
	SEXP substitutionLookupTable,
	SEXP fuzzyMatrix,
	SEXP fuzzyMatrixDim,
	SEXP fuzzyLookupTable
);

SEXP XStringSet_align_distance(
	SEXP string,
	SEXP type,
	SEXP typeCode,
	SEXP gapOpening,
	SEXP gapExtension,
	SEXP useQuality,
	SEXP substitutionArray,
	SEXP substitutionArrayDim,
	SEXP substitutionLookupTable,
	SEXP fuzzyMatrix,
	SEXP fuzzyMatrixDim,
	SEXP fuzzyLookupTable
);


/* align_needwunsQS.c */

SEXP align_needwunsQS(
	SEXP s1,
	SEXP s2,
	SEXP mat,
	SEXP mat_nrow,
	SEXP lkup,
	SEXP gap_cost,
	SEXP gap_code
);


/* strutils.c (belonged originally to old matchprobes package) */

char compbase(char c);
SEXP MP_longestConsecutive(SEXP x, SEXP letter);


/* matchprobes.c (belonged originally to old matchprobes package) */

SEXP MP_matchprobes(SEXP query, SEXP records, SEXP probepos);

