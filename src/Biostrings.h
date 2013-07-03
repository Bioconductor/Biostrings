#include "../inst/include/Biostrings_defines.h"
#include <string.h>

#define DEBUG_BIOSTRINGS 1

#define INIT_STATIC_SYMBOL(NAME) \
{ \
	if (NAME ## _symbol == NULL) \
		NAME ## _symbol = install(# NAME); \
}


/* utils.c */

SEXP debug_utils();

void _init_ByteTrTable_with_lkup(
	ByteTrTable byte2code,
	SEXP lkup
);

SEXP _new_lkup_from_ByteTrTable(const ByteTrTable *byte2code);

void _init_byte2offset_with_INTEGER(
	ByteTrTable byte2offset,
	SEXP bytes,
	int error_on_dup
);

void _init_byte2offset_with_cachedCharSeq(
	ByteTrTable byte2offset,
	const cachedCharSeq *seq,
	int error_on_dup
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
	const cachedCharSeq *seq
);

int _get_twobit_signature_at(
	TwobitEncodingBuffer *teb,
	const cachedCharSeq *seq,
	const int *at,
	int at_length
);


/* io_utils.c */

SEXP new_input_ExternalFilePtr(SEXP filepath);

SEXP new_output_ExternalFilePtr(SEXP filepath, SEXP append);

SEXP ExternalFilePtr_close(SEXP x);

int delete_trailing_LF_or_CRLF(
	const char *buf,
	int size
);


/* RoSeqs_utils.c */

SEXP debug_RoSeqs_utils();

RoSeqs _alloc_RoSeqs(int nelt);


/* XString_class.c */

SEXP debug_XString_class();

const ByteTrTable *get_enc_byte2code(const char *classname);

const ByteTrTable *get_dec_byte2code(const char *classname);

SEXP init_DNAlkups(SEXP enc_lkup, SEXP dec_lkup);

char _DNAencode(char c);

char _DNAdecode(char code);

SEXP init_RNAlkups(SEXP enc_lkup, SEXP dec_lkup);

char _RNAencode(char c);

char _RNAdecode(char code);

void _copy_CHARSXP_to_cachedCharSeq(
	cachedCharSeq *dest,
	SEXP src,
	int start_in_src,
	const int *lkup,
	int lkup_length
);

SEXP _new_CHARSXP_from_cachedCharSeq(
	const cachedCharSeq *x,
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

SEXP debug_XStringSet_class();

int _get_XStringSet_length(SEXP x);

SEXP _get_XStringSet_width(SEXP x);

const char *_get_XStringSet_xsbaseclassname(SEXP x);

cachedXStringSet _cache_XStringSet(SEXP x);

int _get_cachedXStringSet_length(const cachedXStringSet *cached_x);

cachedCharSeq _get_cachedXStringSet_elt(
	const cachedXStringSet *cached_x,
	int i
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


/* xscat.c */

SEXP XString_xscat(SEXP args);

SEXP XStringSet_xscat(SEXP args);


/* XStringSet_io.c */

SEXP debug_XStringSet_io();

SEXP fasta_info(
	SEXP efp_list,
	SEXP nrec,
	SEXP skip,
	SEXP use_names,
	SEXP lkup
);

SEXP read_fasta_in_XStringSet(
	SEXP efp_list,
	SEXP nrec,
	SEXP skip,
	SEXP use_names,
	SEXP elementType,
	SEXP lkup
);

SEXP write_XStringSet_to_fasta(
	SEXP x,
	SEXP efp_list,
	SEXP width,
	SEXP lkup
);

SEXP fastq_geometry(
	SEXP efp_list,
	SEXP nrec,
	SEXP skip
);

SEXP read_fastq_in_XStringSet(
	SEXP efp_list,
	SEXP nrec,
	SEXP skip,
	SEXP use_names,
	SEXP elementType,
	SEXP lkup
);

SEXP write_XStringSet_to_fastq(
	SEXP x,
	SEXP efp_list,
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
	SEXP base_codes,
	SEXP lkup,
	SEXP skipcode
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


/* SparseList_utils.c */

SEXP debug_SparseList_utils();

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

SEXP debug_match_reporting();

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

SEXP _MatchBuf_as_MIndex(const MatchBuf *match_buf);

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

SEXP debug_MIndex_class();

cachedMIndex _cache_MIndex(SEXP x);

int _get_cachedMIndex_length(const cachedMIndex *cached_x);

int _get_cachedMIndex_elt_width0(const cachedMIndex *cached_x, int i);

cachedIRanges _get_cachedMIndex_elt(const cachedMIndex *cached_x, int i);

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

SEXP debug_lowlevel_matching();

int (*_selected_nmismatch_at_Pshift_fun)(
	const cachedCharSeq *P,
	const cachedCharSeq *S,
	int Pshift,
	int max_nmis
);

void _select_nmismatch_at_Pshift_fun(
	int fixedP,
	int fixedS
);

int _nedit_for_Ploffset(
	const cachedCharSeq *P,
	const cachedCharSeq *S,
	int Ploffset,
	int max_nedit,
	int loose_Ploffset,
	int *min_width
);

int _nedit_for_Proffset(
	const cachedCharSeq *P,
	const cachedCharSeq *S,
	int Proffset,
	int max_nedit,
	int loose_Proffset,
	int *min_width
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

SEXP debug_match_pattern_boyermoore();

int _match_pattern_boyermoore(
	const cachedCharSeq *P,
	const cachedCharSeq *S,
	int nfirstmatches,
	int walk_backward
);


/* match_pattern_shiftor.c */

SEXP debug_match_pattern_shiftor();

SEXP bits_per_long();

void _match_pattern_shiftor(
	const cachedCharSeq *P,
	const cachedCharSeq *S,
	int max_nmis,
	int fixedP,
	int fixedS
);


/* match_pattern_indels.c */

SEXP debug_match_pattern_indels();

void _match_pattern_indels(
	const cachedCharSeq *P,
	const cachedCharSeq *S,
	int max_nmis,
	int fixedP,
	int fixedS
);


/* match_pattern.c */

SEXP debug_match_pattern();

void _match_pattern_XString(
	const cachedCharSeq *P,
	const cachedCharSeq *S,
	SEXP max_mismatch,
	SEXP min_mismatch,
	SEXP with_indels,
	SEXP fixed,
	const char *algo
);

void _match_pattern_XStringViews(
	const cachedCharSeq *P,
	const cachedCharSeq *S,
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
	SEXP count_only
);


/* match_BOC.c */

SEXP debug_match_BOC();

SEXP match_BOC_preprocess(
	SEXP s_xp,
	SEXP s_offset,
	SEXP s_length,
	SEXP p_length,
	SEXP code1,
	SEXP code2,
	SEXP code3,
	SEXP code4,
	SEXP buf1_xp,
	SEXP buf2_xp,
	SEXP buf3_xp,
	SEXP pre4buf_xp
);

SEXP match_BOC_exact(
	SEXP p_xp,
	SEXP p_offset,
	SEXP p_length,
	SEXP s_xp,
	SEXP s_offset,
	SEXP s_length,
	SEXP code1,
	SEXP code2,
	SEXP code3,
	SEXP code4,
	SEXP buf1_xp,
	SEXP buf2_xp,
	SEXP buf3_xp,
	SEXP pre4buf_xp,
	SEXP stats,
	SEXP count_only
);


/* match_BOC2.c */

SEXP debug_match_BOC2();

SEXP match_BOC2_preprocess(
	SEXP s_xp,
	SEXP s_offset,
	SEXP s_length,
	SEXP p_length,
	SEXP code1,
	SEXP code2,
	SEXP code3,
	SEXP code4,
	SEXP buf_xp
);

SEXP match_BOC2_exact(
	SEXP p_xp,
	SEXP p_offset,
	SEXP p_length,
	SEXP s_xp,
	SEXP s_offset,
	SEXP s_length,
	SEXP code1,
	SEXP code2,
	SEXP code3,
	SEXP code4,
	SEXP buf_xp,
	SEXP stats,
	SEXP count_only
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

SEXP debug_find_palindromes();

SEXP find_palindromes(
	SEXP s_xp,
	SEXP s_offset,
	SEXP s_length,
	SEXP min_armlength,
	SEXP max_looplength,
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

SEXP debug_BitMatrix();


/* PreprocessedTB_class.c */

SEXP debug_PreprocessedTB_class();

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

SEXP debug_match_pdict_utils();

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
	const cachedCharSeq *S,
	int tb_end,
	int max_nmis,
	int min_nmis,
	int fixedP,
	MatchPDictBuf *matchpdict_buf
);

void _match_pdict_all_flanks(
	SEXP low2high,
	HeadTail *headtail,
	const cachedCharSeq *S,
	int max_nmis,
	int min_nmis,
	MatchPDictBuf *matchpdict_buf
);


/* match_pdict_Twobit.c */

SEXP debug_match_pdict_Twobit();

SEXP build_Twobit(
	SEXP tb,
	SEXP pp_exclude,
	SEXP base_codes
);

void _match_Twobit(
	SEXP pptb,
	const cachedCharSeq *S,
	int fixedS,
	TBMatchBuf *tb_matches
);


/* BAB_class.c */

SEXP debug_BAB_class();

SEXP IntegerBAB_new(SEXP max_nblock);

int *_get_BAB_nblock_ptr(SEXP x);

int *_get_BAB_lastblock_nelt_ptr(SEXP x);

SEXP _get_BAB_blocks(SEXP x);

SEXP _IntegerBAB_addblock(
	SEXP x,
	int block_length
);


/* match_pdict_ACtree2.c */

SEXP debug_match_pdict_ACtree2();

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
	const cachedCharSeq *S,
	int fixedS,
	TBMatchBuf *tb_matches
);

void _match_pdictACtree2(
	SEXP pptb,
	HeadTail *headtail,
	const cachedCharSeq *S,
	int max_nmis,
	int min_nmis,
	int fixedP,
	int fixedS,
	MatchPDictBuf *matchpdict_buf
);


/* match_pdict.c */

SEXP debug_match_pdict();

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

