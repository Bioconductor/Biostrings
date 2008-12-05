#include "../inst/include/Biostrings_defines.h"
#include <string.h>

#define DEBUG_BIOSTRINGS 1


/* utils.c */

SEXP debug_utils();

int fgets_rtrimmed(
	char *s,
	int size,
	FILE *stream
);

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

void _init_byte2offset_with_RoSeq(
	ByteTrTable byte2offset,
	const RoSeq *seq,
	int error_on_dup
);


/* copy_seq.c */

SEXP debug_copy_seq();

void _copy_seq(
	char *dest,
	const char *src,
	size_t n,
	const ByteTrTable *byte2code
);

void _revcopy_seq(
	char *dest,
	const char *src,
	size_t n,
	const ByteTrTable *byte2code
);

void _copy_seq_from_i1i2(
	int i1, int i2,
	char *dest, int dest_length,
	const char *src, int src_length,
	const ByteTrTable *byte2code
);

void _copy_seq_to_i1i2(
	int i1, int i2,
	char *dest, int dest_length,
	const char *src, int src_length,
	const ByteTrTable *byte2code
);

void _copy_seq_from_subset(
	const int *subset, int n,
	char *dest, int dest_length,
	const char *src, int src_length,
	const ByteTrTable *byte2code
);

void _copy_seq_to_subset(
	const int *subset, int n,
	char *dest, int dest_length,
	const char *src, int src_length,
	const ByteTrTable *byte2code
);


/* RoSeq_utils.c */

SEXP debug_RoSeq_utils();

RoSeqs _alloc_RoSeqs(int nelt);

void _narrow_RoSeqs(
	RoSeqs *seqs,
	SEXP start,
	SEXP width
);

SEXP _new_CHARSXP_from_RoSeq(
	const RoSeq *seq,
	SEXP lkup
);

RoSeqs _new_RoSeqs_from_STRSXP(
	int nelt,
	SEXP x
);

SEXP _new_STRSXP_from_RoSeqs(
	const RoSeqs *seqs,
	SEXP lkup
);

SEXP _new_RawPtr_from_RoSeqs(
	const RoSeqs *seqs,
	SEXP lkup
);

SEXP new_RawPtr_from_STRSXP(
	SEXP x,
	SEXP start,
	SEXP width,
	SEXP collapse,
	SEXP lkup
);

RoSeqs _new_RoSeqs_from_CharAEAE(const CharAEAE *char_aeae);

SEXP _new_IRanges_from_RoSeqs(
	const char *classname,
	const RoSeqs *seqs
);

void _write_RoSeq_to_RawPtr(
	SEXP x,
	int offset,
	const RoSeq *seq,
	const ByteTrTable *byte2code
);

void _get_RoSeqs_order(
	const RoSeqs *seqs,
	int *order,
	int base1
);

void _get_RoSeqs_rank(
	const RoSeqs *seqs,
	int *rank
);

void _get_RoSeqs_duplicated(
	const RoSeqs *seqs,
	int *duplicated
);

void _get_RoSeqs_not_duplicated(
	const RoSeqs *seqs,
	int *not_duplicated
);


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

RoSeq _get_XString_asRoSeq(SEXP x);

SEXP new_RawPtr_from_XString(
	SEXP x,
	SEXP start,
	SEXP width,
	SEXP lkup
);

SEXP _new_XString_from_RoSeqs(
	const char *classname,
	const RoSeqs *seqs
);

SEXP _alloc_XString(
	const char *classname,
	int length
);

void _write_RoSeq_to_XString(
	SEXP x,
	int start,
	const RoSeq *seq,
	int encode
);


/* XStringSet_class.c */

SEXP debug_XStringSet_class();

SEXP _get_XStringSet_super(SEXP x);

const char *_get_XStringSet_baseClass(SEXP x);

SEXP _get_XStringSet_ranges(SEXP x);

int _get_XStringSet_length(SEXP x);

SEXP _get_XStringSet_width(SEXP x);

CachedXStringSet _new_CachedXStringSet(SEXP x);

RoSeq _get_CachedXStringSet_elt_asRoSeq(
	CachedXStringSet *x,
	int i
);

RoSeq _get_XStringSet_elt_asRoSeq(
	SEXP x,
	int i
);

RoSeqs _new_RoSeqs_from_XStringSet(
	int nelt,
	SEXP x
);

SEXP _new_XStringSet(
	const char *classname,
	SEXP super,
	SEXP ranges
);

SEXP _new_XStringSet_from_RoSeqs(
	const char *baseClass,
	const RoSeqs *seqs
);

void _set_XStringSet_names(
	SEXP x,
	SEXP names
);

SEXP _alloc_XStringSet(
	const char *baseClass,
	int length,
	int super_length
);

void _write_RoSeq_to_CachedXStringSet_elt(
	CachedXStringSet *x,
	int i,
	const RoSeq *seq,
	int encode
);

void _write_RoSeq_to_XStringSet_elt(
	SEXP x,
	int i,
	const RoSeq *seq,
	int encode
);

SEXP XStringSet_as_STRSXP(
	SEXP x,
	SEXP lkup
);

SEXP XStringSet_order(SEXP x);

SEXP XStringSet_rank(SEXP x);

SEXP XStringSet_duplicated(SEXP x);

SEXP XStringSet_not_duplicated(SEXP x);


/* fasta_io.c */

SEXP debug_fasta_io();

SEXP fasta_info(
	SEXP filepath,
	SEXP use_descs
);

SEXP RawPtr_loadFASTA(
	SEXP rawptr_xp,
	SEXP filepath,
	SEXP collapse,
	SEXP lkup
);


/* char_frequency.c */

SEXP XString_char_frequency(
	SEXP x,
	SEXP codes,
	SEXP with_other
);

SEXP XStringSet_char_frequency(
	SEXP x,
	SEXP codes,
	SEXP with_other,
	SEXP collapse
);

SEXP XStringSet_char_frequency_by_pos(
	SEXP x,
	SEXP codes,
	SEXP with_other
);

SEXP oligonucleotide_frequency(
	SEXP x,
	SEXP base_codes,
	SEXP width,
	SEXP fast_moving_side
);


/* char_translate.c */

SEXP XStringSet_char_translate(
	SEXP x,
	SEXP lkup,
	SEXP reverse
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

SEXP inject_code(
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


/* Dups_utils.c */

SEXP debug_Dups_utils();

SEXP Dups_diff(
	SEXP x_unq2dup,
	SEXP y_dup2unq
);

void _init_dup2unq_buf(int length);

void _report_dup(
	int poffset,
	int P_id
);

SEXP _dup2unq_asINTEGER();


/* match_reporting.c */

SEXP debug_match_reporting();

void _init_match_reporting(SEXP mode);

void _drop_reported_matches();

void _shift_match_on_reporting(int shift);

void _report_match(int start, int width);

SEXP _reported_matches_asSEXP();


/* MIndex_utils.c */

SEXP debug_MIndex_utils();

void _MIndex_init_match_reporting(
	int is_count_only,
	int with_headtail,
	int pdict_L
);

int _MIndex_get_match_reporting_mode();

IntAE *_MIndex_get_match_count();

IntAE *_MIndex_get_match_ends(int key);

IntAE *_MIndex_get_matching_keys();

void _MIndex_report_match(
	int key,
	int end
);

void _MIndex_merge_matches(
	IntAE *global_match_count,
	const IntAEAE *global_match_ends,
	int view_offset
);

SEXP _MIndex_reported_matches_asSEXP(SEXP env);

SEXP ByPos_MIndex_endIndex(
	SEXP x_dup2unq,
	SEXP x_ends,
	SEXP x_width
);

SEXP ByName_MIndex_endIndex(
	SEXP x_ends_envir,
	SEXP x_width,
	SEXP x_names,
	SEXP all_names
);

void ByPos_MIndex_coverage(
	SEXP ends_list,
	SEXP mindex_width,
	SEXP start,
	SEXP ans_xp
);

void ByName_MIndex_coverage(
	SEXP ends_envir,
	SEXP mindex_width,
	SEXP start,
	SEXP ans_xp
);

SEXP ByPos_MIndex_combine(SEXP ends_listlist);


/* match_pattern_at.c */

SEXP debug_match_pattern_at();

int (*_selected_nmismatch_at_Pshift_fun)(
	const RoSeq *P,
	const RoSeq *S,
	int Pshift,
	int max_mm
);

void _select_nmismatch_at_Pshift_fun(
	int fixedP,
	int fixedS
);

int _nedit_for_Ploffset(
	const RoSeq *P,
	const RoSeq *S,
	int Ploffset,
	int max_nedit,
	int loose_Ploffset,
	int *min_width
);

int _nedit_for_Proffset(
	const RoSeq *P,
	const RoSeq *S,
	int Proffset,
	int max_nedit,
	int loose_Proffset,
	int *min_width
);

SEXP match_pattern_at(
	SEXP pattern,
	SEXP subject,
	SEXP at,
	SEXP at_type,
	SEXP max_mismatch,
	SEXP with_indels,
	SEXP fixed,
	SEXP ans_type
);


/* match_pattern_boyermoore.c */

SEXP debug_match_pattern_boyermoore();

void _match_pattern_boyermoore(
	const RoSeq *P,
	const RoSeq *S
);


/* match_pattern_shiftor.c */

SEXP debug_match_pattern_shiftor();

SEXP bits_per_long();

void _match_pattern_shiftor(
	const RoSeq *P,
	const RoSeq *S,
	int max_mm,
	int fixedP,
	int fixedS
);


/* match_pattern_indels.c */

SEXP debug_match_pattern_indels();

void _match_pattern_indels(
	const RoSeq *P,
	const RoSeq *S,
	int max_mm,
	int fixedP,
	int fixedS
);


/* match_pattern.c */

SEXP debug_match_pattern();

SEXP XString_match_pattern(
	SEXP pattern,
	SEXP subject,
	SEXP algorithm,
	SEXP max_mismatch,
	SEXP with_indels,
	SEXP fixed,
	SEXP count_only
);

SEXP XStringViews_match_pattern(
	SEXP pattern,
	SEXP subject,
	SEXP views_start,
	SEXP views_width,
	SEXP algorithm,
	SEXP max_mismatch,
	SEXP with_indels,
	SEXP fixed,
	SEXP count_only
);

SEXP XStringSet_vmatch_pattern(
	SEXP pattern,
	SEXP subject,
	SEXP algorithm,
	SEXP max_mismatch,
	SEXP with_indels,
	SEXP fixed,
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

SEXP PWM_score(
	SEXP pwm,
	SEXP subject,
	SEXP start
);

SEXP match_PWM(
	SEXP pwm,
	SEXP subject,
	SEXP min_score,
	SEXP count_only
);


/* find_palindromes.c */

SEXP debug_find_palindromes();

SEXP find_palindromes(
	SEXP s_xp,
	SEXP s_offset,
	SEXP s_length,
	SEXP min_armlength,
	SEXP max_ngaps,
	SEXP L2R_lkup
);


/* match_pdict_Twobit.c */

SEXP debug_match_pdict_Twobit();

SEXP build_Twobit(
	SEXP tb,
	SEXP dup2unq0,
	SEXP base_codes
);

void _match_Twobit(
	SEXP pptb_comps,
	const RoSeq *S,
	int fixedS
);


/* match_pdict_ACtree.c */

SEXP debug_match_pdict_ACtree();

SEXP free_actree_nodes_buf();

SEXP build_ACtree(
	SEXP tb,
	SEXP dup2unq0,
	SEXP base_codes
);

void _match_ACtree(
	SEXP pptb_comps,
	const RoSeq *S,
	int fixedS
);


/* match_pdict.c */

SEXP debug_match_pdict();

SEXP XString_match_pdict(
	SEXP pptb_comps,
	SEXP pdict_head,
	SEXP pdict_tail,
	SEXP subject,
	SEXP max_mismatch,
	SEXP fixed,
	SEXP count_only,
	SEXP envir
);

SEXP XStringViews_match_pdict(
	SEXP pptb_comps,
	SEXP pdict_head,
	SEXP pdict_tail,
	SEXP subject,
	SEXP views_start,
	SEXP views_width,
	SEXP max_mismatch,
	SEXP fixed,
	SEXP count_only,
	SEXP envir
);


/* align_utils.c */

SEXP PairwiseAlignment_nmatch(
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

SEXP PairwiseAlignment_align_aligned(
	SEXP alignment,
	SEXP gapCode
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

