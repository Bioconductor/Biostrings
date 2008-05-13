#include "../inst/include/Biostrings_defines.h"
#include <string.h>

#define DEBUG_BIOSTRINGS 1

#define CHRTRTABLE_LENGTH 256

typedef int CharToIntTable[CHRTRTABLE_LENGTH];

#define MAX_CHILDREN_PER_ACNODE 4
typedef struct acnode {
	int parent_id;
	int depth;
	int child_id[MAX_CHILDREN_PER_ACNODE];
	int flink;
	int P_id;
} ACNode;


/* copy_seq.c */

SEXP debug_copy_seq();

void _copy_seq(
	char *dest,
	const char *src,
	size_t n,
	const int *chrtrtable
);

void _revcopy_seq(
	char *dest,
	const char *src,
	size_t n,
	const int *chrtrtable
);

void _copy_seq_from_i1i2(
	int i1, int i2,
	char *dest, int dest_length,
	const char *src, int src_length,
	const int *chrtrtable
);

void _copy_seq_to_i1i2(
	int i1, int i2,
	char *dest, int dest_length,
	const char *src, int src_length,
	const int *chrtrtable
);

void _copy_seq_from_subset(
	const int *subset, int n,
	char *dest, int dest_length,
	const char *src, int src_length,
	const int *chrtrtable
);

void _copy_seq_to_subset(
	const int *subset, int n,
	char *dest, int dest_length,
	const char *src, int src_length,
	const int *chrtrtable
);


/* utils.c */

SEXP debug_utils();

char * _Biostrings_alloc_string(int n);

int _Biostrings_memcmp(
	const char *a,
	int ia,
	const char *b,
	int ib,
	int n,
	size_t size
);

void _Biostrings_memcpy_from_i1i2(
	int i1,
	int i2,
	char *dest,
	size_t dest_nmemb,
	const char *src,
	size_t src_nmemb,
	size_t size
);

void _Biostrings_memcpy_from_subset(
	const int *subset,
	int n,
	char *dest,
	size_t dest_nmemb,
	const char *src,
	size_t src_nmemb,
	size_t size
);

void _Biostrings_memcpy_to_i1i2(
	int i1,
	int i2,
	char *dest,
	size_t dest_nmemb,
	const char *src,
	size_t src_nmemb,
	size_t size
);

void _Biostrings_memcpy_to_subset(
	const int *subset,
	int n,
	char *dest,
	size_t dest_nmemb,
	const char *src,
	size_t src_nmemb,
	size_t size
);

void _Biostrings_translate_charcpy_from_i1i2(
	int i1,
	int i2,
	char *dest,
	int dest_length,
	const char *src,
	int src_length,
	const int *lkup,
	int lkup_length
);

void _Biostrings_translate_charcpy_from_subset(
	const int *subset,
	int n,
	char *dest,
	int dest_length,
	const char *src,
	int src_length,
	const int *lkup,
	int lkup_length
);

void _Biostrings_translate_charcpy_to_i1i2(
	int i1,
	int i2,
	char *dest,
	int dest_length,
	const char *src,
	int src_length,
	const int *lkup,
	int lkup_length
);

void _Biostrings_translate_charcpy_to_subset(
	const int *subset,
	int n,
	char *dest,
	int dest_length,
	const char *src,
	int src_length,
	const int *lkup,
	int lkup_length
);

void _Biostrings_reverse_memcpy_from_i1i2(
	int i1,
	int i2,
	char *dest,
	size_t dest_nmemb,
	const char *src,
	size_t src_nmemb,
	size_t size
);

void _Biostrings_reverse_translate_charcpy_from_i1i2(
	int i1,
	int i2,
	char *dest,
	int dest_length,
	const char *src,
	int src_length,
	const int *lkup,
	int lkup_length
);

void _Biostrings_coerce_to_complex_from_i1i2(
	int i1,
	int i2,
	Rcomplex *dest,
	int dest_length,
	const char *src,
	int src_length,
	const Rcomplex *lkup,
	int lkup_length
);

int fgets_rtrimmed(
	char *s,
	int size,
	FILE *stream
);

void get_intorder(
	int len,
	const int *in,
	int *out
);

void _init_chrtrtable(
	const int *codes,
	int len,
	int *chrtrtable
);


/* bufutils.c */

SEXP debug_bufutils();

IntBuf _new_IntBuf(
	int buflength,
	int nelt
);

void _IntBuf_insert_at(
	IntBuf *ibuf,
	int at,
	int val
);

void _IntBuf_delete_at(
	IntBuf *ibuf,
	int at
);

SEXP _IntBuf_asINTEGER(IntBuf *ibuf);

IntBuf _INTEGER_asIntBuf(SEXP x);

IntBuf _CHARACTER_asIntBuf(
	SEXP x,
	int keyshift
);

IntBBuf _new_IntBBuf(
	int buflength,
	int nelt
);

void _IntBBuf_insert_at(
	IntBBuf *ibbuf,
	int at,
	IntBuf ibuf
);

SEXP _IntBBuf_asLIST(
	IntBBuf *ibbuf,
	int mode
);

IntBBuf _LIST_asIntBBuf(SEXP x);

SEXP _IntBBuf_toEnvir(
	IntBBuf *ibbuf,
	SEXP envir,
	int keyshift
);

RangeBuf _new_RangeBuf(
	int buflength,
	int nelt
);

void _RangeBuf_insert_at(
	RangeBuf *rangebuf,
	int at,
	int start,
	int width
);

CharBuf _new_CharBuf(int buflength);

CharBuf _new_CharBuf_from_string(const char *string);

void _CharBuf_insert_at(
	CharBuf *cbuf,
	int at,
	char c
);

SEXP _CharBuf_asRAW(CharBuf *cbuf);

CharBBuf _new_CharBBuf(
	int buflength,
	int nelt
);

void _CharBBuf_insert_at(
	CharBBuf *cbbuf,
	int at,
	CharBuf cbuf
);

void _append_string_to_CharBBuf(
	CharBBuf *cbbuf,
	const char *string
);


/* IRanges_class.c */

SEXP debug_IRanges_class();

SEXP _get_IRanges_start(SEXP x);

SEXP _get_IRanges_width(SEXP x);

int _get_IRanges_length(SEXP x);

const int *_get_IRanges_start0(SEXP x);

const int *_get_IRanges_width0(SEXP x);

void _set_IRanges_names(
	SEXP x,
	SEXP names
);

void _copy_IRanges_slots(
	SEXP x,
	SEXP x0
);

SEXP _new_IRanges(
	const char *class,
	SEXP start,
	SEXP width,
	SEXP names
);

SEXP _new_IRanges_from_RoSeqs(
	const char *class,
	RoSeqs seqs
);

SEXP _alloc_IRanges(
	const char *class,
	int length
);


/* IRanges_utils.c */

SEXP debug_IRanges_utils();

SEXP narrow_IRanges(
	SEXP x,
	SEXP start,
	SEXP end,
	SEXP width
);

SEXP int_to_adjacent_ranges(SEXP x);

SEXP reduce_IRanges(
	SEXP x,
	SEXP with_inframe_start
);


/* XRaw_class.c */

SEXP debug_XRaw_class();

const char *_get_class(SEXP x);

SEXP Biostrings_sexp_address(SEXP s);

SEXP Biostrings_safe_explode(SEXP s);

SEXP Biostrings_xp_show(SEXP xp);

SEXP Biostrings_xp_new();

SEXP Biostrings_XRaw_alloc(
	SEXP xraw_xp,
	SEXP length
);

SEXP _get_XRaw_tag(SEXP x);

int _get_XRaw_length(SEXP x);

SEXP Biostrings_XRaw_get_show_string(SEXP xraw_xp);

SEXP Biostrings_XRaw_length(SEXP xraw_xp);

SEXP _new_XRaw(SEXP tag);

SEXP _new_XRaw_from_RoSeqs(
	RoSeqs seqs,
	SEXP lkup
);

void _write_RoSeq_to_XRaw(
	SEXP x,
	int offset,
	RoSeq seq,
	const int *chrtrtable
);

SEXP _new_CHARSXP_from_RoSeq(
	RoSeq seq,
	SEXP lkup
);

SEXP _new_STRSXP_from_RoSeqs(
	RoSeqs seqs,
	SEXP lkup
);


/* XRaw_utils.c */

SEXP debug_XRaw_utils();

SEXP Biostrings_XRaw_memcmp(
	SEXP xraw1_xp,
	SEXP start1,
	SEXP xraw2_xp,
	SEXP start2,
	SEXP width
);

SEXP Biostrings_XRaw_copy_from_i1i2(
	SEXP dest_xraw_xp,
	SEXP src_xraw_xp,
	SEXP imin,
	SEXP imax
);

SEXP Biostrings_XRaw_copy_from_subset(
	SEXP dest_xraw_xp,
	SEXP src_xraw_xp,
	SEXP subset
);

SEXP Biostrings_XRaw_read_chars_from_i1i2(
	SEXP src_xraw_xp,
	SEXP imin,
	SEXP imax
);

SEXP Biostrings_XRaw_read_chars_from_subset(
	SEXP src_xraw_xp,
	SEXP subset
);

SEXP Biostrings_XRaw_write_chars_to_i1i2(
	SEXP dest_xraw_xp,
	SEXP imin,
	SEXP imax,
	SEXP string
);

SEXP Biostrings_XRaw_write_chars_to_subset(
	SEXP dest_xraw_xp,
	SEXP subset,
	SEXP string
);

SEXP XRaw_read_ints_from_i1i2(
	SEXP src_xraw_xp,
	SEXP imin,
	SEXP imax
);

SEXP XRaw_read_ints_from_subset(
	SEXP src_xraw_xp,
	SEXP subset
);

SEXP XRaw_write_ints_to_i1i2(
	SEXP dest_xraw_xp,
	SEXP imin,
	SEXP imax,
	SEXP val
);

SEXP XRaw_write_ints_to_subset(
	SEXP dest_xraw_xp,
	SEXP subset,
	SEXP val
);

SEXP XRaw_read_enc_chars_from_i1i2(
	SEXP src_xraw_xp,
	SEXP imin,
	SEXP imax,
	SEXP lkup
);

SEXP XRaw_read_enc_chars_from_subset(
	SEXP src_xraw_xp,
	SEXP subset,
	SEXP lkup
);

SEXP XRaw_write_enc_chars_to_i1i2(
	SEXP dest_xraw_xp,
	SEXP imin,
	SEXP imax,
	SEXP string,
	SEXP lkup
);

SEXP XRaw_write_enc_chars_to_subset(
	SEXP dest_xraw_xp,
	SEXP subset,
	SEXP string,
	SEXP lkup
);

SEXP XRaw_read_complexes_from_i1i2(
	SEXP src_xraw_xp,
	SEXP imin,
	SEXP imax,
	SEXP lkup
);

SEXP XRaw_read_complexes_from_subset(
	SEXP src_xraw_xp,
	SEXP subset,
	SEXP lkup
);

SEXP XRaw_loadFASTA(
	SEXP xraw_xp,
	SEXP filepath,
	SEXP collapse,
	SEXP lkup
);


/* XInteger.c */

SEXP debug_XInteger();

SEXP XInteger_alloc(
	SEXP xint_xp,
	SEXP length
);

SEXP XInteger_get_show_string(SEXP xint_xp);

SEXP XInteger_length(SEXP xint_xp);

SEXP XInteger_memcmp(
	SEXP xint1_xp,
	SEXP start1,
	SEXP xint2_xp,
	SEXP start2,
	SEXP width
);

SEXP XInteger_read_ints_from_i1i2(
	SEXP src_xint_xp,
	SEXP imin,
	SEXP imax
);

SEXP XInteger_read_ints_from_subset(
	SEXP src_xint_xp,
	SEXP subset
);

SEXP XInteger_write_ints_to_i1i2(
	SEXP dest_xint_xp,
	SEXP imin,
	SEXP imax,
	SEXP val
);

SEXP XInteger_write_ints_to_subset(
	SEXP dest_xint_xp,
	SEXP subset,
	SEXP val
);


/* XNumeric.c */

SEXP debug_XNumeric();

SEXP XNumeric_alloc(
	SEXP xnum_xp,
	SEXP length
);

SEXP XNumeric_get_show_string(SEXP xnum_xp);

SEXP XNumeric_length(SEXP xnum_xp);

SEXP XNumeric_memcmp(
	SEXP xnum1_xp,
	SEXP start1,
	SEXP xnum2_xp,
	SEXP start2,
	SEXP width
);

SEXP XNumeric_read_nums_from_i1i2(
	SEXP src_xnum_xp,
	SEXP imin,
	SEXP imax
);

SEXP XNumeric_read_nums_from_subset(
	SEXP src_xnum_xp,
	SEXP subset
);

SEXP XNumeric_write_nums_to_i1i2(
	SEXP dest_xnum_xp,
	SEXP imin,
	SEXP imax,
	SEXP val
);

SEXP XNumeric_write_nums_to_subset(
	SEXP dest_xnum_xp,
	SEXP subset,
	SEXP val
);


/* XString_class.c */

SEXP debug_XString_class();

const int *get_enc_chrtrtable(const char *class);

const int *get_dec_chrtrtable(const char *class);

SEXP init_DNAlkups(SEXP enc_lkup, SEXP dec_lkup);

char _DNAencode(char c);

char _DNAdecode(char code);

SEXP init_RNAlkups(SEXP enc_lkup, SEXP dec_lkup);

char _RNAencode(char c);

char _RNAdecode(char code);

SEXP _get_XString_xdata(SEXP x);

RoSeq _get_XString_asRoSeq(SEXP x);

SEXP _new_XString(
	const char *class,
	SEXP xdata,
	int offset,
	int length
);

SEXP _new_XString_from_RoSeqs(
	const char *class,
	RoSeqs seqs
);

SEXP _alloc_XString(
	const char *class,
	int length
);

void _write_RoSeq_to_XString(
	SEXP x,
	int start,
	RoSeq seq,
	int encode
);


/* XStringSet_class.c */

SEXP debug_XStringSet_class();

const char *_get_XStringSet_baseClass(SEXP x);

int _get_XStringSet_length(SEXP x);

CachedXStringSet _new_CachedXStringSet(SEXP x);

RoSeq _get_CachedXStringSet_elt_asRoSeq(
	CachedXStringSet *x,
	int i
);

RoSeq _get_XStringSet_elt_asRoSeq(
	SEXP x,
	int i
);

SEXP _new_XStringSet_from_RoSeqs(
	const char *baseClass,
	RoSeqs seqs
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
	RoSeq seq,
	int encode
);

void _write_RoSeq_to_XStringSet_elt(
	SEXP x,
	int i,
	RoSeq seq,
	int encode
);

SEXP XStringSet_as_STRSXP(
	SEXP x,
	SEXP lkup
);

int _get_XStringList_length(SEXP x);

RoSeq _get_XStringList_elt_asRoSeq(
	SEXP x,
	int i
);

SEXP XStrings_to_nchars(SEXP x_seqs);


/* seqs_to_seqs.c */

SEXP debug_seqs_to_seqs();

void narrow_RoSeqs(
	RoSeqs *seqs,
	const int *safe_starts,
	const int *safe_widths
);

RoSeqs _new_RoSeqs_from_BBuf(CharBBuf cbbuf);

RoSeqs _new_RoSeqs_from_STRSXP(
	int nseq,
	SEXP x
);

RoSeqs _new_RoSeqs_from_XString(
	int nseq,
	SEXP x
);

RoSeqs _new_RoSeqs_from_XStringSet(
	int nseq,
	SEXP x
);

RoSeqs _new_RoSeqs_from_XStringList(
	int nseq,
	SEXP x
);

SEXP copy_subXRaw(
	SEXP x,
	SEXP start,
	SEXP nchar,
	SEXP lkup
);

SEXP new_XRaw_from_STRSXP(
	SEXP x,
	SEXP safe_starts,
	SEXP safe_widths,
	SEXP collapse,
	SEXP lkup
);

SEXP new_XRaw_from_XString(
	SEXP x,
	SEXP safe_starts,
	SEXP safe_widths,
	SEXP lkup
);

SEXP new_XStringList_from_XRaw(
	SEXP x,
	SEXP safe_starts,
	SEXP safe_widths,
	SEXP proto
);

SEXP narrow_XStringList(
	SEXP x,
	SEXP safe_starts,
	SEXP safe_widths,
	SEXP proto
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

SEXP oligonucleotide_frequency(
	SEXP x_xp,
	SEXP x_offset,
	SEXP x_length,
	SEXP base_codes,
	SEXP width,
	SEXP fast_moving_side
);


/* char_translate.c */

SEXP XRaw_translate_copy_from_i1i2(
	SEXP dest_xraw_xp,
	SEXP src_xraw_xp,
	SEXP imin,
	SEXP imax,
	SEXP lkup
);

SEXP XRaw_translate_copy_from_subset(
	SEXP dest_xraw_xp,
	SEXP src_xraw_xp,
	SEXP subset,
	SEXP lkup
);

SEXP XRaw_reverse_copy_from_i1i2(
	SEXP dest_xraw_xp,
	SEXP src_xraw_xp,
	SEXP imin,
	SEXP imax
);

SEXP XRaw_reverse_translate_copy_from_i1i2(
	SEXP dest_xraw_xp,
	SEXP src_xraw_xp,
	SEXP imin,
	SEXP imax,
	SEXP lkup
);

SEXP XStringSet_char_translate(
	SEXP x,
	SEXP lkup,
	SEXP reverse
);


/* replace_locs.c */

SEXP XString_replace_locs_bySTRSXP(
	SEXP x,
	SEXP loc,
	SEXP letter,
	SEXP lkup,
	SEXP if_not_extending,
	SEXP verbose
);

SEXP XString_inplace_replace_locs_bySTRSXP(
	SEXP x,
	SEXP loc,
	SEXP letter,
	SEXP lkup
);


/* views_buffer.c */

SEXP debug_views_buffer();

void _Biostrings_reset_viewsbuf(int mrmode);

void _init_match_reporting(int mrmode);

void _set_match_shift(int shift);

int _Biostrings_append_view(
	int start,
	int end,
	const char *desc
);

int _Biostrings_report_match(
	int Lpos,
	int Rpos
);

int _report_match(
	int start,
	int end
);

SEXP _Biostrings_viewsbuf_count_asINTEGER();

SEXP _Biostrings_viewsbuf_start_asINTEGER();

SEXP _Biostrings_viewsbuf_end_asINTEGER();

SEXP _Biostrings_viewsbuf_desc_asCHARACTER();

SEXP _Biostrings_viewsbuf_asLIST();

SEXP _reported_matches_asSEXP();


/* normalize_views.c */

SEXP Biostrings_normalize_views(
	SEXP start,
	SEXP end
);


/* match_pattern.c */

SEXP debug_match_pattern();

int _is_matching(
	RoSeq P,
	RoSeq S,
	int Pshift,
	int max_mm,
	int fixedP,
	int fixedS
);

SEXP is_matching(
	SEXP pattern,
	SEXP subject,
	SEXP start,
	SEXP max_mismatch,
	SEXP fixed
);

SEXP XString_match_pattern(
	SEXP pattern,
	SEXP subject,
	SEXP algorithm,
	SEXP max_mismatch,
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
	SEXP fixed,
	SEXP count_only
);


/* match_pattern_boyermoore.c */

SEXP debug_match_pattern_boyermoore();

void _match_pattern_boyermoore(
	RoSeq P,
	RoSeq S
);


/* match_pattern_shiftor.c */

SEXP debug_match_pattern_shiftor();

SEXP bits_per_long();

void _match_pattern_shiftor(
	RoSeq P,
	RoSeq S,
	int max_mm,
	int fixedP,
	int fixedS
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


/* ACtree_utils.c */

SEXP debug_ACtree_utils();

SEXP CWdna_free_actree_nodes_buf();

SEXP CWdna_pp_STRSXP(
	SEXP dict,
	SEXP start,
	SEXP end
);

SEXP CWdna_pp_XStringSet(
	SEXP dict,
	SEXP start,
	SEXP end
);


/* MIndex_utils.c */

SEXP debug_MIndex_utils();

SEXP shiftListOfInts(
	SEXP x,
	SEXP shift
);

SEXP extract_endIndex(
	SEXP ends_envir,
	SEXP shift,
	SEXP names,
	SEXP all_names
);

SEXP ByPos_MIndex_coverage(
	SEXP ends_list,
	SEXP mindex_width,
	SEXP start,
	SEXP end
);

SEXP ByName_MIndex_coverage(
	SEXP ends_envir,
	SEXP mindex_width,
	SEXP start,
	SEXP end
);


/* match_TBdna.c */

SEXP debug_match_TBdna();

SEXP match_TBdna(
	SEXP actree_nodes_xp,
	SEXP actree_base_codes,
	SEXP pdict_dups,
	SEXP pdict_head_XStringSet,
	SEXP pdict_tail_XStringSet,
	SEXP subject_XString,
	SEXP max_mismatch,
	SEXP fixed,
	SEXP count_only,
	SEXP envir
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


/* align_pairwiseAlignment.c */

SEXP align_pairwiseAlignment(
	SEXP string1,
	SEXP string2,
	SEXP quality1,
	SEXP quality2,
	SEXP typeCode,
	SEXP scoreOnly,
	SEXP gapOpening,
	SEXP gapExtension,
	SEXP qualityLookupTable,
	SEXP qualityMatchMatrix,
	SEXP qualityMismatchMatrix,
	SEXP qualityMatrixDim,
	SEXP constantLookupTable,
	SEXP constantMatrix,
	SEXP constantMatrixDim
);
