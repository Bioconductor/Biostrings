#include "../inst/include/Biostrings_defines.h"
#include <string.h>

#define DEBUG_BIOSTRINGS 1


/* utils.c */

SEXP Biostrings_debug_utils();

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
	int *subset,
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
	int *subset,
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
	int *subset,
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
	int *subset,
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

void _init_code2offset_lkup(
	const int *codes,
	int len,
	int *lkup
);


/* bufutils.c */

SEXP Biostrings_debug_bufutils();

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


/* IRanges.c */

SEXP Biostrings_debug_IRanges();

int _get_IRanges_length(SEXP x);

const int *_get_IRanges_start(SEXP x);

const int *_get_IRanges_width(SEXP x);

SEXP _new_IRanges(
	SEXP start,
	SEXP width,
	SEXP names
);

SEXP _new_IRanges_from_RoSeqs(RoSeqs seqs);

SEXP _replace_IRanges_names(
	SEXP x,
	SEXP names
);

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


/* XRaw_utils.c */

SEXP Biostrings_debug_XRaw_utils();

const char *get_class(SEXP x);

SEXP Biostrings_sexp_address(SEXP s);

SEXP Biostrings_safe_explode(SEXP s);

SEXP Biostrings_xp_show(SEXP xp);

SEXP Biostrings_xp_new();

SEXP Biostrings_XRaw_alloc(
	SEXP xraw_xp,
	SEXP length
);

int _get_XRaw_length(SEXP x);

SEXP Biostrings_XRaw_get_show_string(SEXP xraw_xp);

SEXP Biostrings_XRaw_length(SEXP xraw_xp);

SEXP _new_XRaw(SEXP tag);

SEXP _new_XRaw_from_RoSeqs(
	RoSeqs seqs,
	SEXP lkup
);

SEXP _new_STRSXP_from_RoSeqs(
	RoSeqs seqs,
	SEXP lkup
);


/* XRaw_fillread.c */

SEXP Biostrings_debug_XRaw_fillread();

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

SEXP Biostrings_debug_XInteger();

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
	SEXP src_xraw_xp,
	SEXP imin,
	SEXP imax
);

SEXP XInteger_read_ints_from_subset(
	SEXP src_xraw_xp,
	SEXP subset
);

SEXP XInteger_write_ints_to_i1i2(
	SEXP dest_xraw_xp,
	SEXP imin,
	SEXP imax,
	SEXP val
);

SEXP XInteger_write_ints_to_subset(
	SEXP dest_xraw_xp,
	SEXP subset,
	SEXP val
);


/* XString_utils.c */

SEXP Biostrings_debug_XString_utils();

SEXP init_DNAlkups(SEXP enc_lkup, SEXP dec_lkup);

char _DNAencode(char c);

char _DNAdecode(char code);

SEXP init_RNAlkups(SEXP enc_lkup, SEXP dec_lkup);

char _RNAencode(char c);

char _RNAdecode(char code);

const char *_get_XString_charseq(
	SEXP x,
	int *length
);

SEXP _new_XString(
	const char *class,
	SEXP data,
	int offset,
	int length
);

SEXP _new_XString_from_RoSeqs(
	const char *class,
	RoSeqs seqs
);

int _get_XStringSet_length(SEXP x);

const char *_get_XStringSet_charseq(
	SEXP x,
	int i,
	int *nchar
);

SEXP _new_XStringSet(
	const char *baseClass,
	RoSeqs seqs
);

void _set_XStringSet_names(
	SEXP x,
	RoSeqs names
);

int _get_XStringList_length(SEXP x);

const char *_get_XStringList_charseq(
	SEXP x,
	int i,
	int *nchar
);

SEXP XStrings_to_nchars(SEXP x_seqs);


/* seqs_to_seqs.c */

SEXP Biostrings_debug_seqs_to_seqs();

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


/* reverseComplement.c */

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


/* char_frequency.c */

SEXP char_frequency(
	SEXP x_xp,
	SEXP x_offset,
	SEXP x_length
);

SEXP oligonucleotide_frequency(
	SEXP x_xp,
	SEXP x_offset,
	SEXP x_length,
	SEXP base_codes,
	SEXP width,
	SEXP fast_moving_side
);


/* views_buffer.c */

SEXP Biostrings_debug_views_buffer();

void _Biostrings_reset_viewsbuf(int mrmode);

void _init_match_reporting(int mrmode);

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


/* match_naive.c */

SEXP match_naive_debug();

int _is_matching(
	const char *P,
	int nP,
	const char *S,
	int nS,
	int Pshift,
	int max_mm,
	int fixedP,
	int fixedS
);

SEXP is_matching(
	SEXP pattern_XString,
	SEXP subject_XString,
	SEXP start,
	SEXP max_mismatch,
	SEXP fixed
);

SEXP match_naive_exact(
	SEXP pattern_XString,
	SEXP subject_XString,
	SEXP count_only
);

SEXP match_naive_inexact(
	SEXP pattern_XString,
	SEXP subject_XString,
	SEXP max_mismatch,
	SEXP fixed,
	SEXP count_only
);


/* match_boyermoore.c */

SEXP match_boyermoore_debug();

SEXP match_boyermoore(
	SEXP p_xp,
	SEXP p_offset,
	SEXP p_length,
	SEXP s_xp,
	SEXP s_offset,
	SEXP s_length,
	SEXP count_only
);


/* match_shiftor.c */

SEXP match_shiftor_debug();

SEXP bits_per_long();

SEXP match_shiftor(
	SEXP p_xp,
	SEXP p_offset,
	SEXP p_length,
	SEXP s_xp,
	SEXP s_offset,
	SEXP s_length,
	SEXP max_mismatch,
	SEXP fixed,
	SEXP count_only
);


/* find_palindromes.c */

SEXP find_palindromes_debug();

SEXP find_palindromes(
	SEXP s_xp,
	SEXP s_offset,
	SEXP s_length,
	SEXP min_armlength,
	SEXP max_ngaps,
	SEXP L2R_lkup
);


/* match_BOC.c */

SEXP match_BOC_debug();

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

SEXP match_BOC2_debug();

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


/* match_TBdna.c */

SEXP match_TBdna_debug();

SEXP CWdna_free_actree_nodes_buf();

SEXP CWdna_pp_charseqs(
	SEXP dict,
	SEXP start,
	SEXP end
);

SEXP CWdna_pp_XStringSet(
	SEXP dict,
	SEXP start,
	SEXP end
);

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

