#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <string.h>

#define DEBUG_BIOSTRINGS 1

/*
 * Buffer structures used for temporary storage of arrays (or arrays of arrays)
 * of ints, chars, strings, etc...
 * Only buffers of ints and chars are suported for now.
 */

typedef struct ibuf {
        int *vals;
        int maxcount;
        int count;
} IBuf; // Buffer for an array of integers

typedef struct ibbuf {
        IBuf *ibufs;
        int maxcount;
        int count;
} IBBuf; // Buffer for an array of arrays of integers

typedef struct cbuf {
        char *vals;
        int maxcount;
        int count;
} CBuf; // Buffer for an array of chars


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
		const int *in,
		int *out,
		int len
);

void _init_code2offset_lkup(
		const int *codes,
		int len,
		int *lkup
);


/* bufutils.c */

SEXP Biostrings_debug_bufutils();

void _IBuf_init(
		IBuf *ibuf,
		int maxcount,
		int count
);

void _IBuf_get_more_room(
		IBuf *ibuf
);

void _IBuf_insert_at(
		IBuf *ibuf,
		int at,
		int val
);

void _IBuf_delete_at(
		IBuf *ibuf,
		int at
);

SEXP _IBuf_asINTEGER(
		IBuf *ibuf
);

IBuf _INTEGER_asIBuf(
		SEXP x
);

IBuf _CHARACTER_asIBuf(
		SEXP x,
		int keyshift
);

void _IBBuf_init(
		IBBuf *ibbuf,
		int maxcount,
		int count
);

void _IBBuf_get_more_room(
		IBBuf *ibbuf
);

void _IBBuf_insert_at(
		IBBuf *ibbuf,
		int at,
		IBuf ibuf
);

SEXP _IBBuf_asLIST(
		IBBuf *ibbuf,
		int mode
);

IBBuf _LIST_asIBBuf(
		SEXP x
);

SEXP _IBBuf_toEnvir(
		IBBuf *ibbuf,
		SEXP envir,
		int keyshift
);

void _CBuf_init(
		CBuf *cbuf,
		int maxcount
);

void _CBuf_get_more_room(
		CBuf *cbuf
);

void _CBuf_insert_at(
		CBuf *cbuf,
		int at,
		char val
);

SEXP _CBuf_asRAW(
		CBuf *cbuf
);


/* views_buffer.c */

SEXP Biostrings_debug_views_buffer();

void _Biostrings_reset_viewsbuf(
		int reporting_mode
);

int _Biostrings_append_view(
		int start,
		int end,
		const char *desc
);

int _Biostrings_report_match(
		int Lpos,
		int Rpos
);

SEXP _Biostrings_viewsbuf_count_asINTEGER();

SEXP _Biostrings_viewsbuf_start_asINTEGER();

SEXP _Biostrings_viewsbuf_end_asINTEGER();

SEXP _Biostrings_viewsbuf_desc_asCHARACTER();

SEXP _Biostrings_viewsbuf_asLIST();


/* XRaw.c */

SEXP Biostrings_debug_XRaw();

SEXP Biostrings_sexp_address(SEXP s);

SEXP Biostrings_xp_show(SEXP xp);

SEXP Biostrings_xp_new();

SEXP Biostrings_safe_explode(SEXP s);

SEXP Biostrings_XRaw_alloc(
		SEXP xraw_xp,
		SEXP length
);

SEXP Biostrings_XRaw_get_show_string(SEXP xraw_xp);

SEXP Biostrings_XRaw_length(SEXP xraw_xp);

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


/* BString_utils.c */

const char *get_BString_seq(
		SEXP bstring,
		int *seq_len
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
		SEXP pattern_BString,
		SEXP subject_BString,
		SEXP start,
		SEXP max_mismatch,
		SEXP fixed
);

SEXP match_naive_exact(
		SEXP pattern_BString,
		SEXP subject_BString,
		SEXP count_only
);

SEXP match_naive_inexact(
		SEXP pattern_BString,
		SEXP subject_BString,
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


/* match_TPdna.c */

SEXP match_TPdna_debug();

SEXP ULdna_free_actree_nodes_buf();

SEXP ULdna_pp_StrVect(
		SEXP dict,
		SEXP width
);

SEXP ULdna_pp_BStringList(
		SEXP dict,
		SEXP width
);

SEXP ULdna_pp_views(
		SEXP dict_subj_xp,
		SEXP dict_subj_offset,
		SEXP dict_subj_length,
		SEXP dict_start,
		SEXP dict_end,
		SEXP width
);

SEXP match_TPdna(
		SEXP actree_nodes_xp,
		SEXP actree_base_codes,
		SEXP pdict_dups,
		SEXP pdict_tail_bstrings,
		SEXP subject_BString,
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
		SEXP s1_xp,
		SEXP s1_offset,
		SEXP s1_length,
		SEXP s2_xp,
		SEXP s2_offset,
		SEXP s2_length,
		SEXP mat,
		SEXP mat_nrow,
		SEXP lkup,
		SEXP gap_cost,
		SEXP gap_code
);

