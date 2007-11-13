#include <Rdefines.h>
#include <R_ext/Rdynload.h>
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

void _Biostrings_reset_views_buffer();

int *_Biostrings_get_views_start();
int *_Biostrings_get_views_end();
char **_Biostrings_get_views_desc();

int _Biostrings_report_view(
		int start,
		int end,
		const char *desc
);

int _Biostrings_report_match(
		int Lpos,
		int Rpos
);

int fgets_rtrimmed(
		char *s,
		int size,
		FILE *stream
);


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

SEXP XInteger_alloc(
		SEXP xint_xp,
		SEXP length
);

SEXP XInteger_show(SEXP xint_xp);

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


/* normalize_views.c */

SEXP Biostrings_normalize_views(
		SEXP start,
		SEXP end
);


/* match_naive.c */

SEXP match_naive_debug();

SEXP match_naive_exact(
		SEXP p_xp,
		SEXP p_offset,
		SEXP p_length,
		SEXP s_xp,
		SEXP s_offset,
		SEXP s_length,
		SEXP count_only
);

SEXP match_naive_fuzzy(
		SEXP p_xp,
		SEXP p_offset,
		SEXP p_length,
		SEXP s_xp,
		SEXP s_offset,
		SEXP s_length,
		SEXP mismatch,
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
		SEXP mismatch,
		SEXP fixed,
		SEXP count_only
);


/* match_BOC.c */

SEXP match_BOC_debug();

SEXP match_BOC_preprocess(
		SEXP s_xp,
		SEXP s_offset,
		SEXP s_length,
		SEXP p_length,
		SEXP code1,
		SEXP buf1_xp,
		SEXP code2,
		SEXP buf2_xp,
		SEXP code3,
		SEXP buf3_xp,
		SEXP code4
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

