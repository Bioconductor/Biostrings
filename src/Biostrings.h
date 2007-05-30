#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <string.h>

#define DEBUG_BIOSTRINGS 1
#define R_ALLOC_STRING(n)  (char *) R_alloc((long) (n) + 1L, sizeof(char))

/* utils.c */

SEXP utils_debug();

int Biostrings_memcmp(const char *a, int ia, const char *b, int ib, int n, size_t size);

void Biostrings_memcpy_from_i1i2(int i1, int i2,
		char *dest, size_t dest_nmemb,
		const char *src, size_t src_nmemb, size_t size);
void Biostrings_memcpy_from_subset(int *subset, int n,
		char *dest, size_t dest_nmemb,
		const char *src, size_t src_nmemb, size_t size);

void Biostrings_memcpy_to_i1i2(int i1, int i2,
		char *dest, size_t dest_nmemb,
		const char *src, size_t src_nmemb, size_t size);
void Biostrings_memcpy_to_subset(int *subset, int n,
		char *dest, size_t dest_nmemb,
		const char *src, size_t src_nmemb, size_t size);

void Biostrings_translate_charcpy_from_i1i2(int i1, int i2,
		char *dest, int dest_length,
		const char *src, int src_length,
		int *lkup, int lkup_length);
void Biostrings_translate_charcpy_from_subset(int *subset, int n,
		char *dest, int dest_length,
		const char *src, int src_length,
		int *lkup, int lkup_length);

void Biostrings_translate_charcpy_to_i1i2(int i1, int i2,
		char *dest, int dest_length,
		const char *src, int src_length,
		int *lkup, int lkup_length);
void Biostrings_translate_charcpy_to_subset(int *subset, int n,
		char *dest, int dest_length,
		const char *src, int src_length,
		int *lkup, int lkup_length);

void Biostrings_reverse_memcpy_from_i1i2(int i1, int i2,
		char *dest, size_t dest_nmemb,
		const char *src, size_t src_nmemb, size_t size);
void Biostrings_reverse_translate_charcpy_from_i1i2(int i1, int i2,
		char *dest, int dest_length,
		const char *src, int src_length,
		int *lkup, int lkup_length);

void Biostrings_coerce_to_complex_from_i1i2(int i1, int i2,
		Rcomplex *dest, int dest_length,
		const char *src, int src_length,
		Rcomplex *lkup, int lkup_length);

int Biostrings_estimateExpectedMatchCount(int nP, int nS, int nalphabet);
SEXP Biostrings_expandMatchIndex(SEXP index, int ndone, int nleft);
int *Biostrings_resetMatchPosBuffer();
int Biostrings_reportMatch(int pos);


/* XRaw.c */

SEXP XRaw_debug();

SEXP sexp_address(SEXP s);
SEXP xp_show(SEXP xp);
SEXP xp_new();
SEXP safe_explode(SEXP s);

SEXP XRaw_alloc(SEXP xraw_xp, SEXP length);
SEXP XRaw_get_show_string(SEXP xraw_xp);
SEXP XRaw_length(SEXP xraw_xp);
SEXP XRaw_memcmp(SEXP xraw1_xp, SEXP start1,
		SEXP xraw2_xp, SEXP start2, SEXP width);

SEXP XRaw_copy_from_i1i2(SEXP dest_xraw_xp, SEXP src_xraw_xp, SEXP imin, SEXP imax);
SEXP XRaw_copy_from_subset(SEXP dest_xraw_xp, SEXP src_xraw_xp, SEXP subset);

SEXP XRaw_read_chars_from_i1i2(SEXP src_xraw_xp, SEXP imin, SEXP imax);
SEXP XRaw_read_chars_from_subset(SEXP src_xraw_xp, SEXP subset);
SEXP XRaw_write_chars_to_i1i2(SEXP dest_xraw_xp, SEXP imin, SEXP imax, SEXP string);
SEXP XRaw_write_chars_to_subset(SEXP dest_xraw_xp, SEXP subset, SEXP string);

SEXP XRaw_read_ints_from_i1i2(SEXP src_xraw_xp, SEXP imin, SEXP imax);
SEXP XRaw_read_ints_from_subset(SEXP src_xraw_xp, SEXP subset);
SEXP XRaw_write_ints_to_i1i2(SEXP dest_xraw_xp, SEXP imin, SEXP imax, SEXP val);
SEXP XRaw_write_ints_to_subset(SEXP dest_xraw_xp, SEXP subset, SEXP val);

SEXP XRaw_read_enc_chars_from_i1i2(SEXP src_xraw_xp, SEXP imin, SEXP imax, SEXP lkup);
SEXP XRaw_read_enc_chars_from_subset(SEXP src_xraw_xp, SEXP subset, SEXP lkup);
SEXP XRaw_write_enc_chars_to_i1i2(SEXP dest_xraw_xp, SEXP imin, SEXP imax,
		SEXP string, SEXP lkup);
SEXP XRaw_write_enc_chars_to_subset(SEXP dest_xraw_xp, SEXP subset,
		SEXP string, SEXP lkup);

SEXP XRaw_read_complexes_from_i1i2(SEXP src_xraw_xp, SEXP imin, SEXP imax, SEXP lkup);
SEXP XRaw_read_complexes_from_subset(SEXP src_xraw_xp, SEXP subset, SEXP lkup);


/* XInteger.c */

SEXP XInteger_alloc(SEXP xint_xp, SEXP length);
SEXP XInteger_show(SEXP xint_xp);
SEXP XInteger_length(SEXP xint_xp);
SEXP XInteger_memcmp(SEXP xint1_xp, SEXP start1,
		SEXP xint2_xp, SEXP start2, SEXP width);

SEXP XInteger_read_ints_from_i1i2(SEXP src_xraw_xp, SEXP imin, SEXP imax);
SEXP XInteger_read_ints_from_subset(SEXP src_xraw_xp, SEXP subset);
SEXP XInteger_write_ints_to_i1i2(SEXP dest_xraw_xp, SEXP imin, SEXP imax, SEXP val);
SEXP XInteger_write_ints_to_subset(SEXP dest_xraw_xp, SEXP subset, SEXP val);


/* reverseComplement.c */

SEXP XRaw_translate_copy_from_i1i2(SEXP dest_xraw_xp, SEXP src_xraw_xp, SEXP imin, SEXP imax, SEXP lkup);
SEXP XRaw_translate_copy_from_subset(SEXP dest_xraw_xp, SEXP src_xraw_xp, SEXP subset, SEXP lkup);
SEXP XRaw_reverse_copy_from_i1i2(SEXP dest_xraw_xp, SEXP src_xraw_xp, SEXP imin, SEXP imax);
SEXP XRaw_reverse_translate_copy_from_i1i2(SEXP dest_xraw_xp, SEXP src_xraw_xp, SEXP imin, SEXP imax, SEXP lkup);


/* char_frequency.c */

SEXP char_frequency(SEXP x_xp, SEXP x_offset, SEXP x_length);


/* match_naive.c */

SEXP match_naive_debug();

SEXP match_naive(SEXP p_xp, SEXP p_offset, SEXP p_length,
		SEXP s_xp, SEXP s_offset, SEXP s_length,
		SEXP count_only);


/* match_boyermoore.c */

SEXP match_boyermoore_debug();

SEXP match_boyermoore(SEXP p_xp, SEXP p_offset, SEXP p_length,
		SEXP s_xp, SEXP s_offset, SEXP s_length,
		SEXP count_only);


/* match_shiftor.c */

SEXP match_shiftor_debug();
SEXP bits_per_long();

SEXP match_shiftor(SEXP p_xp, SEXP p_offset, SEXP p_length,
		SEXP s_xp, SEXP s_offset, SEXP s_length,
		SEXP mismatch, SEXP fixed, SEXP count_only);

