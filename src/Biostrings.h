#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <string.h>

#define DEBUG_BIOSTRINGS 1


/* utils.c */

SEXP utils_debug();

int Biostrings_memcmp(char *a, int ia, char *b, int ib, int n, size_t size);

void Biostrings_memcpy_from_i1i2(int i1, int i2,
		char *dest, size_t dest_nmemb,
		char *src, size_t src_nmemb, size_t size);
void Biostrings_memcpy_from_subset(int *subset, int n,
		char *dest, size_t dest_nmemb,
		char *src, size_t src_nmemb, size_t size);

void Biostrings_memcpy_to_i1i2(int i1, int i2,
		char *dest, size_t dest_nmemb,
		char *src, size_t src_nmemb, size_t size);
void Biostrings_memcpy_to_subset(int *subset, int n,
		char *dest, size_t dest_nmemb,
		char *src, size_t src_nmemb, size_t size);

void Biostrings_translate_charcpy_from_i1i2(int i1, int i2,
		char *dest, int dest_length,
		char *src, int src_length,
		int *lkup, int lkup_length);
void Biostrings_translate_charcpy_from_subset(int *subset, int n,
		char *dest, int dest_length,
		char *src, int src_length,
		int *lkup, int lkup_length);

void Biostrings_translate_charcpy_to_i1i2(int i1, int i2,
		char *dest, int dest_length,
		char *src, int src_length,
		int *lkup, int lkup_length);
void Biostrings_translate_charcpy_to_subset(int *subset, int n,
		char *dest, int dest_length,
		char *src, int src_length,
		int *lkup, int lkup_length);

void Biostrings_reverse_memcpy_from_i1i2(int i1, int i2,
		char *dest, size_t dest_nmemb,
		char *src, size_t src_nmemb, size_t size);
void Biostrings_reverse_translate_charcpy_from_i1i2(int i1, int i2,
		char *dest, int dest_length,
		char *src, int src_length,
		int *lkup, int lkup_length);

void Biostrings_coerce_to_complex_from_i1i2(int i1, int i2,
		Rcomplex *dest, int dest_length,
		char *src, int src_length,
		Rcomplex *lkup, int lkup_length);

int Biostrings_estimateExpectedMatchCount(int nP, int nS, int nalphabet);
SEXP Biostrings_expandMatchIndex(SEXP index, int ndone, int nleft);
int *Biostrings_resetMatchPosBuffer();
int Biostrings_reportMatch(int pos);


/* CharBuffer.c */

SEXP CharBuffer_debug();

SEXP sexp_address(SEXP s);
SEXP xp_show(SEXP xp);
SEXP xp_new();
SEXP safe_explode(SEXP s);

SEXP CharBuffer_alloc(SEXP cb_xp, SEXP length);
SEXP CharBuffer_get_show_string(SEXP cb_xp);
SEXP CharBuffer_length(SEXP cb_xp);
SEXP CharBuffer_memcmp(SEXP cb1_xp, SEXP start1,
		SEXP cb2_xp, SEXP start2, SEXP width);

SEXP CharBuffer_copy_from_i1i2(SEXP dest_xp, SEXP src_xp, SEXP imin, SEXP imax);
SEXP CharBuffer_copy_from_subset(SEXP dest_xp, SEXP src_xp, SEXP subset);

SEXP CharBuffer_read_chars_from_i1i2(SEXP src_xp, SEXP imin, SEXP imax);
SEXP CharBuffer_read_chars_from_subset(SEXP src_xp, SEXP subset);
SEXP CharBuffer_write_chars_to_i1i2(SEXP dest_xp, SEXP imin, SEXP imax, SEXP string);
SEXP CharBuffer_write_chars_to_subset(SEXP dest_xp, SEXP subset, SEXP string);

SEXP CharBuffer_read_ints_from_i1i2(SEXP src_xp, SEXP imin, SEXP imax);
SEXP CharBuffer_read_ints_from_subset(SEXP src_xp, SEXP subset);
SEXP CharBuffer_write_ints_to_i1i2(SEXP dest_xp, SEXP imin, SEXP imax, SEXP val);
SEXP CharBuffer_write_ints_to_subset(SEXP dest_xp, SEXP subset, SEXP val);

SEXP CharBuffer_read_enc_chars_from_i1i2(SEXP src_xp, SEXP imin, SEXP imax, SEXP lkup);
SEXP CharBuffer_read_enc_chars_from_subset(SEXP src_xp, SEXP subset, SEXP lkup);
SEXP CharBuffer_write_enc_chars_to_i1i2(SEXP dest_xp, SEXP imin, SEXP imax,
		SEXP string, SEXP lkup);
SEXP CharBuffer_write_enc_chars_to_subset(SEXP dest_xp, SEXP subset,
		SEXP string, SEXP lkup);

SEXP CharBuffer_read_complexes_from_i1i2(SEXP src_xp, SEXP imin, SEXP imax, SEXP lkup);
SEXP CharBuffer_read_complexes_from_subset(SEXP src_xp, SEXP subset, SEXP lkup);


/* IntBuffer.c */

SEXP IntBuffer_alloc(SEXP ib_xp, SEXP length);
SEXP IntBuffer_show(SEXP ib_xp);
SEXP IntBuffer_length(SEXP ib_xp);
SEXP IntBuffer_memcmp(SEXP ib1_xp, SEXP start1,
		SEXP ib2_xp, SEXP start2, SEXP width);

SEXP IntBuffer_read_ints_from_i1i2(SEXP src_xp, SEXP imin, SEXP imax);
SEXP IntBuffer_read_ints_from_subset(SEXP src_xp, SEXP subset);
SEXP IntBuffer_write_ints_to_i1i2(SEXP dest_xp, SEXP imin, SEXP imax, SEXP val);
SEXP IntBuffer_write_ints_to_subset(SEXP dest_xp, SEXP subset, SEXP val);


/* reverseComplement.c */

SEXP CharBuffer_translate_copy_from_i1i2(SEXP dest_xp, SEXP src_xp, SEXP imin, SEXP imax, SEXP lkup);
SEXP CharBuffer_translate_copy_from_subset(SEXP dest_xp, SEXP src_xp, SEXP subset, SEXP lkup);
SEXP CharBuffer_reverse_copy_from_i1i2(SEXP dest_xp, SEXP src_xp, SEXP imin, SEXP imax);
SEXP CharBuffer_reverse_translate_copy_from_i1i2(SEXP dest_xp, SEXP src_xp, SEXP imin, SEXP imax, SEXP lkup);


/* char_frequency.c */

SEXP CharBuffer_char_frequency(SEXP x_xp, SEXP x_offset, SEXP x_length);


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

