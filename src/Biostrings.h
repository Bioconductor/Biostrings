#include <Rdefines.h>
#include <R_ext/Rdynload.h>

#define DEBUG_BIOSTRINGS 1

/* utils.c */

int _memcmp(char *a, int ia, char *b, int ib, int n, size_t size);
void _memcpy_from_range(int i1, int i2,
	       char *dest, size_t dest_nmemb,
	       char *src, size_t src_nmemb, size_t size);
void _memcpy_to_range(int i1, int i2,
	       char *dest, size_t dest_nmemb,
	       char *src, size_t src_nmemb, size_t size);
void _memcpy_from_subset(int *subset, int n,
	       char *dest, size_t dest_nmemb,
	       char *src, size_t src_nmemb, size_t size);
void _memcpy_to_subset(int *subset, int n,
	       char *dest, size_t dest_nmemb,
	       char *src, size_t src_nmemb, size_t size);


/* bbuf.c */

SEXP sexp_address(SEXP s);
SEXP xp_show(SEXP xp);
SEXP xp_new();
SEXP safe_explode(SEXP s);

SEXP bbuf_alloc(SEXP bb_xp, SEXP length);
SEXP bbuf_show(SEXP bb_xp);
SEXP bbuf_length(SEXP bb_xp);
SEXP bbuf_memcmp(SEXP bb1_xp, SEXP first1,
                 SEXP bb2_xp, SEXP first2, SEXP width);

SEXP bbuf_copy(SEXP dest_xp, SEXP imin, SEXP imax, SEXP src_xp);
SEXP bbuf_copyii(SEXP dest_xp, SEXP ii, SEXP src_xp);

SEXP bbuf_read_chars(SEXP bb_xp, SEXP imin, SEXP imax);
SEXP bbuf_readii_chars(SEXP bb_xp, SEXP ii);
SEXP bbuf_write_chars(SEXP bb_xp, SEXP imin, SEXP imax, SEXP val);
SEXP bbuf_writeii_chars(SEXP bb_xp, SEXP ii, SEXP val);

SEXP bbuf_read_ints(SEXP bb_xp, SEXP imin, SEXP imax);
SEXP bbuf_readii_ints(SEXP bb_xp, SEXP ii);
SEXP bbuf_write_ints(SEXP bb_xp, SEXP imin, SEXP imax, SEXP val);
SEXP bbuf_writeii_ints(SEXP bb_xp, SEXP ii, SEXP val);

SEXP bbuf_read_enc_chars(SEXP bb_xp, SEXP imin, SEXP imax, SEXP dec_xp);
SEXP bbuf_readii_enc_chars(SEXP bb_xp, SEXP ii, SEXP dec_xp);
SEXP bbuf_write_enc_chars(SEXP bb_xp, SEXP imin, SEXP imax,
		SEXP val, SEXP enc_xp);
SEXP bbuf_writeii_enc_chars(SEXP bb_xp, SEXP ii, SEXP val, SEXP enc_xp);


/* ibuf.c */

SEXP ibuf_alloc(SEXP ib_xp, SEXP length);
SEXP ibuf_show(SEXP ib_xp);
SEXP ibuf_length(SEXP ib_xp);
SEXP ibuf_memcmp(SEXP ib1_xp, SEXP first1,
                 SEXP ib2_xp, SEXP first2, SEXP width);

SEXP ibuf_read_ints(SEXP ib_xp, SEXP imin, SEXP imax);
SEXP ibuf_readii_ints(SEXP ib_xp, SEXP ii);
SEXP ibuf_write_ints(SEXP ib_xp, SEXP imin, SEXP imax, SEXP val);
SEXP ibuf_writeii_ints(SEXP ib_xp, SEXP ii, SEXP val);


/* shiftor.c */

SEXP shiftor_debug();
SEXP shiftor(SEXP p_xp, SEXP p_pos, SEXP p_length,
             SEXP s_xp, SEXP s_pos, SEXP s_length,
             SEXP mismatch, SEXP fixed, SEXP count_only);
