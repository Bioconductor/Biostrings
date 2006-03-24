#include "Biostrings.h"

static const R_CallMethodDef callMethods[] = {

/* utils.c */
	{"utils_debug", (DL_FUNC) &utils_debug, 0},

/* bbuf.c */
	{"bbuf_debug", (DL_FUNC) &bbuf_debug, 0},
	{"sexp_address", (DL_FUNC) &sexp_address, 1},
	{"xp_show", (DL_FUNC) &xp_show, 1},
	{"xp_new", (DL_FUNC) &xp_new, 0},
	{"safe_explode", (DL_FUNC) &safe_explode, 1},

	{"bbuf_alloc", (DL_FUNC) &bbuf_alloc, 2},
	{"bbuf_show", (DL_FUNC) &bbuf_show, 1},
	{"bbuf_length", (DL_FUNC) &bbuf_length, 1},
	{"bbuf_memcmp", (DL_FUNC) &bbuf_memcmp, 5},

	{"bbuf_copy", (DL_FUNC) &bbuf_copy, 4},
	{"bbuf_copyii", (DL_FUNC) &bbuf_copyii, 3},

	{"bbuf_read_chars", (DL_FUNC) &bbuf_read_chars, 3},
	{"bbuf_readii_chars", (DL_FUNC) &bbuf_readii_chars, 2},
	{"bbuf_write_chars", (DL_FUNC) &bbuf_write_chars, 4},
	{"bbuf_writeii_chars", (DL_FUNC) &bbuf_writeii_chars, 3},

	{"bbuf_read_ints", (DL_FUNC) &bbuf_read_ints, 3},
	{"bbuf_readii_ints", (DL_FUNC) &bbuf_readii_ints, 2},
	{"bbuf_write_ints", (DL_FUNC) &bbuf_write_ints, 4},
	{"bbuf_writeii_ints", (DL_FUNC) &bbuf_writeii_ints, 3},

	{"bbuf_read_enc_chars", (DL_FUNC) &bbuf_read_enc_chars, 4},
	{"bbuf_readii_enc_chars", (DL_FUNC) &bbuf_readii_enc_chars, 3},
	{"bbuf_write_enc_chars", (DL_FUNC) &bbuf_write_enc_chars, 5},
	{"bbuf_writeii_enc_chars", (DL_FUNC) &bbuf_writeii_enc_chars, 4},

/* ibuf.c */
	{"ibuf_alloc", (DL_FUNC) &ibuf_alloc, 2},
	{"ibuf_show", (DL_FUNC) &ibuf_show, 1},
	{"ibuf_length", (DL_FUNC) &ibuf_length, 1},
	{"ibuf_memcmp", (DL_FUNC) &ibuf_memcmp, 5},

	{"ibuf_read_ints", (DL_FUNC) &ibuf_read_ints, 3},
	{"ibuf_readii_ints", (DL_FUNC) &ibuf_readii_ints, 2},
	{"ibuf_write_ints", (DL_FUNC) &ibuf_write_ints, 4},
	{"ibuf_writeii_ints", (DL_FUNC) &ibuf_writeii_ints, 3},

/* shiftor.c */
	{"shiftor_debug", (DL_FUNC) &shiftor_debug, 0},
	{"shiftor", (DL_FUNC) &shiftor, 9},

	{NULL, NULL, 0}
};

void R_init_Biostrings(DllInfo *info)
{
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
