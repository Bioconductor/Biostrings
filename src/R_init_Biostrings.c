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

	{"bbuf_copy_from_i1i2", (DL_FUNC) &bbuf_copy_from_i1i2, 4},
	{"bbuf_copy_from_subset", (DL_FUNC) &bbuf_copy_from_subset, 3},

	{"bbuf_read_chars_from_i1i2", (DL_FUNC) &bbuf_read_chars_from_i1i2, 3},
	{"bbuf_read_chars_from_subset", (DL_FUNC) &bbuf_read_chars_from_subset, 2},
	{"bbuf_write_chars_to_i1i2", (DL_FUNC) &bbuf_write_chars_to_i1i2, 4},
	{"bbuf_write_chars_to_subset", (DL_FUNC) &bbuf_write_chars_to_subset, 3},

	{"bbuf_read_ints_from_i1i2", (DL_FUNC) &bbuf_read_ints_from_i1i2, 3},
	{"bbuf_read_ints_from_subset", (DL_FUNC) &bbuf_read_ints_from_subset, 2},
	{"bbuf_write_ints_to_i1i2", (DL_FUNC) &bbuf_write_ints_to_i1i2, 4},
	{"bbuf_write_ints_to_subset", (DL_FUNC) &bbuf_write_ints_to_subset, 3},

	{"bbuf_read_enc_chars_from_i1i2", (DL_FUNC) &bbuf_read_enc_chars_from_i1i2, 4},
	{"bbuf_read_enc_chars_from_subset", (DL_FUNC) &bbuf_read_enc_chars_from_subset, 3},
	{"bbuf_write_enc_chars_to_i1i2", (DL_FUNC) &bbuf_write_enc_chars_to_i1i2, 5},
	{"bbuf_write_enc_chars_to_subset", (DL_FUNC) &bbuf_write_enc_chars_to_subset, 4},

/* ibuf.c */
	{"ibuf_alloc", (DL_FUNC) &ibuf_alloc, 2},
	{"ibuf_show", (DL_FUNC) &ibuf_show, 1},
	{"ibuf_length", (DL_FUNC) &ibuf_length, 1},
	{"ibuf_memcmp", (DL_FUNC) &ibuf_memcmp, 5},

	{"ibuf_read_ints_from_i1i2", (DL_FUNC) &ibuf_read_ints_from_i1i2, 3},
	{"ibuf_read_ints_from_subset", (DL_FUNC) &ibuf_read_ints_from_subset, 2},
	{"ibuf_write_ints_to_i1i2", (DL_FUNC) &ibuf_write_ints_to_i1i2, 4},
	{"ibuf_write_ints_to_subset", (DL_FUNC) &ibuf_write_ints_to_subset, 3},

/* reverseComplement.c */
        {"bbuf_translate_copy_from_i1i2", (DL_FUNC) &bbuf_translate_copy_from_i1i2, 5},
        {"bbuf_translate_copy_from_subset", (DL_FUNC) &bbuf_translate_copy_from_subset, 4},
        {"bbuf_reverse_copy_from_i1i2", (DL_FUNC) &bbuf_reverse_copy_from_i1i2, 4},
        {"bbuf_reverse_translate_copy_from_i1i2", (DL_FUNC) &bbuf_reverse_translate_copy_from_i1i2, 5},

/* alphabetFrequency.c */
	{"alphabetFrequency", (DL_FUNC) &alphabetFrequency, 3},

/* shiftor.c */
	{"shiftor_debug", (DL_FUNC) &shiftor_debug, 0},

	{"shiftor", (DL_FUNC) &shiftor, 9},

	{NULL, NULL, 0}
};

void R_init_Biostrings(DllInfo *info)
{
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
