#include "Biostrings.h"

static const R_CallMethodDef callMethods[] = {

/* utils.c */
	{"utils_debug", (DL_FUNC) &utils_debug, 0},

/* ByteBuffer.c */
	{"ByteBuffer_debug", (DL_FUNC) &ByteBuffer_debug, 0},

	{"sexp_address", (DL_FUNC) &sexp_address, 1},
	{"xp_show", (DL_FUNC) &xp_show, 1},
	{"xp_new", (DL_FUNC) &xp_new, 0},
	{"safe_explode", (DL_FUNC) &safe_explode, 1},

	{"ByteBuffer_alloc", (DL_FUNC) &ByteBuffer_alloc, 2},
	{"ByteBuffer_show", (DL_FUNC) &ByteBuffer_show, 1},
	{"ByteBuffer_length", (DL_FUNC) &ByteBuffer_length, 1},
	{"ByteBuffer_memcmp", (DL_FUNC) &ByteBuffer_memcmp, 5},

	{"ByteBuffer_copy_from_i1i2", (DL_FUNC) &ByteBuffer_copy_from_i1i2, 4},
	{"ByteBuffer_copy_from_subset", (DL_FUNC) &ByteBuffer_copy_from_subset, 3},

	{"ByteBuffer_read_chars_from_i1i2", (DL_FUNC) &ByteBuffer_read_chars_from_i1i2, 3},
	{"ByteBuffer_read_chars_from_subset", (DL_FUNC) &ByteBuffer_read_chars_from_subset, 2},
	{"ByteBuffer_write_chars_to_i1i2", (DL_FUNC) &ByteBuffer_write_chars_to_i1i2, 4},
	{"ByteBuffer_write_chars_to_subset", (DL_FUNC) &ByteBuffer_write_chars_to_subset, 3},

	{"ByteBuffer_read_ints_from_i1i2", (DL_FUNC) &ByteBuffer_read_ints_from_i1i2, 3},
	{"ByteBuffer_read_ints_from_subset", (DL_FUNC) &ByteBuffer_read_ints_from_subset, 2},
	{"ByteBuffer_write_ints_to_i1i2", (DL_FUNC) &ByteBuffer_write_ints_to_i1i2, 4},
	{"ByteBuffer_write_ints_to_subset", (DL_FUNC) &ByteBuffer_write_ints_to_subset, 3},

	{"ByteBuffer_read_enc_chars_from_i1i2", (DL_FUNC) &ByteBuffer_read_enc_chars_from_i1i2, 4},
	{"ByteBuffer_read_enc_chars_from_subset", (DL_FUNC) &ByteBuffer_read_enc_chars_from_subset, 3},
	{"ByteBuffer_write_enc_chars_to_i1i2", (DL_FUNC) &ByteBuffer_write_enc_chars_to_i1i2, 5},
	{"ByteBuffer_write_enc_chars_to_subset", (DL_FUNC) &ByteBuffer_write_enc_chars_to_subset, 4},

/* IntBuffer.c */
	{"IntBuffer_alloc", (DL_FUNC) &IntBuffer_alloc, 2},
	{"IntBuffer_show", (DL_FUNC) &IntBuffer_show, 1},
	{"IntBuffer_length", (DL_FUNC) &IntBuffer_length, 1},
	{"IntBuffer_memcmp", (DL_FUNC) &IntBuffer_memcmp, 5},

	{"IntBuffer_read_ints_from_i1i2", (DL_FUNC) &IntBuffer_read_ints_from_i1i2, 3},
	{"IntBuffer_read_ints_from_subset", (DL_FUNC) &IntBuffer_read_ints_from_subset, 2},
	{"IntBuffer_write_ints_to_i1i2", (DL_FUNC) &IntBuffer_write_ints_to_i1i2, 4},
	{"IntBuffer_write_ints_to_subset", (DL_FUNC) &IntBuffer_write_ints_to_subset, 3},

/* reverseComplement.c */
        {"ByteBuffer_translate_copy_from_i1i2", (DL_FUNC) &ByteBuffer_translate_copy_from_i1i2, 5},
        {"ByteBuffer_translate_copy_from_subset", (DL_FUNC) &ByteBuffer_translate_copy_from_subset, 4},
        {"ByteBuffer_reverse_copy_from_i1i2", (DL_FUNC) &ByteBuffer_reverse_copy_from_i1i2, 4},
        {"ByteBuffer_reverse_translate_copy_from_i1i2", (DL_FUNC) &ByteBuffer_reverse_translate_copy_from_i1i2, 5},

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
