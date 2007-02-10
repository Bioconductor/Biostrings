#include "Biostrings.h"

static const R_CallMethodDef callMethods[] = {

/* utils.c */
	{"utils_debug", (DL_FUNC) &utils_debug, 0},

/* CharBuffer.c */
	{"CharBuffer_debug", (DL_FUNC) &CharBuffer_debug, 0},

	{"sexp_address", (DL_FUNC) &sexp_address, 1},
	{"xp_show", (DL_FUNC) &xp_show, 1},
	{"xp_new", (DL_FUNC) &xp_new, 0},
	{"safe_explode", (DL_FUNC) &safe_explode, 1},

	{"CharBuffer_alloc", (DL_FUNC) &CharBuffer_alloc, 2},
	{"CharBuffer_get_show_string", (DL_FUNC) &CharBuffer_get_show_string, 1},
	{"CharBuffer_length", (DL_FUNC) &CharBuffer_length, 1},
	{"CharBuffer_memcmp", (DL_FUNC) &CharBuffer_memcmp, 5},

	{"CharBuffer_copy_from_i1i2", (DL_FUNC) &CharBuffer_copy_from_i1i2, 4},
	{"CharBuffer_copy_from_subset", (DL_FUNC) &CharBuffer_copy_from_subset, 3},

	{"CharBuffer_read_chars_from_i1i2", (DL_FUNC) &CharBuffer_read_chars_from_i1i2, 3},
	{"CharBuffer_read_chars_from_subset", (DL_FUNC) &CharBuffer_read_chars_from_subset, 2},
	{"CharBuffer_write_chars_to_i1i2", (DL_FUNC) &CharBuffer_write_chars_to_i1i2, 4},
	{"CharBuffer_write_chars_to_subset", (DL_FUNC) &CharBuffer_write_chars_to_subset, 3},

	{"CharBuffer_read_ints_from_i1i2", (DL_FUNC) &CharBuffer_read_ints_from_i1i2, 3},
	{"CharBuffer_read_ints_from_subset", (DL_FUNC) &CharBuffer_read_ints_from_subset, 2},
	{"CharBuffer_write_ints_to_i1i2", (DL_FUNC) &CharBuffer_write_ints_to_i1i2, 4},
	{"CharBuffer_write_ints_to_subset", (DL_FUNC) &CharBuffer_write_ints_to_subset, 3},

	{"CharBuffer_read_enc_chars_from_i1i2", (DL_FUNC) &CharBuffer_read_enc_chars_from_i1i2, 4},
	{"CharBuffer_read_enc_chars_from_subset", (DL_FUNC) &CharBuffer_read_enc_chars_from_subset, 3},
	{"CharBuffer_write_enc_chars_to_i1i2", (DL_FUNC) &CharBuffer_write_enc_chars_to_i1i2, 5},
	{"CharBuffer_write_enc_chars_to_subset", (DL_FUNC) &CharBuffer_write_enc_chars_to_subset, 4},

	{"CharBuffer_read_complexes_from_i1i2", (DL_FUNC) &CharBuffer_read_complexes_from_i1i2, 4},
	{"CharBuffer_read_complexes_from_subset", (DL_FUNC) &CharBuffer_read_complexes_from_subset, 3},

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
        {"CharBuffer_translate_copy_from_i1i2", (DL_FUNC) &CharBuffer_translate_copy_from_i1i2, 5},
        {"CharBuffer_translate_copy_from_subset", (DL_FUNC) &CharBuffer_translate_copy_from_subset, 4},
        {"CharBuffer_reverse_copy_from_i1i2", (DL_FUNC) &CharBuffer_reverse_copy_from_i1i2, 4},
        {"CharBuffer_reverse_translate_copy_from_i1i2", (DL_FUNC) &CharBuffer_reverse_translate_copy_from_i1i2, 5},

/* char_frequency.c */
	{"CharBuffer_char_frequency", (DL_FUNC) &CharBuffer_char_frequency, 3},

/* match_naive.c */
	{"match_naive_debug", (DL_FUNC) &match_naive_debug, 0},
	{"match_naive", (DL_FUNC) &match_naive, 7},

/* match_boyermoore.c */
	{"match_boyermoore_debug", (DL_FUNC) &match_boyermoore_debug, 0},
	{"match_boyermoore", (DL_FUNC) &match_boyermoore, 7},

/* match_shiftor.c */
	{"match_shiftor_debug", (DL_FUNC) &match_shiftor_debug, 0},
	{"bits_per_long", (DL_FUNC) &bits_per_long, 0},
	{"match_shiftor", (DL_FUNC) &match_shiftor, 9},

	{NULL, NULL, 0}
};

void R_init_Biostrings(DllInfo *info)
{
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

