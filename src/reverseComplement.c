#include "Biostrings.h"


/*
 * See CharBuffer_copy_from_i1i2() in CharBuffer.c for a description of the first 4 arguments.
 */
SEXP CharBuffer_translate_copy_from_i1i2(SEXP dest_xp, SEXP src_xp, SEXP imin, SEXP imax, SEXP hash_xp)
{
	SEXP dest, src, hash;
	int i1, i2;
	char hash_hole;

	dest = R_ExternalPtrTag(dest_xp);
	i1 = INTEGER(imin)[0] - 1;
	i2 = INTEGER(imax)[0] - 1;
	src = R_ExternalPtrTag(src_xp);
	/* hash = STRING_ELT(hash_xp, 0); */
	hash = R_ExternalPtrTag(hash_xp);
	hash_hole = CHAR(hash)[0];
        Biostrings_translate_charcpy_from_i1i2(i1, i2,
		CHAR(dest), LENGTH(dest), CHAR(src), LENGTH(src),
                CHAR(hash) + 1, LENGTH(hash) - 1, hash_hole, 0);
	return dest_xp;
}

/*
 * See CharBuffer_copy_from_subset() in CharBuffer.c for a description of the first 3 arguments.
 */
SEXP CharBuffer_translate_copy_from_subset(SEXP dest_xp, SEXP src_xp, SEXP subset, SEXP hash_xp)
{
	SEXP dest, src, hash;
	char hash_hole;

	dest = R_ExternalPtrTag(dest_xp);
	src = R_ExternalPtrTag(src_xp);
	hash = R_ExternalPtrTag(hash_xp);
	hash_hole = CHAR(hash)[0];
	Biostrings_translate_charcpy_from_subset(INTEGER(subset), LENGTH(subset),
		CHAR(dest), LENGTH(dest), CHAR(src), LENGTH(src),
		CHAR(hash) + 1, LENGTH(hash) - 1, hash_hole, 0);
	return dest_xp;
}

/*
 * See CharBuffer_copy_from_i1i2() in CharBuffer.c for a description of the arguments.
 */
SEXP CharBuffer_reverse_copy_from_i1i2(SEXP dest_xp, SEXP src_xp, SEXP imin, SEXP imax)
{
	SEXP dest, src;
	int i1, i2;

	dest = R_ExternalPtrTag(dest_xp);
	i1 = INTEGER(imin)[0] - 1;
	i2 = INTEGER(imax)[0] - 1;
	src = R_ExternalPtrTag(src_xp);
	Biostrings_reverse_memcpy_from_i1i2(i1, i2,
		CHAR(dest), LENGTH(dest),
		CHAR(src), LENGTH(src), sizeof(char));
	return dest_xp;
}

/*
 * See CharBuffer_copy_from_i1i2() in CharBuffer.c for a description of the first 4 arguments.
 */
SEXP CharBuffer_reverse_translate_copy_from_i1i2(SEXP dest_xp, SEXP src_xp, SEXP imin, SEXP imax, SEXP hash_xp)
{
	SEXP dest, src, hash;
	int i1, i2;
	char hash_hole;

	dest = R_ExternalPtrTag(dest_xp);
	i1 = INTEGER(imin)[0] - 1;
	i2 = INTEGER(imax)[0] - 1;
	src = R_ExternalPtrTag(src_xp);
	/* hash = STRING_ELT(hash_xp, 0); */
	hash = R_ExternalPtrTag(hash_xp);
	hash_hole = CHAR(hash)[0];
        Biostrings_reverse_translate_charcpy_from_i1i2(i1, i2,
		CHAR(dest), LENGTH(dest), CHAR(src), LENGTH(src),
                CHAR(hash) + 1, LENGTH(hash) - 1, hash_hole, 0);
	return dest_xp;
}

