#include "Biostrings.h"


/*
 * See CharBuffer_copy_from_i1i2() in CharBuffer.c for a description of the first 4 arguments.
 */
SEXP CharBuffer_translate_copy_from_i1i2(SEXP dest_xp, SEXP src_xp, SEXP imin, SEXP imax, SEXP lkup)
{
	SEXP dest, src;
	int i1, i2;

	dest = R_ExternalPtrTag(dest_xp);
	i1 = INTEGER(imin)[0] - 1;
	i2 = INTEGER(imax)[0] - 1;
	src = R_ExternalPtrTag(src_xp);
        Biostrings_translate_charcpy_from_i1i2(i1, i2,
		CHAR(dest), LENGTH(dest), CHAR(src), LENGTH(src),
                INTEGER(lkup), LENGTH(lkup));
	return dest_xp;
}

/*
 * See CharBuffer_copy_from_subset() in CharBuffer.c for a description of the first 3 arguments.
 */
SEXP CharBuffer_translate_copy_from_subset(SEXP dest_xp, SEXP src_xp, SEXP subset, SEXP lkup)
{
	SEXP dest, src;

	dest = R_ExternalPtrTag(dest_xp);
	src = R_ExternalPtrTag(src_xp);
	Biostrings_translate_charcpy_from_subset(INTEGER(subset), LENGTH(subset),
		CHAR(dest), LENGTH(dest), CHAR(src), LENGTH(src),
		INTEGER(lkup), LENGTH(lkup));
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
SEXP CharBuffer_reverse_translate_copy_from_i1i2(SEXP dest_xp, SEXP src_xp, SEXP imin, SEXP imax, SEXP lkup)
{
	SEXP dest, src;
	int i1, i2;

	dest = R_ExternalPtrTag(dest_xp);
	i1 = INTEGER(imin)[0] - 1;
	i2 = INTEGER(imax)[0] - 1;
	src = R_ExternalPtrTag(src_xp);
        Biostrings_reverse_translate_charcpy_from_i1i2(i1, i2,
		CHAR(dest), LENGTH(dest), CHAR(src), LENGTH(src),
                INTEGER(lkup), LENGTH(lkup));
	return dest_xp;
}

