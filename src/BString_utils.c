/*
 * Functions for low-level manipulation of BString objects.
 *
 * The functions defined in this file are NOT .Call methods (but they are used
 * by .Call methods defined in other .c files) so THEY DON'T NEED TO BE
 * REGISTERED in R_init_Biostrings.c. They are prefixed with a "_" (underscore)
 * to emphasize the fact that they are used internally within the Biostrings
 * shared lib.
 */
#include "Biostrings.h"


const char *get_BString_seq(SEXP bstring, int *seq_len)
{
	SEXP xp;
	int offset;

	xp = GET_SLOT(GET_SLOT(bstring, install("data")), install("xp"));
	offset = INTEGER(GET_SLOT(bstring, install("offset")))[0];
	*seq_len = INTEGER(GET_SLOT(bstring, install("length")))[0];
	return (const char *) (RAW(R_ExternalPtrTag(xp)) + offset);
}

