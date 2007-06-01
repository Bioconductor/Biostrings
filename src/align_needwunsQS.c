#include "Biostrings.h"

#include <stdio.h>

/*
 * 's1_xp', 's1_offset', 's1_length': left BString object
 * 's2_xp', 's2_offset', 's2_length': right BString object
 * 'mat': scoring matrix (integer vector)
 * 'gapcost': gap cost
 * align_needwunsQS() returns a list of 2 BString objects.
 * Note that the 2 BString objects to align should contain no gaps.
 */
SEXP align_needwunsQS(SEXP s1_xp, SEXP s1_offset, SEXP s1_length,
		SEXP s2_xp, SEXP s2_offset, SEXP s2_length,
		SEXP mat, SEXP gapcost)
{
	SEXP ans;

	return ans;
}

