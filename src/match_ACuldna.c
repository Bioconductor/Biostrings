/****************************************************************************
 *              Aho-Corasick for uniform-length DNA dictionary              *
 *                           Author: Herve Pages                            *
 *                                                                          *
 * Note: a uniform-length dictionary is a non-empty set of non-empty        *
 * strings of the same length.                                              *
 ****************************************************************************/
#include "Biostrings.h"
#include <S.h> /* for Srealloc() */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/****************************************************************************/
static int debug = 0;

SEXP match_ACuldna_debug()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'match_ACuldna.c'\n", debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'match_ACuldna.c'\n");
#endif
	return R_NilValue;
}


/****************************************************************************
 * Initialization of the Aho-Corasick 4-ary tree
 * =============================================
 */

static void ACuldna_init()
{
	return;
}


/****************************************************************************
 * Exact matching
 * ======================================================
 */

static int ACuldna_exact_search()
{
	return 0;
}


/****************************************************************************
 * .Call entry points: "ACuldna_init_with_StrVect"
 *                 and "ACuldna_init_with_BStringList"
 *
 * Arguments:
 *   'dict': a string vector (aka character vector) containing the
 *           uniform-length dictionary for ACuldna_init_with_StrVect.
 *           A list of (pattern@data@xp, pattern@offset, pattern@length)
 *           triplets containing the uniform-length dictionary for
 *           ACuldna_init_with_BStringList.
 *
 * Return an R list with the following elements:
 *   - AC_tree: XInteger object containing the Aho-Corasick 4-ary tree built
 *         from 'dict'.
 *   - base_codes: integer vector containing the 4 character codes (ASCII)
 *         attached to the 4 child slots of any node in the AC_tree object.
 *   - dups: an unnamed (and eventually empty) list of integer vectors
 *         containing the indices of the duplicated words found in 'dict'.
 *
 ****************************************************************************/

SEXP ACuldna_init_with_StrVect(SEXP dict)
{
	int subj_offset, subj_length, pat_length, c1, c2, c3, c4;
	const Rbyte *subj;
	SEXP buf, ans, ans_names, ans_elt;

	error("Not ready yet!\n");
	ACuldna_init();

	PROTECT(ans = NEW_LIST(3));
	/* set the names */
	PROTECT(ans_names = NEW_CHARACTER(3));
	SET_STRING_ELT(ans_names, 0, mkChar("AC_tree"));
	SET_STRING_ELT(ans_names, 1, mkChar("base_codes"));
	SET_STRING_ELT(ans_names, 2, mkChar("dups"));
	SET_NAMES(ans, ans_names);
	UNPROTECT(1);
	/* set the "AC_tree" element */
	PROTECT(ans_elt = NEW_NUMERIC(4));
	SET_ELEMENT(ans, 0, ans_elt);
	UNPROTECT(1);
	/* set the "base_codes" element */
	PROTECT(ans_elt = NEW_INTEGER(4));
	SET_ELEMENT(ans, 1, ans_elt);
	UNPROTECT(1);
	/* set the "dups" element */
	PROTECT(ans_elt = NEW_INTEGER(4));
	SET_ELEMENT(ans, 2, ans_elt);
	UNPROTECT(1);

	UNPROTECT(1);
	return ans;
}

SEXP ACuldna_init_with_BStringList(SEXP dict)
{
	SEXP ans;

	error("Not ready yet!\n");
	return ans;
}

