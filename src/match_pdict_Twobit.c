/****************************************************************************
 *                           The Twobit algorithm                           *
 *                   for constant width DNA dictionaries                    *
 *                                                                          *
 *                           Author: Herve Pages                            *
 *                                                                          *
 * Note: a constant width dictionary is a non-empty set of non-empty        *
 * words of the same length.                                                *
 ****************************************************************************/
#include "Biostrings.h"

static int debug = 0;


SEXP debug_match_pdict_Twobit()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'match_pdict_Twobit.c'\n",
		debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'match_pdict_Twobit.c'\n");
#endif
	return R_NilValue;
}



/****************************************************************************
 *                                                                          *
 *                             A. PREPROCESSING                             *
 *                                                                          *
 ****************************************************************************/

/****************************************************************************
 * Building the Twobit_PDict object
 * --------------------------------
 */

static void pp_pattern(int poffset)
{
	int width, n;
	const char *pattern;
	char c;

	width = _CroppedDict_width();
	pattern = _CroppedDict_pattern(poffset);
	for (n = 0; n < width; n++) {
		c = pattern[n];
	}
	//report_dup(poffset, node->P_id);
	return;
}

static void build_Twobit_PDict()
{
	int length, width, poffset;

	length = _CroppedDict_length();
	width = _CroppedDict_width();
	init_dup2unq_buf(length);
	for (poffset = 0; poffset < length; poffset++)
		pp_pattern(poffset);
	return;
}


/****************************************************************************
 * Turning our local data structures into an R list (SEXP)
 * -------------------------------------------------------
 *
 * Twobit_PDict_asLIST() returns an R list with the following elements:
 *   - geom: a list describing the geometry of the cropped dictionary;
 *   - dup2unq: an integer vector containing the mapping between duplicated and
 *         primary reads.
 */

static SEXP Twobit_PDict_asLIST()
{
	SEXP ans, ans_names, ans_elt;

	PROTECT(ans = NEW_LIST(2));

	/* set the names */
	PROTECT(ans_names = NEW_CHARACTER(2));
	SET_STRING_ELT(ans_names, 0, mkChar("geom"));
	SET_STRING_ELT(ans_names, 1, mkChar("dup2unq"));
	SET_NAMES(ans, ans_names);
	UNPROTECT(1);

	/* set the "geom" element */
	PROTECT(ans_elt = _CroppedDict_geom_asLIST());
	SET_ELEMENT(ans, 0, ans_elt);
	UNPROTECT(1);

	/* set the "dup2unq" element */
	PROTECT(ans_elt = _dup2unq_asINTEGER());
	SET_ELEMENT(ans, 1, ans_elt);
	UNPROTECT(1);

	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * .Call entry points for preprocessing
 * ------------------------------------
 *
 * Arguments:
 *   'dict': a character vector or DNAStringSet object containing the input
 *           sequences
 *   'tb_start': single >= 1, <= -1 or NA integer
 *   'tb_end': single >= 1, <= -1 or NA integer
 *
 * See Twobit_PDict_asLIST() for a description of the returned SEXP.
 */

SEXP build_Twobit_PDict_from_CHARACTER(SEXP dict, SEXP tb_start, SEXP tb_end)
{
	_init_CroppedDict_with_CHARACTER(dict,
			INTEGER(tb_start)[0], INTEGER(tb_end)[0]);
	build_Twobit_PDict();
	return Twobit_PDict_asLIST();
}

SEXP build_Twobit_PDict_from_XStringSet(SEXP dict, SEXP tb_start, SEXP tb_end)
{
	_init_CroppedDict_with_XStringSet(dict,
			INTEGER(tb_start)[0], INTEGER(tb_end)[0]);
	build_Twobit_PDict();
	return Twobit_PDict_asLIST();
}



/****************************************************************************
 *                                                                          *
 *                             B. MATCH FINDING                             *
 *                                                                          *
 ****************************************************************************/

void _match_Twobit_PDict(SEXP pdict_data, const RoSeq *S, int fixedS)
{
#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] ENTERING _match_Twobit_PDict()\n");
#endif

#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] LEAVING _match_Twobit_PDict()\n");
#endif
	return;
}

