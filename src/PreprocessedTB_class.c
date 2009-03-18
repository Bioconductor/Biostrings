/****************************************************************************
 *               Basic manipulation of PreprocessedTB objects               *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"
#include "IRanges_interface.h"

static int debug = 0;

SEXP debug_PreprocessedTB_class()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in file %s\n",
		debug ? "on" : "off", __FILE__);
#else
	Rprintf("Debug mode not available in file %s\n", __FILE__);
#endif
	return R_NilValue;
}


/****************************************************************************
 * Accessor functions for PreprocessedTB objects.
 */

int _get_PreprocessedTB_length(SEXP x)
{
	SEXP tb;

	tb = GET_SLOT(x, install("tb"));
	return _get_XStringSet_length(tb);
}

int _get_PreprocessedTB_width(SEXP x)
{
	SEXP tb;

	tb = GET_SLOT(x, install("tb"));
	return INTEGER(get_IRanges_width(_get_XStringSet_ranges(tb)))[0];
}

/* Be careful that this function does NOT copy the returned slot! (this
 * prevents it from being made .Call() entry point) */
SEXP _get_PreprocessedTB_unq2dup(SEXP x)
{
	SEXP dups;

	dups = GET_SLOT(x, install("dups"));
	return GET_SLOT(dups, install("unq2dup"));
}


/****************************************************************************
 * Accessor functions for Twobit objects.
 *
 * Be careful that these functions do NOT copy the returned slot! (this
 * prevents them from being made .Call() entry points)
 */

SEXP _get_Twobit_sign2pos_tag(SEXP x)
{
	return get_XSequence_tag(GET_SLOT(x, install("sign2pos")));
}

SEXP _get_Twobit_base_codes(SEXP x)
{
	return GET_SLOT(x, install("base_codes"));
}


/****************************************************************************
 * Accessor functions for ACtree objects.
 *
 * Be careful that these functions do NOT copy the returned slot! (this
 * prevents them from being made .Call() entry points)
 */

SEXP _get_ACtree_nodes_tag(SEXP x)
{
	return get_XSequence_tag(GET_SLOT(x, install("nodes")));
}

SEXP _get_ACtree_base_codes(SEXP x)
{
	return GET_SLOT(x, install("base_codes"));
}


/****************************************************************************
 * Accessor functions for ACtree2 objects.
 *
 * Be careful that these functions do NOT copy the returned slot! (this
 * prevents them from being made .Call() entry points)
 */

SEXP _get_ACtree2_nodebuf_ptr(SEXP x)
{
	return GET_SLOT(x, install("nodebuf_ptr"));
}

SEXP _get_ACtree2_nodeextbuf_ptr(SEXP x)
{
	return GET_SLOT(x, install("nodeextbuf_ptr"));
}

SEXP _get_ACtree2_base_codes(SEXP x)
{
	return GET_SLOT(x, install("base_codes"));
}

