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
 * C-level slot getters for PreprocessedTB objects.
 *
 * Be careful that these functions do NOT duplicate the returned slot.
 * Thus they cannot be made .Call() entry points!
 */

static SEXP
	tb_symbol = NULL,
	dups_symbol = NULL,
	base_codes_symbol = NULL;

static SEXP get_PreprocessedTB_tb(SEXP x)
{
	INIT_STATIC_SYMBOL(tb)
	return GET_SLOT(x, tb_symbol);
}

static SEXP get_PreprocessedTB_dups(SEXP x)
{
	INIT_STATIC_SYMBOL(dups)
	return GET_SLOT(x, dups_symbol);
}

SEXP _get_PreprocessedTB_base_codes(SEXP x)
{
	INIT_STATIC_SYMBOL(base_codes)
	return GET_SLOT(x, base_codes_symbol);
}

/* Not strict "slot getters" but very much like. */

int _get_PreprocessedTB_length(SEXP x)
{
	return _get_XStringSet_length(get_PreprocessedTB_tb(x));
}

int _get_PreprocessedTB_width(SEXP x)
{
	SEXP tb;

	tb = get_PreprocessedTB_tb(x);
	return INTEGER(get_IRanges_width(_get_XStringSet_ranges(tb)))[0];
}

SEXP _get_PreprocessedTB_low2high(SEXP x)
{
	return get_H2LGrouping_low2high(get_PreprocessedTB_dups(x));
}


/****************************************************************************
 * C-level slot getters for Twobit objects.
 *
 * Be careful that these functions do NOT duplicate the returned slot.
 * Thus they cannot be made .Call() entry points!
 */

static SEXP sign2pos_symbol = NULL;

static SEXP get_Twobit_sign2pos(SEXP x)
{
	INIT_STATIC_SYMBOL(sign2pos)
	return GET_SLOT(x, sign2pos_symbol);
}

/* Not a strict "slot getter" but very much like. */
SEXP _get_Twobit_sign2pos_tag(SEXP x)
{
	return get_XSequence_tag(get_Twobit_sign2pos(x));
}


/****************************************************************************
 * C-level slot getters for ACtree objects.
 *
 * Be careful that these functions do NOT duplicate the returned slot.
 * Thus they cannot be made .Call() entry points!
 */

static SEXP nodes_symbol = NULL;

static SEXP get_ACtree_nodes(SEXP x)
{
	INIT_STATIC_SYMBOL(nodes)
	return GET_SLOT(x, nodes_symbol);
}

/* Not a strict "slot getter" but very much like. */
SEXP _get_ACtree_nodes_tag(SEXP x)
{
	return get_XSequence_tag(get_ACtree_nodes(x));
}


/****************************************************************************
 * C-level slot getters for ACtree2 objects.
 *
 * Be careful that these functions do NOT duplicate the returned slot.
 * Thus they cannot be made .Call() entry points!
 */

static SEXP
	nodebuf_ptr_symbol = NULL,
	nodeextbuf_ptr_symbol = NULL;

SEXP _get_ACtree2_nodebuf_ptr(SEXP x)
{
	INIT_STATIC_SYMBOL(nodebuf_ptr)
	return GET_SLOT(x, nodebuf_ptr_symbol);
}

SEXP _get_ACtree2_nodeextbuf_ptr(SEXP x)
{
	INIT_STATIC_SYMBOL(nodeextbuf_ptr)
	return GET_SLOT(x, nodeextbuf_ptr_symbol);
}


/****************************************************************************
 * Buffer of duplicates.
 */

static IntAE ppdups_buf;

void _init_ppdups_buf(int length)
{
	ppdups_buf = new_IntAE(length, length, NA_INTEGER);
	return;
}

void _report_ppdup(int poffset, int P_id)
{
	ppdups_buf.elts[poffset] = P_id;
	return;
}

SEXP _get_ppdups_buf_asINTEGER()
{
	return IntAE_asINTEGER(&ppdups_buf);
}

