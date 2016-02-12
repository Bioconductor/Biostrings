/****************************************************************************
 *               Basic manipulation of PreprocessedTB objects               *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "Biostrings.h"
#include "XVector_interface.h"
#include "IRanges_interface.h"
#include "S4Vectors_interface.h"


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

SEXP _get_PreprocessedTB_tb(SEXP x)
{
	INIT_STATIC_SYMBOL(tb)
	return GET_SLOT(x, tb_symbol);
}

SEXP _get_PreprocessedTB_dups(SEXP x)
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
	return _get_XStringSet_length(_get_PreprocessedTB_tb(x));
}

int _get_PreprocessedTB_width(SEXP x)
{
	SEXP tb;

	tb = _get_PreprocessedTB_tb(x);
	return INTEGER(_get_XStringSet_width(tb))[0];
}

SEXP _get_PreprocessedTB_low2high(SEXP x)
{
	return get_H2LGrouping_low2high(_get_PreprocessedTB_dups(x));
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
	return get_XVector_tag(get_Twobit_sign2pos(x));
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

static IntAE *ppdups_buf;

void _init_ppdups_buf(int length)
{
	ppdups_buf = new_IntAE(length, length, NA_INTEGER);
	return;
}

void _report_ppdup(int poffset, int P_id)
{
	ppdups_buf->elts[poffset] = P_id;
	return;
}

SEXP _get_ppdups_buf_asINTEGER()
{
	return new_INTEGER_from_IntAE(ppdups_buf);
}

