/****************************************************************************
 *                  Basic manipulation of IRanges objects                   *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"

static int debug = 0;

SEXP debug_IRanges_class()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'IRanges_class.c'\n",
	        debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'IRanges_class.c'\n");
#endif
	return R_NilValue;
}

SEXP _get_IRanges_start(SEXP x)
{
	return GET_SLOT(x, install("start"));
}

SEXP _get_IRanges_width(SEXP x)
{
	return GET_SLOT(x, install("width"));
}

int _get_IRanges_length(SEXP x)
{
	return LENGTH(_get_IRanges_start(x));
}

const int *_get_IRanges_start0(SEXP x)
{
	return INTEGER(_get_IRanges_start(x));
}

const int *_get_IRanges_width0(SEXP x)
{
	return INTEGER(_get_IRanges_width(x));
}

/*
 * Does NOT duplicate 'x'. The @NAMES slot is modified in place!
 */
void _set_IRanges_names(SEXP x, SEXP names)
{
	SEXP names_slot;

	if (names == R_NilValue) {
		PROTECT(names_slot = NEW_CHARACTER(1));
		SET_STRING_ELT(names_slot, 0, NA_STRING);
		SET_SLOT(x, mkChar("NAMES"), names_slot);
		UNPROTECT(1);
	} else {
		if (LENGTH(names) != _get_IRanges_length(x))
			error("number of names and number of elements differ");
		SET_SLOT(x, mkChar("NAMES"), names);
	}
	return;
}

/*
 * Note that 'start' and 'width' must NOT contain NAs.
 * set_IRanges_slots() trusts the caller and does NOT check this!
 */
static void set_IRanges_slots(SEXP x, SEXP start, SEXP width, SEXP names)
{
	if (LENGTH(width) != LENGTH(start))
		error("number of starts and number of widths differ");
	SET_SLOT(x, mkChar("start"), start);
	SET_SLOT(x, mkChar("width"), width);
	_set_IRanges_names(x, names);
	return;
}

void _copy_IRanges_slots(SEXP x, SEXP x0)
{
	SET_SLOT(x, mkChar("start"), duplicate(GET_SLOT(x0, install("start"))));
	SET_SLOT(x, mkChar("width"), duplicate(GET_SLOT(x0, install("width"))));
	SET_SLOT(x, mkChar("NAMES"), duplicate(GET_SLOT(x0, install("NAMES"))));
	return;
}

/*
 * Never try to make this a .Call() entry point!
 * Its arguments are NOT duplicated so it would be a disaster if they were
 * coming from the user space.
 */
SEXP _new_IRanges(const char *class, SEXP start, SEXP width, SEXP names)
{
	SEXP class_def, ans;

	class_def = MAKE_CLASS(class);
	PROTECT(ans = NEW_OBJECT(class_def));
	set_IRanges_slots(ans, start, width, names);
	UNPROTECT(1);
	return ans;
}

SEXP _new_IRanges_from_RoSeqs(const char *class, RoSeqs seqs)
{
	const RoSeq *seq;
	SEXP start, width, ans;
	int *start_elt, *width_elt, *start_prev_elt, i;

#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _new_IRanges_from_RoSeqs(): BEGIN\n");
	}
#endif
	seq = seqs.elts;
	PROTECT(start = NEW_INTEGER(seqs.nelt));
	PROTECT(width = NEW_INTEGER(seqs.nelt));
	start_elt = INTEGER(start);
	width_elt = INTEGER(width);
	if (seqs.nelt >= 1) {
		*(start_elt++) = 1;
		*(width_elt++) = seq->nelt;
	}
	if (seqs.nelt >= 2)
		for (i = 1, start_prev_elt = INTEGER(start); i < seqs.nelt; i++) {
			*(start_elt++) = *(start_prev_elt++) + (seq++)->nelt;
			*(width_elt++) = seq->nelt;
		}
	PROTECT(ans = _new_IRanges(class, start, width, R_NilValue));
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _new_IRanges_from_RoSeqs(): END\n");
	}
#endif
	UNPROTECT(3);
	return ans;
}

/*
 * Allocation WITHOUT initialization.
 * The 'start' and 'width' slots are not initialized (they contain junk).
 */
SEXP _alloc_IRanges(const char *class, int length)
{
        SEXP start, width, ans;

        PROTECT(start = NEW_INTEGER(length));
        PROTECT(width = NEW_INTEGER(length));
        PROTECT(ans = _new_IRanges(class, start, width, R_NilValue));
        UNPROTECT(3);
        return ans;
}

