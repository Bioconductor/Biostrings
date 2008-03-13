/****************************************************************************
 *      Converting a set of sequences from one internal representation      *
 *                               into another                               *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"
#include <S.h> /* for Salloc() and Srealloc() */

static int debug = 0;

SEXP Biostrings_debug_seqs_to_seqs()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'seqs_to_seqs.c'\n", debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'seqs_to_seqs.c'\n");
#endif
	return R_NilValue;
}

/*
 * IMPORTANT: All the functions in this file assume that their 'safe_starts'
 * and 'safe_widths' arguments have been obtained by processing the
 * user-specified Start/End/Nchar values thru R function restrict().
 * In other words those arguments are assumed to be safe i.e. they should
 * describe a set of valid locations in their first argument 'x'.
 */

static CharAArr new_CharAArr(int nelt)
{
	CharAArr arr;

	arr.elts = Salloc((long) nelt, CharArr);
	arr.nelt = nelt;
	return arr;
}

/*
 * Should never raise an error.
 */
static int start2offset(int safe_start)
{
	if (safe_start < 1)
		error("Biostrings internal error in start2offset(): safe_start < 1");
	return --safe_start;
}


/****************************************************************************
 * Converting a set of sequences into a CharAArr struct (array of arrays of
 * const chars).
 * For all these functions 'safe_starts' must either be NULL of have at least
 * 'nelt' elements. Same for 'safe_widths'.
 */

CharAArr _new_CharAArr_from_BBuf(CharBBuf cbbuf)
{
	CharAArr arr;
	CharArr *elt1;
	CharBuf *elt2;
	int i;

	arr = new_CharAArr(cbbuf.nelt);
	for (i = 0, elt1 = arr.elts, elt2 = cbbuf.elts;
	     i < cbbuf.nelt;
	     i++, elt1++, elt2++) {
		elt1->elts = elt2->elts;
		elt1->nelt = elt2->nelt;
	}
	return arr;
}

CharAArr _new_CharAArr_from_STRSXP(int nelt, SEXP x,
		const int *safe_starts, const int *safe_widths)
{
	CharAArr arr;
	CharArr *elt1;
	SEXP elt2;
	int i;
	const int *safe_start, *safe_width;

	if (nelt > LENGTH(x))
		error("_new_CharAArr_from_STRSXP(): 'nelt' must be <= 'LENGTH(x)'");
	arr = new_CharAArr(nelt);
	safe_start = safe_starts;
	safe_width = safe_widths;
	for (i = 0, elt1 = arr.elts; i < nelt; i++, elt1++) {
		elt2 = STRING_ELT(x, i);
		if (elt2 == NA_STRING)
			error("input sequence %d is NA", i+1);
		elt1->elts = CHAR(elt2);
		if (safe_start != NULL) {
			elt1->elts += start2offset(*safe_start);
			safe_start++;
		}
		if (safe_width != NULL) {
			elt1->nelt = *safe_width;
			safe_width++;
		} else {
			elt1->nelt = LENGTH(elt2);
		}
	}
	return arr;
}

CharAArr _new_CharAArr_from_XString(int nelt, SEXP x,
		const int *safe_starts, const int *safe_widths)
{
	CharAArr arr;
	CharArr *elt1;
	int i;
	const int *safe_start, *safe_width;

	arr = new_CharAArr(nelt);
	safe_start = safe_starts;
	safe_width = safe_widths;
	for (i = 0, elt1 = arr.elts; i < nelt; i++, elt1++) {
		elt1->elts = _get_XString_charseq(x, &(elt1->nelt));
		if (safe_start != NULL) {
			elt1->elts += start2offset(*safe_start);
			safe_start++;
		}
		if (safe_width != NULL) {
			elt1->nelt = *safe_width;
			safe_width++;
		}
	}
	return arr;
}

CharAArr _new_CharAArr_from_XStringSet(int nelt, SEXP x,
		const int *safe_starts, const int *safe_widths)
{
	CharAArr arr;
	CharArr *elt1;
	int i;
	const int *safe_start, *safe_width;

	if (nelt > _get_XStringSet_length(x))
		error("_new_CharAArr_from_XStringSet(): 'nelt' must be <= '_get_XStringSet_length(x)'");
	arr = new_CharAArr(nelt);
	safe_start = safe_starts;
	safe_width = safe_widths;
	for (i = 0, elt1 = arr.elts; i < nelt; i++, elt1++) {
		elt1->elts = _get_XStringSet_charseq(x, i, &(elt1->nelt));
		if (safe_start != NULL) {
			elt1->elts += start2offset(*safe_start);
			safe_start++;
		}
		if (safe_width != NULL) {
			elt1->nelt = *safe_width;
			safe_width++;
		}
	}
	return arr;
}

CharAArr _new_CharAArr_from_XStringList(int nelt, SEXP x,
		const int *safe_starts, const int *safe_widths)
{
	CharAArr arr;
	CharArr *elt1;
	int i;
	const int *safe_start, *safe_width;

	if (nelt > _get_XStringList_length(x))
		error("_new_CharAArr_from_XStringList(): 'nelt' must be <= '_get_XStringList_length(x)'");
	arr = new_CharAArr(nelt);
	safe_start = safe_starts;
	safe_width = safe_widths;
	for (i = 0, elt1 = arr.elts; i < nelt; i++, elt1++) {
		elt1->elts = _get_XStringList_charseq(x, i, &(elt1->nelt));
		if (safe_start != NULL) {
			elt1->elts += start2offset(*safe_start);
			safe_start++;
		}
		if (safe_width != NULL) {
			elt1->nelt = *safe_width;
			safe_width++;
		}
	}
	return arr;
}


/****************************************************************************
 * Converting a set of sequences into an XRaw object.
 */

/* NOT a Call() entry point! */
SEXP mkXRaw(SEXP tag)
{
	SEXP ans;

	PROTECT(ans = NEW_OBJECT(MAKE_CLASS("XRaw")));
	SET_SLOT(ans, mkChar("xp"), R_MakeExternalPtr(NULL, tag, R_NilValue));
	UNPROTECT(1);
        return ans;
}

static SEXP new_RAW_from_CharAArr(CharAArr arr, SEXP lkup)
{
	SEXP ans;
	int ans_length, i;
	CharArr *arr_elt;
	char *dest;

	ans_length = 0;
	for (i = 0, arr_elt = arr.elts; i < arr.nelt; i++, arr_elt++)
		ans_length += arr_elt->nelt;
	PROTECT(ans = NEW_RAW(ans_length));
	dest = (char *) RAW(ans);
	for (i = 0, arr_elt = arr.elts; i < arr.nelt; i++, arr_elt++) {
		if (lkup == R_NilValue) {
			_Biostrings_memcpy_to_i1i2(0, arr_elt->nelt - 1,
				dest, arr_elt->nelt,
				arr_elt->elts, arr_elt->nelt, sizeof(char));
		} else {
			_Biostrings_translate_charcpy_to_i1i2(0, arr_elt->nelt - 1,
				dest, arr_elt->nelt,
				arr_elt->elts, arr_elt->nelt,
				INTEGER(lkup), LENGTH(lkup));
		}
		dest += arr_elt->nelt;
	}
	UNPROTECT(1);
	return ans;
}

/*
 * --- .Call ENTRY POINT ---
 */
SEXP copy_subXRaw(SEXP x, SEXP start, SEXP nchar, SEXP lkup)
{
	SEXP ans;

	error("copy_subXRaw() not ready yet");
	return R_NilValue;
}

/*
 * --- .Call ENTRY POINT ---
 * TODO: Support the 'collapse' argument
 */
SEXP STRSXP_to_XRaw(SEXP x, SEXP safe_starts, SEXP safe_widths, SEXP collapse, SEXP lkup)
{
	int nseq;
	CharAArr arr;
	SEXP tag, ans;

	nseq = LENGTH(safe_starts);
	if (collapse == R_NilValue) {
		if (nseq != 1)
			error("'collapse' must be specified when the number "
			      "of input sequences is not exactly 1");
	} else {
		if (LENGTH(collapse) != 1 || LENGTH(STRING_ELT(collapse, 0)) != 0)
			error("'collapse' can only be NULL or the empty string for now");
	}
	arr = _new_CharAArr_from_STRSXP(nseq, x, INTEGER(safe_starts), INTEGER(safe_widths));
	PROTECT(tag = new_RAW_from_CharAArr(arr, lkup));
	ans = mkXRaw(tag);
	UNPROTECT(1);
	return ans;
}

/*
 * --- .Call ENTRY POINT ---
 */
SEXP XString_to_XRaw(SEXP x, SEXP safe_starts, SEXP safe_widths, SEXP lkup)
{
	int nseq;
	CharAArr arr;
	SEXP tag, ans;

	nseq = LENGTH(safe_starts);
	arr = _new_CharAArr_from_XString(nseq, x, INTEGER(safe_starts), INTEGER(safe_widths));
	PROTECT(tag = new_RAW_from_CharAArr(arr, lkup));
	ans = mkXRaw(tag);
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * Converting a set of sequences into an XStringSet object.
 */

SEXP _new_XStringSet_from_seqsnames(const char *baseClass,
		CharAArr seqs, CharAArr names)
{
	SEXP ans;

	error("_new_XStringSet_from_seqsnames() is not ready yet");
	return ans;
}


/****************************************************************************
 * Converting a set of sequences into an XStringList object.
 */

static SEXP mkXStringList(const char *class, SEXP seqs)
{
	SEXP class_def, ans;

	class_def = MAKE_CLASS(class);
	PROTECT(ans = NEW_OBJECT(class_def));
	SET_SLOT(ans, mkChar("seqs"), seqs);
	UNPROTECT(1);
	return ans;
}

/*
 * --- .Call ENTRY POINT ---
 */
SEXP XRaw_to_XStringList(SEXP x, SEXP safe_starts, SEXP safe_widths, SEXP proto)
{
	int nseq, x_length, i;
	SEXP ans, ans_seqs, ans_seq, data;
	const int *safe_start, *safe_width;
	char classbuf[14]; // longest string will be "DNAStringList"

	if (LENGTH(safe_starts) != LENGTH(safe_widths))
		error("'safe_starts' and 'safe_widths' must have the same length");
	nseq = LENGTH(safe_starts);
	x_length = LENGTH(R_ExternalPtrTag(GET_SLOT(x, install("xp"))));
	PROTECT(ans_seqs = NEW_LIST(nseq));
	for (i = 0, safe_start = INTEGER(safe_starts), safe_width = INTEGER(safe_widths);
	     i < nseq;
	     i++, safe_start++, safe_width++) {
		PROTECT(ans_seq = duplicate(proto));
		PROTECT(data = duplicate(x));
		SET_SLOT(ans_seq, mkChar("data"), data);
		SET_SLOT(ans_seq, mkChar("offset"), ScalarInteger(start2offset(*safe_start)));
		SET_SLOT(ans_seq, mkChar("length"), ScalarInteger(*safe_width));
		SET_ELEMENT(ans_seqs, i, ans_seq);
		UNPROTECT(2);
	}
	snprintf(classbuf, sizeof(classbuf), "%sList", get_class(proto));
	ans = mkXStringList(classbuf, ans_seqs);
	UNPROTECT(1);
	return ans;
}

/*
 * --- .Call ENTRY POINT ---
 */
SEXP narrow_XStringList(SEXP x, SEXP safe_starts, SEXP safe_widths, SEXP proto)
{
	int nseq, i, seq_length, offset;
	SEXP x_seqs, x_seq, ans, ans_seqs, ans_seq, data;
	const int *safe_start, *safe_width;
	char classbuf[14]; // longest string will be "DNAStringList"
	const char *class;

	nseq = _get_XStringList_length(x);
	if (LENGTH(safe_starts) != nseq || LENGTH(safe_widths) != nseq)
		error("invalid length of 'safe_starts' or 'safe_widths'");
	x_seqs = GET_SLOT(x, install("seqs"));
	PROTECT(ans_seqs = NEW_LIST(nseq));
	for (i = 0, safe_start = INTEGER(safe_starts), safe_width = INTEGER(safe_widths);
	     i < nseq;
	     i++, safe_start++, safe_width++) {
	        x_seq = VECTOR_ELT(x_seqs, i);
		_get_XString_charseq(x_seq, &seq_length);
		if (proto == R_NilValue) {
			PROTECT(ans_seq = duplicate(x_seq));
		} else {
			PROTECT(ans_seq = duplicate(proto));
			PROTECT(data = duplicate(GET_SLOT(x_seq, install("data"))));
			SET_SLOT(ans_seq, mkChar("data"), data);
			UNPROTECT(1);
		}
		offset = INTEGER(GET_SLOT(x_seq, install("offset")))[0] + *safe_start - 1;
		SET_SLOT(ans_seq, mkChar("offset"), ScalarInteger(offset));
		SET_SLOT(ans_seq, mkChar("length"), ScalarInteger(*safe_width));
		SET_ELEMENT(ans_seqs, i, ans_seq);
		UNPROTECT(1);
	}
	if (proto == R_NilValue) {
		class = get_class(x);
	} else {
		snprintf(classbuf, sizeof(classbuf), "%sList", get_class(proto));
		class = classbuf;
	}
	ans = mkXStringList(class, ans_seqs);
	UNPROTECT(1);
	return ans;
}

