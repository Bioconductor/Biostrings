/****************************************************************************
 *                          Fast MIndex utilities                           *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"
#include "IRanges_interface.h"

static int debug = 0;

SEXP debug_MIndex_utils()
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
 * Match reporting facilities used by C code performing matching of multiple
 * sequences against a reference sequence.
 * -------------------------------------------------------------------
 */

static int string2code(const char *s)
{
	if (strcmp(s, "MATCHES_AS_NULL") == 0)
		return MATCHES_AS_NULL;
	if (strcmp(s, "MATCHES_AS_WHICH") == 0)
		return MATCHES_AS_WHICH;
	if (strcmp(s, "MATCHES_AS_COUNTS") == 0)
		return MATCHES_AS_COUNTS;
	if (strcmp(s, "MATCHES_AS_ENDS") == 0)
		return MATCHES_AS_ENDS;
	error("\"%s\": unsupported \"matches as\" value", s);
	return -1; /* keeps gcc -Wall happy */
}

Seq2MatchBuf _new_Seq2MatchBuf(SEXP matches_as, int nseq)
{
	static Seq2MatchBuf buf;

	buf.matches_as = string2code(CHAR(STRING_ELT(matches_as, 0)));
	if (buf.matches_as != MATCHES_AS_NULL) {
		buf.matching_keys = new_IntAE(0, 0, 0);
		buf.match_counts = new_IntAE(nseq, nseq, 0);
		buf.match_ends = new_IntAEAE(nseq, nseq);
	}
	return buf;
}

void _flush_Seq2MatchBuf(Seq2MatchBuf *buf)
{
	int i;
	const int *key;

	if (buf->matches_as == MATCHES_AS_NULL)
		return;
	for (i = 0, key = buf->matching_keys.elts;
	     i < buf->matching_keys.nelt;
	     i++, key++)
	{
		buf->match_counts.elts[*key] = 0;
		buf->match_ends.elts[*key].nelt = 0;
	}
	buf->matching_keys.nelt = 0;
	return;
}

void _Seq2MatchBuf_report_match(Seq2MatchBuf *buf, int key, int end)
{
	IntAE *ends_buf;

	if (buf->matches_as == MATCHES_AS_NULL)
		return;
	if (buf->match_counts.elts[key]++ == 0)
		IntAE_insert_at(&(buf->matching_keys),
				buf->matching_keys.nelt, key);
	ends_buf = buf->match_ends.elts + key;
	IntAE_insert_at(ends_buf, ends_buf->nelt, end);
	return;
}

void _Seq2MatchBuf_append_and_flush(Seq2MatchBuf *buf1, Seq2MatchBuf *buf2, int view_offset)
{
	int i;
	const int *key;
	IntAE *ends_buf1, *ends_buf2;

	if (buf1->matches_as == MATCHES_AS_NULL)
		return;
	for (i = 0, key = buf2->matching_keys.elts;
	     i < buf2->matching_keys.nelt;
	     i++, key++)
	{
		if (buf1->match_counts.elts[*key] == 0)
			IntAE_insert_at(&(buf1->matching_keys),
					buf1->matching_keys.nelt, *key);
		buf1->match_counts.elts[*key] += buf2->match_counts.elts[*key];
		ends_buf1 = buf1->match_ends.elts + *key;
		ends_buf2 = buf2->match_ends.elts + *key;
		IntAE_append_shifted_vals(ends_buf1,
			ends_buf2->elts, ends_buf2->nelt, view_offset);
	}
	_flush_Seq2MatchBuf(buf2);
	return;
}

SEXP _Seq2MatchBuf_which_asINTEGER(Seq2MatchBuf *buf)
{
	IntAE_sum_val(&(buf->matching_keys), 1);
	IntAE_qsort(&(buf->matching_keys));
	return IntAE_asINTEGER(&(buf->matching_keys));
}

SEXP _Seq2MatchBuf_counts_asINTEGER(Seq2MatchBuf *buf)
{
	return IntAE_asINTEGER(&(buf->match_counts));
}

SEXP _Seq2MatchBuf_ends_asLIST(Seq2MatchBuf *buf)
{
	return IntAEAE_asLIST(&(buf->match_ends), 1);
}

SEXP _Seq2MatchBuf_as_MIndex(Seq2MatchBuf *buf)
{
	error("_Seq2MatchBuf_as_MIndex(): IMPLEMENT ME!");
	return R_NilValue;
}

SEXP _Seq2MatchBuf_as_SEXP(Seq2MatchBuf *buf, SEXP env)
{
	switch (buf->matches_as) {
	    case MATCHES_AS_NULL:
		return R_NilValue;
	    case MATCHES_AS_WHICH:
		return _Seq2MatchBuf_which_asINTEGER(buf);
	    case MATCHES_AS_COUNTS:
		return _Seq2MatchBuf_counts_asINTEGER(buf);
	    case MATCHES_AS_ENDS:
		if (env != R_NilValue)
			return IntAEAE_toEnvir(&(buf->match_ends), env, 1);
		return _Seq2MatchBuf_ends_asLIST(buf);
	    case MATCHES_AS_MINDEX:
		return _Seq2MatchBuf_as_MIndex(buf);
	}
	error("Biostrings internal error in _Seq2MatchBuf_as_SEXP(): "
	      "unexpected 'buf->matches_as' value %d", buf->matches_as);
	return R_NilValue;
}



/****************************************************************************
 * OTHER MIndex FACILITIES
 */

/*
 * Does *inplace* addition of 'val' to all the elements in 'x' (INTSXP).
 * Never use it if 'x' is coming from the user space!
 */
static void add_val_to_INTEGER(SEXP x, int val)
{
	int i, *x_elt;

	for (i = 0, x_elt = INTEGER(x); i < LENGTH(x); i++, x_elt++)
		*x_elt += val;
	return;
}

/*
 * --- .Call ENTRY POINT ---
 * If 'x_width' is NULL => returns the endIndex (list).
 * Otherwise 'x_width' must be an integer vector of same length as the
 * 'x_ends' list and the startIndex is returned.
 */
SEXP ByPos_MIndex_endIndex(SEXP x_high2low, SEXP x_ends, SEXP x_width)
{
	SEXP ans, ans_elt;
	int i, k1;

	PROTECT(ans = duplicate(x_ends));
	for (i = 0; i < LENGTH(ans); i++) {
		if (LENGTH(x_high2low) != 0
		 && (k1 = INTEGER(x_high2low)[i]) != NA_INTEGER) {
			PROTECT(ans_elt = duplicate(VECTOR_ELT(ans, k1 - 1)));
			SET_ELEMENT(ans, i, ans_elt);
			UNPROTECT(1);
			continue;
		}
		if (x_width == R_NilValue)
			continue;
		ans_elt = VECTOR_ELT(ans, i);
		if (!IS_INTEGER(ans_elt)) // could be NULL
			continue;
		add_val_to_INTEGER(ans_elt, 1 - INTEGER(x_width)[i]);
	}
	UNPROTECT(1);
	return ans;
}

/*
 * --- .Call ENTRY POINT ---
 * All the keys in 'x_ends_envir' must be representing integers left-padded with 0s
 * so they have the same length. This works properly:
     library(Biostrings)
     ends_envir <- new.env(parent=emptyenv())
     ends_envir[['0000000010']] <- -2:1
     ends_envir[['0000000004']] <- 9:6
     .Call("SparseMIndex_endIndex",
           ends_envir, NULL, letters[1:10], TRUE,
           PACKAGE="Biostrings")
     .Call("SparseMIndex_endIndex",
           ends_envir, NULL, letters[1:10], FALSE,
           PACKAGE="Biostrings")
 * but this doesn't:
     ends_envir[['3']] <- 33L
     .Call("SparseMIndex_endIndex",
           ends_envir, NULL, letters[1:10], FALSE,
           PACKAGE="Biostrings")
 */
SEXP SparseMIndex_endIndex(SEXP x_ends_envir, SEXP x_width, SEXP x_names, SEXP all_names)
{
	SEXP ans, ans_elt, ans_names, symbols, end;
	int i, j;
	IntAE poffsets, poffsets_order;

	PROTECT(symbols = R_lsInternal(x_ends_envir, 1));
	poffsets = CHARACTER_asIntAE(symbols, -1);
	if (LOGICAL(all_names)[0]) {
		PROTECT(ans = NEW_LIST(LENGTH(x_names)));
		for (i = 0; i < poffsets.nelt; i++) {
			j = poffsets.elts[i];
			end = _get_val_from_env(STRING_ELT(symbols, i), x_ends_envir, 1);
			PROTECT(ans_elt = duplicate(end));
			if (x_width != R_NilValue)
				add_val_to_INTEGER(ans_elt, 1 - INTEGER(x_width)[j]);
			SET_ELEMENT(ans, j, ans_elt);
			UNPROTECT(1);
		}
		SET_NAMES(ans, duplicate(x_names));
		UNPROTECT(1);
	} else {
		//poffsets_order = new_IntAE(poffsets.nelt, 0, 0);
		//get_int_array_order(poffsets.elts, poffsets.nelt, poffsets_order.elts);
		//poffsets_order.nelt = poffsets.nelt; /* = poffsets_order.buflength */
		PROTECT(ans = NEW_LIST(poffsets.nelt));
		PROTECT(ans_names = NEW_CHARACTER(poffsets.nelt));
		for (i = 0; i < poffsets.nelt; i++) {
			//j = poffsets_order.elts[i];
			j = i;
			end = _get_val_from_env(STRING_ELT(symbols, j), x_ends_envir, 1);
			PROTECT(ans_elt = duplicate(end));
			if (x_width != R_NilValue)
				add_val_to_INTEGER(ans_elt, 1 - INTEGER(x_width)[i]);
			SET_ELEMENT(ans, i, ans_elt);
			UNPROTECT(1);
			SET_STRING_ELT(ans_names, i, duplicate(STRING_ELT(x_names, poffsets.elts[j])));
		}
		SET_NAMES(ans, ans_names);
		UNPROTECT(2);
	}
	UNPROTECT(1);
	return ans;
}

/*
 * --- .Call ENTRY POINT ---
 */
SEXP ByPos_MIndex_combine(SEXP ends_listlist)
{
	int NTB, ans_length, i, j;
	SEXP ans, ans_elt, ends;
	IntAE ends_buf;

	NTB = LENGTH(ends_listlist);
	if (NTB == 0)
		error("nothing to combine");
	ans_length = LENGTH(VECTOR_ELT(ends_listlist, 0));
	for (j = 1; j < NTB; j++)
		if (LENGTH(VECTOR_ELT(ends_listlist, j)) != ans_length)
			error("cannot combine MIndex objects of different lengths");
	ends_buf = new_IntAE(0, 0, 0);
	PROTECT(ans = NEW_LIST(ans_length));
	for (i = 0; i < ans_length; i++) {
		ends_buf.nelt = 0;
		for (j = 0; j < NTB; j++) {
			ends = VECTOR_ELT(VECTOR_ELT(ends_listlist, j), i);
			if (ends == R_NilValue)
				continue;
			IntAE_append(&ends_buf, INTEGER(ends), LENGTH(ends));
		}
		if (ends_buf.nelt == 0)
			continue;
		IntAE_qsort(&ends_buf);
		IntAE_delete_adjdups(&ends_buf);
		PROTECT(ans_elt = IntAE_asINTEGER(&ends_buf));
		SET_ELEMENT(ans, i, ans_elt);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	return ans;
}

