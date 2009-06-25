/****************************************************************************
 *                                                                          *
 *       Low-level utility functions used by the matchPDict() C code        *
 *                           Author: Herve Pages                            *
 *                                                                          *
 ****************************************************************************/
#include "Biostrings.h"
#include "IRanges_interface.h"

#include <limits.h> /* for UINT_MAX */
#include <time.h> /* for clock() and CLOCKS_PER_SEC */


static int debug = 0;

SEXP debug_match_pdict_utils()
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



/*****************************************************************************
 * Low-level manipulation of the MatchPDictBuf buffer
 * --------------------------------------------------
 * The MatchPDictBuf struct is used for storing the matches found by the
 * matchPDict() function.
 */

static int string2code(const char *s)
{
	if (strcmp(s, "MATCHES_AS_NULL") == 0)
		return MATCHES_AS_NULL;
	if (strcmp(s, "MATCHES_AS_WHICH") == 0)
		return MATCHES_AS_WHICH;
	if (strcmp(s, "MATCHES_AS_COUNTS") == 0)
		return MATCHES_AS_COUNTS;
	if (strcmp(s, "MATCHES_AS_STARTS") == 0)
		return MATCHES_AS_ENDS;
	if (strcmp(s, "MATCHES_AS_ENDS") == 0)
		return MATCHES_AS_ENDS;
	error("\"%s\": unsupported \"matches as\" value", s);
	return -1; /* keeps gcc -Wall happy */
}


/****************************************************************************
 * Manipulation of the TBMatchBuf struct.
 */

TBMatchBuf _new_TBMatchBuf(int tb_length, int tb_width,
		const int *head_widths, const int *tail_widths)
{
	static TBMatchBuf buf;

	buf.is_init = 1;
	buf.tb_width = tb_width;
	buf.head_widths = head_widths;
	buf.tail_widths = tail_widths;
	buf.matching_keys = new_IntAE(0, 0, 0);
	buf.match_ends = new_IntAEAE(tb_length, tb_length);
	return buf;
}

void _TBMatchBuf_report_match(TBMatchBuf *buf, int key, int end)
{
	IntAE *end_buf;

	if (!buf->is_init)
		return;
	end_buf = buf->match_ends.elts + key;
	if (end_buf->nelt == 0)
		IntAE_insert_at(&(buf->matching_keys),
				buf->matching_keys.nelt, key);
	IntAE_insert_at(end_buf, end_buf->nelt, end);
	return;
}

void _TBMatchBuf_flush(TBMatchBuf *buf)
{
	int i;
	const int *key;

	if (!buf->is_init)
		return;
	for (i = 0, key = buf->matching_keys.elts;
	     i < buf->matching_keys.nelt;
	     i++, key++)
	{
		buf->match_ends.elts[*key].nelt = 0;
	}
	buf->matching_keys.nelt = 0;
	return;
}


/****************************************************************************
 * Manipulation of the Seq2MatchBuf struct.
 */

Seq2MatchBuf _new_Seq2MatchBuf(SEXP matches_as, int nseq)
{
	int code, count_only;
	static Seq2MatchBuf buf;

	code = string2code(CHAR(STRING_ELT(matches_as, 0)));
	count_only = code == MATCHES_AS_WHICH ||
		     code == MATCHES_AS_COUNTS;
	buf.matching_keys = new_IntAE(0, 0, 0);
	buf.match_counts = new_IntAE(nseq, nseq, 0);
	if (count_only) {
		/* By setting 'buflength' to -1 we indicate that these
		   buffers must not be used */
		buf.match_starts.buflength = -1;
		buf.match_widths.buflength = -1;
	} else {
		buf.match_starts = new_IntAEAE(nseq, nseq);
		buf.match_widths = new_IntAEAE(nseq, nseq);
	}
	return buf;
}

void _Seq2MatchBuf_flush(Seq2MatchBuf *buf)
{
	int i;
	const int *key;

	for (i = 0, key = buf->matching_keys.elts;
	     i < buf->matching_keys.nelt;
	     i++, key++)
	{
		buf->match_counts.elts[*key] = 0;
		if (buf->match_starts.buflength != -1)
			buf->match_starts.elts[*key].nelt = 0;
		if (buf->match_widths.buflength != -1)
			buf->match_widths.elts[*key].nelt = 0;
	}
	buf->matching_keys.nelt = 0;
	return;
}

SEXP _Seq2MatchBuf_which_asINTEGER(Seq2MatchBuf *buf)
{
	SEXP ans;
	int i;

	IntAE_qsort(&(buf->matching_keys));
	PROTECT(ans = IntAE_asINTEGER(&(buf->matching_keys)));
	for (i = 0; i < LENGTH(ans); i++)
		INTEGER(ans)[i]++;
	UNPROTECT(1);
	return ans;
}

SEXP _Seq2MatchBuf_counts_asINTEGER(Seq2MatchBuf *buf)
{
	return IntAE_asINTEGER(&(buf->match_counts));
}

SEXP _Seq2MatchBuf_starts_asLIST(Seq2MatchBuf *buf)
{
	if (buf->match_starts.buflength == -1)
		error("Biostrings internal error: _Seq2MatchBuf_starts_asLIST() "
		      "was called in the wrong context");
	return IntAEAE_asLIST(&(buf->match_starts), 1);
}

static SEXP _Seq2MatchBuf_starts_toEnvir(Seq2MatchBuf *buf, SEXP env)
{
	if (buf->match_starts.buflength == -1)
		error("Biostrings internal error: _Seq2MatchBuf_starts_toEnvir() "
		      "was called in the wrong context");
	return IntAEAE_toEnvir(&(buf->match_starts), env, 1);
}

SEXP _Seq2MatchBuf_ends_asLIST(Seq2MatchBuf *buf)
{
	if (buf->match_starts.buflength == -1
	 || buf->match_widths.buflength == -1)
		error("Biostrings internal error: _Seq2MatchBuf_ends_asLIST() "
		      "was called in the wrong context");
	IntAEAE_sum_and_shift(&(buf->match_starts), &(buf->match_widths), -1);
	return IntAEAE_asLIST(&(buf->match_starts), 1);
}

static SEXP _Seq2MatchBuf_ends_toEnvir(Seq2MatchBuf *buf, SEXP env)
{
	if (buf->match_starts.buflength == -1
	 || buf->match_widths.buflength == -1)
		error("Biostrings internal error: _Seq2MatchBuf_ends_toEnvir() "
		      "was called in the wrong context");
	IntAEAE_sum_and_shift(&(buf->match_starts), &(buf->match_widths), -1);
	return IntAEAE_toEnvir(&(buf->match_starts), env, 1);
}

SEXP _Seq2MatchBuf_as_MIndex(Seq2MatchBuf *buf)
{
	error("_Seq2MatchBuf_as_MIndex(): IMPLEMENT ME!");
	return R_NilValue;
}

SEXP _Seq2MatchBuf_as_SEXP(int matches_as, Seq2MatchBuf *buf, SEXP env)
{
	switch (matches_as) {
	    case MATCHES_AS_NULL:
		return R_NilValue;
	    case MATCHES_AS_WHICH:
		return _Seq2MatchBuf_which_asINTEGER(buf);
	    case MATCHES_AS_COUNTS:
		return _Seq2MatchBuf_counts_asINTEGER(buf);
	    case MATCHES_AS_STARTS:
		if (env != R_NilValue)
			return _Seq2MatchBuf_starts_toEnvir(buf, env);
		return _Seq2MatchBuf_starts_asLIST(buf);
	    case MATCHES_AS_ENDS:
		if (env != R_NilValue)
			return _Seq2MatchBuf_ends_toEnvir(buf, env);
		return _Seq2MatchBuf_ends_asLIST(buf);
	    case MATCHES_AS_MINDEX:
		return _Seq2MatchBuf_as_MIndex(buf);
	}
	error("Biostrings internal error in _Seq2MatchBuf_as_SEXP(): "
	      "unsupported 'matches_as' value %d", matches_as);
	return R_NilValue;
}


/****************************************************************************
 * Manipulation of the MatchPDictBuf struct.
 */

MatchPDictBuf _new_MatchPDictBuf(SEXP matches_as, int nseq, int tb_width,
		const int *head_widths, const int *tail_widths)
{
	static MatchPDictBuf buf;

	buf.matches_as = string2code(CHAR(STRING_ELT(matches_as, 0)));
	if (buf.matches_as == MATCHES_AS_NULL) {
		buf.tb_matches.is_init = 0;
	} else {
		buf.tb_matches = _new_TBMatchBuf(nseq, tb_width, head_widths, tail_widths);
		buf.matches = _new_Seq2MatchBuf(matches_as, nseq);
	}
	return buf;
}

void _MatchPDictBuf_report_match(MatchPDictBuf *buf, int key, int tb_end)
{
	IntAE *matching_keys, *count_buf, *start_buf, *width_buf;
	int start, width;

	if (buf->matches_as == MATCHES_AS_NULL)
		return;
	matching_keys = &(buf->matches.matching_keys);
	count_buf = &(buf->matches.match_counts);
	if (count_buf->elts[key]++ == 0)
		IntAE_insert_at(matching_keys, matching_keys->nelt, key);
	width = buf->tb_matches.tb_width;
	start = tb_end - width + 1;
	if (buf->tb_matches.head_widths != NULL) {
		start -= buf->tb_matches.head_widths[key];
		width += buf->tb_matches.head_widths[key];
	}
	if (buf->tb_matches.tail_widths != NULL)
		width += buf->tb_matches.tail_widths[key];
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _MatchPDictBuf_report_match():\n");
		Rprintf("  key=%d  tb_end=%d  start=%d  width=%d\n",
			key, tb_end, start, width);
	}
#endif
	if (buf->matches.match_starts.buflength != -1) {
		start_buf = buf->matches.match_starts.elts + key;
		IntAE_insert_at(start_buf, start_buf->nelt, start);
	}
	if (buf->matches.match_widths.buflength != -1) {
		width_buf = buf->matches.match_widths.elts + key;
		IntAE_insert_at(width_buf, width_buf->nelt, width);
	}
	return;
}

void _MatchPDictBuf_flush(MatchPDictBuf *buf)
{
	if (buf->matches_as == MATCHES_AS_NULL)
		return;
	_TBMatchBuf_flush(&(buf->tb_matches));
	_Seq2MatchBuf_flush(&(buf->matches));
	return;
}

void _MatchPDictBuf_append_and_flush(Seq2MatchBuf *buf1, MatchPDictBuf *buf2,
		int view_offset)
{
	Seq2MatchBuf *buf2_matches;
	int i;
	const int *key;
	IntAE *start_buf1, *start_buf2, *width_buf1, *width_buf2;

	if (buf2->matches_as == MATCHES_AS_NULL)
		return;
	buf2_matches = &(buf2->matches);
	if (buf1->match_counts.nelt != buf2_matches->match_counts.nelt
	 || (buf1->match_starts.buflength == -1) != (buf2_matches->match_starts.buflength == -1)
	 || (buf1->match_widths.buflength == -1) != (buf2_matches->match_widths.buflength == -1))
		error("Biostrings internal error in _MatchPDictBuf_append_and_flush(): "
		      "'buf1' and 'buf2' are incompatible");
	for (i = 0, key = buf2_matches->matching_keys.elts;
	     i < buf2_matches->matching_keys.nelt;
	     i++, key++)
	{
		if (buf1->match_counts.elts[*key] == 0)
			IntAE_insert_at(&(buf1->matching_keys),
					buf1->matching_keys.nelt, *key);
		buf1->match_counts.elts[*key] += buf2_matches->match_counts.elts[*key];
		if (buf1->match_starts.buflength != -1) {
			start_buf1 = buf1->match_starts.elts + *key;
			start_buf2 = buf2_matches->match_starts.elts + *key;
			IntAE_append_shifted_vals(start_buf1,
				start_buf2->elts, start_buf2->nelt, view_offset);
		}
		if (buf1->match_widths.buflength != -1) {
			width_buf1 = buf1->match_widths.elts + *key;
			width_buf2 = buf2_matches->match_widths.elts + *key;
			IntAE_append_shifted_vals(width_buf1,
				width_buf2->elts, width_buf2->nelt, view_offset);
		}
	}
	_MatchPDictBuf_flush(buf2);
	return;
}



/*****************************************************************************
 * Brute force flank comparison routines
 * -------------------------------------
 */

static void collect_keys(int key0, SEXP low2high, IntAE *keys)
{
	SEXP dups;
	int i, *key;

	keys->nelt = 1;
	if (keys->nelt > keys->buflength)
		error("Biostrings internal error in collect_keys(): "
		      "keys->nelt > keys->buflength");
	keys->elts[0] = key0;
	dups = VECTOR_ELT(low2high, key0);
	if (dups == R_NilValue)
		return;
	keys->nelt += LENGTH(dups);
	if (keys->nelt > keys->buflength)
		error("Biostrings internal error in collect_keys(): "
		      "keys->nelt > keys->buflength");
	memcpy(keys->elts + 1, INTEGER(dups), LENGTH(dups) * sizeof(int));
	for (i = 1, key = keys->elts + 1; i < keys->nelt; i++, key++)
		*key -= 1;
	return;
}

static int nmismatch_in_HT(const RoSeq *H, const RoSeq *T,
		const RoSeq *S, int Hshift, int Tshift, int max_mm)
{
	int nmismatch;

	nmismatch = _selected_nmismatch_at_Pshift_fun(H, S, Hshift, max_mm);
	if (nmismatch > max_mm)
		return nmismatch;
	max_mm -= nmismatch;
	nmismatch += _selected_nmismatch_at_Pshift_fun(T, S, Tshift, max_mm);
	return nmismatch;
}

static void match_HT(const RoSeq *H, const RoSeq *T,
		const RoSeq *S, int tb_end, int max_mm,
		MatchPDictBuf *matchpdict_buf, int key)
{
	int HTdeltashift, nmismatch;

	HTdeltashift = H->nelt + matchpdict_buf->tb_matches.tb_width;
	nmismatch = nmismatch_in_HT(H, T,
			S, tb_end - HTdeltashift, tb_end, max_mm);
	if (nmismatch <= max_mm)
		_MatchPDictBuf_report_match(matchpdict_buf, key, tb_end);
	return;
}

static void match_headtail_for_loc(const HeadTail *headtail,
		const RoSeq *S, int tb_end, int max_mm,
		MatchPDictBuf *matchpdict_buf)
{
	int i;
	const int *key;
	const RoSeq *H, *T;

	for (i = 0, key = headtail->keys_buf.elts;
	     i < headtail->keys_buf.nelt;
	     i++, key++)
	{
		H = headtail->head.elts + *key;
		T = headtail->tail.elts + *key;
		match_HT(H, T, S, tb_end, max_mm, matchpdict_buf, *key);
	}
	return;
}

static void match_headtail_for_key(const HeadTail *headtail, int key,
		const RoSeq *S, const IntAE *tb_end_buf, int max_mm,
		MatchPDictBuf *matchpdict_buf)
{
	const RoSeq *H, *T;
	int j;
	const int *tb_end;

	H = headtail->head.elts + key;
	T = headtail->tail.elts + key;
	for (j = 0, tb_end = tb_end_buf->elts;
	     j < tb_end_buf->nelt;
	     j++, tb_end++)
	{
		match_HT(H, T, S, *tb_end, max_mm, matchpdict_buf, key);
	}
	return;
}

static void match_headtail_by_loc(const HeadTail *headtail,
		const RoSeq *S, const IntAE *tb_end_buf, int max_mm,
		MatchPDictBuf *matchpdict_buf)
{
	int j;
	const int *tb_end;

	for (j = 0, tb_end = tb_end_buf->elts;
	     j < tb_end_buf->nelt;
	     j++, tb_end++)
	{
		match_headtail_for_loc(headtail,
				S, *tb_end, max_mm,
				matchpdict_buf);
	}
	return;
}

static void match_headtail_by_key(HeadTail *headtail,
		const RoSeq *S, const IntAE *tb_end_buf, int max_mm,
		MatchPDictBuf *matchpdict_buf)
{
	int i;
	const int *key;

	for (i = 0, key = headtail->keys_buf.elts;
	     i < headtail->keys_buf.nelt;
	     i++, key++)
	{
		match_headtail_for_key(headtail, *key,
				S, tb_end_buf, max_mm,
				matchpdict_buf);
	}
	return;
}



/*****************************************************************************
 * Preprocessing and fast matching of the head and tail of a PDict object
 * ----------------------------------------------------------------------
 *
 * Note that, unlike for the Trusted Band, this is not persistent
 * preprocessing, i.e. the result of this preprocessing is not stored in
 * the PDict object so it has to be done again each time matchPDict() is
 * called.
 * TODO: Estimate the cost of this preprocessing and decide whether it's
 * worth to make it persistent. Not a trivial task!
 */

static PPHeadTail new_PPHeadTail(SEXP base_codes, int bmbuf_nrow,
		int max_Hwidth, int max_Twidth, int max_mm)
{
	PPHeadTail ppheadtail;
	int i;

	ppheadtail.is_init = 1;
	if (LENGTH(base_codes) != 4)
		error("Biostrings internal error in _new_HeadTail(): "
			"LENGTH(base_codes) != 4");
	_init_byte2offset_with_INTEGER(ppheadtail.byte2offset, base_codes, 1);
	if (max_Hwidth > 0)
		for (i = 0; i < 4; i++)
			ppheadtail.head_bmbuf[i] = _new_BitMatrix(bmbuf_nrow,
							max_Hwidth, 0UL);
	if (max_Twidth > 0)
		for (i = 0; i < 4; i++)
			ppheadtail.tail_bmbuf[i] = _new_BitMatrix(bmbuf_nrow,
							max_Twidth, 0UL);
	ppheadtail.nmis_bmbuf = _new_BitMatrix(bmbuf_nrow, max_mm + 1, 0UL);
	//Rprintf("new_PPHeadTail():\n");
	//Rprintf("  nb of rows in each BitMatrix buffer=%d\n", bmbuf_nrow);
	return ppheadtail;
}

HeadTail _new_HeadTail(SEXP pdict_head, SEXP pdict_tail,
		SEXP pptb, SEXP max_mismatch, SEXP fixed)
{
	HeadTail headtail;
	int tb_length, max_mm, fixedP, fixedS,
	    key, max_Hwidth, max_Twidth, max_HTwidth, HTwidth, keys_buflength;
	SEXP low2high, dups, base_codes;
	RoSeqs head, tail;
	RoSeq *H, *T;

	tb_length = _get_PreprocessedTB_length(pptb);
	low2high = _get_PreprocessedTB_low2high(pptb);
	max_mm = INTEGER(max_mismatch)[0];
	fixedP = LOGICAL(fixed)[0];
	fixedS = LOGICAL(fixed)[1];
	if (pdict_head == R_NilValue) {
		head = _alloc_RoSeqs(tb_length);
		for (key = 0, H = head.elts; key < tb_length; key++, H++)
			H->nelt = 0;
	} else {
		head = _new_RoSeqs_from_XStringSet(tb_length, pdict_head);
	}
	if (pdict_tail == R_NilValue) {
		tail = _alloc_RoSeqs(tb_length);
		for (key = 0, T = tail.elts; key < tb_length; key++, T++)
			T->nelt = 0;
	} else {
		tail = _new_RoSeqs_from_XStringSet(tb_length, pdict_tail);
	}
	max_Hwidth = max_Twidth = max_HTwidth = keys_buflength = 0;
	for (key = 0, H = head.elts, T = tail.elts; key < tb_length; key++, H++, T++) {
		if (H->nelt > max_Hwidth)
			max_Hwidth = H->nelt;
		if (T->nelt > max_Twidth)
			max_Twidth = T->nelt;
		HTwidth = H->nelt + T->nelt;
		if (HTwidth > max_HTwidth)
			max_HTwidth = HTwidth;
		dups = VECTOR_ELT(low2high, key);
		if (dups != R_NilValue && LENGTH(dups) > keys_buflength)
			keys_buflength = LENGTH(dups);
	}
	keys_buflength++;
	headtail.head = head;
	headtail.tail = tail;
	headtail.max_Hwidth = max_Hwidth;
	headtail.max_Twidth = max_Twidth;
	headtail.max_HTwidth = max_HTwidth;
	headtail.keys_buf = new_IntAE(keys_buflength, keys_buflength, 0);
	//Rprintf("_new_HeadTail():\n");
	//Rprintf("  tb_length=%d max_mm=%d\n", tb_length, max_mm);
	//Rprintf("  max_Hwidth=%d max_Twidth=%d max_HTwidth=%d\n",
	//	max_Hwidth, max_Twidth, max_HTwidth);
	//Rprintf("  keys_buflength=%d\n", keys_buflength);
	if (max_mm < max_HTwidth && max_mm <= 4 && fixedP && fixedS) {
		/* The base codes for the head and tail are assumed to be the
		   same as for the Trusted Band */
		base_codes = _get_PreprocessedTB_base_codes(pptb);
		headtail.ppheadtail = new_PPHeadTail(base_codes, keys_buflength,
						max_Hwidth, max_Twidth, max_mm);
	} else {
		headtail.ppheadtail.is_init = 0;
	}
	return headtail;
}

static void init_headortail_bmbuf(BitMatrix *bmbuf, int nrow)
{
	int i;

	//Rprintf("init_headortail_bmbuf(): nrow=%d\n", nrow);
	for (i = 0; i < 4; i++) {
		if (nrow > bmbuf[i].nword_per_col * NBIT_PER_BITWORD)
			error("Biostrings internal error in init_headortail_bmbuf(): "
			      "not enough rows in 'bmbuf[%d]'", i);
		bmbuf[i].nrow = nrow;
		/* Set all bits to 1
		   FIXME: no need to do this for the entire buffer! */
		_BitMatrix_set_val(bmbuf + i, UINT_MAX);
	}
	return;
}

static void init_nmis_bmbuf(BitMatrix *bmbuf, int nrow)
{
	if (nrow > bmbuf->nword_per_col * NBIT_PER_BITWORD)
		error("Biostrings internal error in init_nmis_bmbuf(): "
		      "not enough rows in 'bmbuf'");
	bmbuf->nrow = nrow;
	/* Set all bits to 0
	   FIXME: no need to do this for the entire buffer! */
	_BitMatrix_set_val(bmbuf, 0UL);
	return;
}

static void preprocess_H(const RoSeq *H, ByteTrTable byte2offset,
		BitMatrix *bmbuf0, int i)
{
	int j, offset;
	const char *c;
	BitMatrix *bmbuf;

	for (j = 0, c = H->elts + H->nelt - 1; j < H->nelt; j++, c--) {
		offset = byte2offset[(unsigned char) *c];
		if (offset == NA_INTEGER)
			error("preprocess_H(): don't know how to handle "
			      "non-base letters in the preprocessed head or "
			      "tail of a PDict object yet, sorry ==> FIXME");
		bmbuf = bmbuf0 + offset;
		_BitMatrix_set_bit(bmbuf, i, j, 0);
	}
	for (offset = 0; offset < 4; offset++) {
		bmbuf = bmbuf0 + offset;
		for (j = H->nelt; j < bmbuf->ncol; j++)
			_BitMatrix_set_bit(bmbuf, i, j, 0);
	}
	return;
}

static void preprocess_T(const RoSeq *T, ByteTrTable byte2offset,
		BitMatrix *bmbuf0, int i)
{
	int j, offset;
	const char *c;
	BitMatrix *bmbuf;

	for (j = 0, c = T->elts; j < T->nelt; j++, c++) {
		offset = byte2offset[(unsigned char) *c];
		if (offset == NA_INTEGER)
			error("preprocess_T(): don't know how to handle "
			      "non-base letters in the preprocessed head or "
			      "tail of a PDict object yet, sorry ==> FIXME");
		bmbuf = bmbuf0 + offset;
		_BitMatrix_set_bit(bmbuf, i, j, 0);
	}
	for (offset = 0; offset < 4; offset++) {
		bmbuf = bmbuf0 + offset;
		for (j = T->nelt; j < bmbuf->ncol; j++)
			_BitMatrix_set_bit(bmbuf, i, j, 0);
	}
	return;
}

static void preprocess_head(const RoSeqs *head, const IntAE *keys,
		ByteTrTable byte2offset, BitMatrix *bmbuf0)
{
	int i, *key;

	init_headortail_bmbuf(bmbuf0, keys->nelt);
	for (i = 0, key = keys->elts; i < keys->nelt; i++, key++)
		preprocess_H(head->elts + *key, byte2offset, bmbuf0, i);
	return;
}

static void preprocess_tail(const RoSeqs *tail, const IntAE *keys,
		ByteTrTable byte2offset, BitMatrix *bmbuf0)
{
	int i, *key;

	init_headortail_bmbuf(bmbuf0, keys->nelt);
	for (i = 0, key = keys->elts; i < keys->nelt; i++, key++)
		preprocess_T(tail->elts + *key, byte2offset, bmbuf0, i);
	return;
}

static void match_ppheadtail_for_loc(HeadTail *headtail,
		const RoSeq *S, int tb_end, int max_mm,
		MatchPDictBuf *matchpdict_buf)
{
	const BitMatrix *head_bmbuf, *tail_bmbuf;
	int j1, j2, i, offset;
	char s;
	BitCol bitcol;

	// Match the heads
	head_bmbuf = headtail->ppheadtail.head_bmbuf;
	for (j1 = 0, j2 = tb_end - matchpdict_buf->tb_matches.tb_width - 1;
	     j1 < headtail->max_Hwidth;
	     j1++, j2--)
	{
		// 'j2' should be a safe location in 'S' because we call
		// match_ppheadtail_for_loc() only when 'tb_end' is guaranteed
		// not to be too close to 'S' boundaries.
		s = S->elts[j2];
		offset = headtail->ppheadtail.byte2offset[(unsigned char) s];
		bitcol = _BitMatrix_get_col(head_bmbuf + offset, j1);
		_BitMatrix_grow1rows(&(headtail->ppheadtail.nmis_bmbuf), &bitcol);
	}
	// Match the tails
	tail_bmbuf = headtail->ppheadtail.tail_bmbuf;
	for (j1 = 0, j2 = tb_end;
	     j1 < headtail->max_Twidth;
	     j1++, j2++)
	{
		// 'j2' should be a safe location in 'S' because we call
		// match_ppheadtail_for_loc() only when 'tb_end' is guaranteed
		// not to be too close from 'S' boundaries.
		s = S->elts[j2];
		offset = headtail->ppheadtail.byte2offset[(unsigned char) s];
		bitcol = _BitMatrix_get_col(tail_bmbuf + offset, j1);
		_BitMatrix_grow1rows(&(headtail->ppheadtail.nmis_bmbuf), &bitcol);
	}
	// Report the matches
	bitcol = _BitMatrix_get_col(&(headtail->ppheadtail.nmis_bmbuf), max_mm);
	for (i = 0; i < bitcol.nbit; i++) {
		if (!_BitCol_get_bit(&bitcol, i))
			_MatchPDictBuf_report_match(matchpdict_buf,
				headtail->keys_buf.elts[i], tb_end);
	}
	return;
}

static void match_ppheadtail(HeadTail *headtail,
		const RoSeq *S, const IntAE *tb_end_buf, int max_mm,
		MatchPDictBuf *matchpdict_buf)
{
	int min_safe_tb_end, max_safe_tb_end, j;
	const int *tb_end;

	if (headtail->max_Hwidth > 0)
		preprocess_head(&(headtail->head), &(headtail->keys_buf),
			headtail->ppheadtail.byte2offset,
			headtail->ppheadtail.head_bmbuf);
	if (headtail->max_Twidth > 0)
		preprocess_tail(&(headtail->tail), &(headtail->keys_buf),
			headtail->ppheadtail.byte2offset,
			headtail->ppheadtail.tail_bmbuf);
	min_safe_tb_end = headtail->max_Hwidth
	                + matchpdict_buf->tb_matches.tb_width;
	max_safe_tb_end = S->nelt - headtail->max_Twidth;
	for (j = 0, tb_end = tb_end_buf->elts;
	     j < tb_end_buf->nelt;
	     j++, tb_end++)
	{
		if (*tb_end < min_safe_tb_end || max_safe_tb_end < *tb_end) {
			match_headtail_for_loc(headtail,
					S, *tb_end, max_mm,
					matchpdict_buf);
			continue;
		}
		// From now 'tb_end' is guaranteed to be "safe" i.e. not too
		// close to 'S' boundaries.
		init_nmis_bmbuf(&(headtail->ppheadtail.nmis_bmbuf),
					headtail->keys_buf.nelt);
		match_ppheadtail_for_loc(headtail,
				S, *tb_end, max_mm,
				matchpdict_buf);
	}
	return;
}



/*****************************************************************************
 * _match_pdict_flanks_at() and _match_pdict_all_flanks()
 * ------------------------------------------------------
 */

void _match_pdict_flanks_at(int key0, SEXP low2high,
		HeadTail *headtail,
		const RoSeq *S, int tb_end, int max_mm, int fixedP,
		MatchPDictBuf *matchpdict_buf)
{
/*
	static ncalls = 0;

	ncalls++;
	Rprintf("_match_pdict_flanks_at(): ncalls=%d key0=%d tb_end=%d\n",
		ncalls, key0, tb_end);
*/
	collect_keys(key0, low2high, &(headtail->keys_buf));
	match_headtail_for_loc(headtail,
		S, tb_end, max_mm,
		matchpdict_buf);
	return;
}

/* If 'headtail' is empty (i.e. headtail->max_HTwidth == 0) then
   _match_pdict_all_flanks() just propagates the matches to the duplicates */
void _match_pdict_all_flanks(SEXP low2high,
		HeadTail *headtail,
		const RoSeq *S, int max_mm,
		MatchPDictBuf *matchpdict_buf)
{
	const IntAE *tb_matching_keys, *tb_end_buf;
	int i, key0;

	long unsigned int ndup, nloci, NFC; // NFC = Number of Flank Comparisons
	static long unsigned int total_NFC, subtotal_NFC;

#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] ENTERING _match_pdict_all_flanks()\n");
#endif
	tb_matching_keys = &(matchpdict_buf->tb_matches.matching_keys);
	for (i = 0; i < tb_matching_keys->nelt; i++) {
		key0 = tb_matching_keys->elts[i];
		collect_keys(key0, low2high, &(headtail->keys_buf));
		//Rprintf("_match_pdict_all_flanks(): i=%d key0=%d\n", i, key0);
		tb_end_buf = matchpdict_buf->tb_matches.match_ends.elts + key0;
		ndup = (long unsigned int) headtail->keys_buf.nelt;
		nloci = (long unsigned int) tb_end_buf->nelt;
		NFC = ndup * nloci;
		total_NFC += NFC;

		// Use brute force
		match_headtail_by_key(headtail,
			S, tb_end_buf, max_mm,
			matchpdict_buf);

/*
		if (headtail->ppheadtail.is_init
		 && headtail->keys_buf.nelt >= 24
		 && tb_end_buf->nelt >= 10) {
			// Use the BitMatrix horse-power
			//clock_t time0;
			//double dt;
			//static double total_time;
			//int j;

			subtotal_NFC += NFC;
			//Rprintf("START using the BitMatrix horse-power: "
			//	"ndup=%lu nloci=%lu subtotal_NFC=%lu max_mm=%d\n",
			//	ndup, nloci, subtotal_NFC, max_mm);
			//time0 = clock();
			//for (j = 0; j < 1000; j++)
				match_ppheadtail(headtail,
						S, tb_end_buf, max_mm,
						matchpdict_buf);
			//dt = (double) (clock() - time0) / CLOCKS_PER_SEC;
			//total_time += dt;
			//Rprintf("END using the BitMatrix horse-power: "
			//	"time=%g total_time=%g (in seconds)\n",
			//	dt, total_time);
		} else {
			// Use brute force
			match_headtail_by_key(headtail,
				S, tb_end_buf, max_mm,
				matchpdict_buf);
		}
*/

/*
BENCHMARK
=========
Uncomment below (and comment out the use of brute force above), reinstall
and use with the following code:

library(BSgenome.Dmelanogaster.UCSC.dm3)
chr3R <- unmasked(Dmelanogaster$chr3R)
chr3R_50000 <- subseq(chr3R, end=50000)


library(drosophila2probe)
dict0 <- DNAStringSet(drosophila2probe$sequence)
pdict5 <- PDict(dict0, tb.start=11, tb.end=15)

mi5 <- matchPDict(pdict5, chr3R_50000, max.mismatch=3)

		if (headtail->ppheadtail.is_init
		 && headtail->keys_buf.nelt % 64 >= 48
		 && tb_end_buf->nelt >= 120) {
			clock_t time0;
			double dt1, dt2;

			Rprintf("_match_pdict_all_flanks(): "
				"headtail->keys_buf.nelt=%d "
				"tb_end_buf->nelt=%d\n",
				headtail->keys_buf.nelt,
				tb_end_buf->nelt);

			time0 = clock();
			for (int j = 0; j < 100; j++)
				match_ppheadtail(headtail,
					S, tb_end_buf, max_mm,
					matchpdict_buf);
			dt1 = (double) (clock() - time0) / CLOCKS_PER_SEC;
			time0 = clock();
			for (int j = 0; j < 100; j++)
				match_headtail_by_key(headtail,
					S, tb_end_buf, max_mm,
					matchpdict_buf);
			dt2 = (double) (clock() - time0) / CLOCKS_PER_SEC;
			Rprintf("  --> dt1=%.3f dt2=%.3f\n", dt1, dt2);
		}
*/
	}
	//Rprintf("_match_pdict_all_flanks(): "
	//	"total_NFC=%lu subtotal_NFC=%lu ratio=%.2f\n",
	//	total_NFC, subtotal_NFC, (double) subtotal_NFC / total_NFC);
#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] LEAVING _match_pdict_all_flanks()\n");
#endif
	return;
}

