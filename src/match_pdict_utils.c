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

static PPHeadTail new_PPHeadTail(SEXP base_codes, SEXP low2high,
		int max_Hwidth, int max_Twidth, int max_mm)
{
	PPHeadTail ppheadtail;
	int tb_length, bmbuf_nrow, key, i;
	SEXP dups;

	ppheadtail.is_init = 1;
	if (LENGTH(base_codes) != 4)
		error("Biostrings internal error in _new_HeadTail(): "
			"LENGTH(base_codes) != 4");
	_init_byte2offset_with_INTEGER(ppheadtail.byte2offset, base_codes, 1);
	tb_length = LENGTH(low2high);
	bmbuf_nrow = 0;
	for (key = 0; key < tb_length; key++) {
		dups = VECTOR_ELT(low2high, key);
		if (dups == R_NilValue)
			continue;
		if (LENGTH(dups) > bmbuf_nrow)
			bmbuf_nrow = LENGTH(dups);
	}
	bmbuf_nrow++;
	ppheadtail.row2key = new_IntAE(bmbuf_nrow, bmbuf_nrow, 0);
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
	    key, max_Hwidth, max_Twidth, max_HTwidth, HTwidth;
	RoSeqs head, tail;
	RoSeq *H, *T;
	SEXP base_codes, low2high;

	tb_length = _get_PreprocessedTB_length(pptb);
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
	max_Hwidth = max_Twidth = max_HTwidth = 0;
	for (key = 0, H = head.elts, T = tail.elts; key < tb_length; key++, H++, T++) {
		if (H->nelt > max_Hwidth)
			max_Hwidth = H->nelt;
		if (T->nelt > max_Twidth)
			max_Twidth = T->nelt;
		HTwidth = H->nelt + T->nelt;
		if (HTwidth > max_HTwidth)
			max_HTwidth = HTwidth;
	}
	headtail.head = head;
	headtail.tail = tail;
	headtail.max_Hwidth = max_Hwidth;
	headtail.max_Twidth = max_Twidth;
	headtail.max_HTwidth = max_HTwidth;
	//Rprintf("_new_HeadTail():\n");
	//Rprintf("  tb_length=%d max_mm=%d\n", tb_length, max_mm);
	//Rprintf("  max_Hwidth=%d max_Twidth=%d max_HTwidth=%d\n",
	//	max_Hwidth, max_Twidth, max_HTwidth);
	if (max_mm < max_HTwidth && max_mm <= 4 && fixedP && fixedS) {
		/* The base codes for the head and tail are assumed to be the
		   same as for the Trusted Band */
		base_codes = _get_PreprocessedTB_base_codes(pptb);
		low2high = _get_PreprocessedTB_low2high(pptb);
		headtail.ppheadtail = new_PPHeadTail(base_codes, low2high,
						max_Hwidth, max_Twidth, max_mm);
	} else {
		headtail.ppheadtail.is_init = 0;
	}
	return headtail;
}

/* Assumes that LENGTH(dups) + 1 <= row2key->buflength  (NOT CHECKED!) */
static void init_ppheadtail_row2key(int key, SEXP dups, IntAE *row2key)
{
	int i, *elt;

	row2key->elts[0] = key;
	memcpy(row2key->elts + 1, INTEGER(dups), LENGTH(dups) * sizeof(int));
	row2key->nelt = LENGTH(dups) + 1;
	for (i = 1, elt = row2key->elts; i < row2key->nelt; i++, elt++)
		*elt -= 1;
	return;
}

static void init_ppheadtail_bmbuf(BitMatrix *bmbuf, int nrow)
{
	int i;

	//Rprintf("init_ppheadtail_bmbuf(): nrow=%d\n", nrow);
	for (i = 0; i < 4; i++) {
		if (nrow > bmbuf[i].nword_per_col * NBIT_PER_BITWORD)
			error("Biostrings internal error in init_ppheadtail_bmbuf(): "
			      "not enough rows in 'bmbuf'");
		bmbuf[i].nrow = nrow;
		/* Set all bits to 1
		   FIXME: no need to do this for the entire buffer! */
		_BitMatrix_set_val(bmbuf + i, UINT_MAX);
	}
	return;
}

static void preprocess_HorT(const RoSeq *HorT, ByteTrTable byte2offset,
		BitMatrix *bmbuf, int row)
{
	int j, offset;
	const char *c;
	BitMatrix *bmbuf0;

	for (j = 0, c = HorT->elts; j < HorT->nelt; j++, c++) {
		offset = byte2offset[(unsigned char) *c];
		if (offset == NA_INTEGER)
			error("preprocess_HorT(): don't know how to handle "
			      "non-base letters in the preprocessed head or "
			      "tail of a PDict object yet, sorry ==> FIXME");
		bmbuf0 = bmbuf + offset;
		_BitMatrix_set_bit(bmbuf0, row, j, 0);
	}
	for (offset = 0; offset < 4; offset++) {
		bmbuf0 = bmbuf + offset;
		for (j = HorT->nelt; j < bmbuf0->ncol; j++)
			_BitMatrix_set_bit(bmbuf0, row, j, 0);
	}
	return;
}

static void preprocess_head(const RoSeqs *head, PPHeadTail *ppheadtail)
{
	int i, *key;

	init_ppheadtail_bmbuf(ppheadtail->head_bmbuf, ppheadtail->row2key.nelt);
	for (i = 0, key = ppheadtail->row2key.elts;
	     i < ppheadtail->row2key.nelt;
	     i++, key++)
	{
		preprocess_HorT(head->elts + *key, ppheadtail->byte2offset,
				ppheadtail->head_bmbuf, i);
	}
	return;
}

static void preprocess_tail(const RoSeqs *tail, PPHeadTail *ppheadtail)
{
	int i, *key;

	init_ppheadtail_bmbuf(ppheadtail->tail_bmbuf, ppheadtail->row2key.nelt);
	for (i = 0, key = ppheadtail->row2key.elts;
	     i < ppheadtail->row2key.nelt;
	     i++, key++)
	{
		preprocess_HorT(tail->elts + *key, ppheadtail->byte2offset,
				ppheadtail->tail_bmbuf, i);
	}
	return;
}



/*****************************************************************************
 * _match_pdict_flanks()
 * ---------------------
 */

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

static void match_HT(int key,
		const HeadTail *headtail,
		const RoSeq *S, int tb_end, int max_mm, int fixedP,
		MatchPDictBuf *matchpdict_buf)
{
	const RoSeq *H, *T;
	int HTdeltashift, nmismatch;

	H = headtail->head.elts + key;
	T = headtail->tail.elts + key;
	HTdeltashift = H->nelt + matchpdict_buf->tb_matches.tb_width;
	nmismatch = nmismatch_in_HT(H, T,
			S, tb_end - HTdeltashift, tb_end, max_mm);
	if (nmismatch <= max_mm)
		_MatchPDictBuf_report_match(matchpdict_buf, key, tb_end);
	return;
}

void _match_pdict_flanks(int key, SEXP low2high,
		const HeadTail *headtail,
		const RoSeq *S, int tb_end, int max_mm, int fixedP,
		MatchPDictBuf *matchpdict_buf)
{
	SEXP dups;
	int *dup, i;
/*
	static ncalls = 0;

	ncalls++;
	Rprintf("_match_pdict_flanks(): ncalls=%d key=%d tb_end=%d\n", ncalls, key, tb_end);
*/
	match_HT(key, headtail,
		S, tb_end, max_mm, fixedP,
		matchpdict_buf);
	dups = VECTOR_ELT(low2high, key);
	if (dups == R_NilValue)
		return;
	for (i = 0, dup = INTEGER(dups); i < LENGTH(dups); i++, dup++)
		match_HT(*dup - 1, headtail,
			S, tb_end, max_mm, fixedP,
			matchpdict_buf);
	return;
}



/*****************************************************************************
 * _match_pdict_all_flanks()
 * -------------------------
 */

static void match_dup_headtail(int key,
		const HeadTail *headtail,
		const RoSeq *S, const IntAE *tb_end_buf, int max_mm,
		MatchPDictBuf *matchpdict_buf)
{
	const RoSeq *H, *T;
	int HTdeltashift, i, Tshift, nmismatch;

	H = headtail->head.elts + key;
	T = headtail->tail.elts + key;
	HTdeltashift = H->nelt + matchpdict_buf->tb_matches.tb_width;
	for (i = 0; i < tb_end_buf->nelt; i++) {
		Tshift = tb_end_buf->elts[i];
		nmismatch = nmismatch_in_HT(H, T,
			S, Tshift - HTdeltashift, Tshift, max_mm);
		if (nmismatch <= max_mm)
			_MatchPDictBuf_report_match(matchpdict_buf, key, Tshift);
	}
	return;
}

static void match_dups_headtail(int key, SEXP low2high,
		HeadTail *headtail,
		const RoSeq *S, const IntAE *tb_end_buf, int max_mm,
		MatchPDictBuf *matchpdict_buf)
{
	SEXP dups;
	const int *dup;
	int i, dup_key;

	dups = VECTOR_ELT(low2high, key);
/*
	if (headtail->ppheadtail.is_init
	 && dups != R_NilValue && LENGTH(dups) >= 96
	 && tb_end_buf->nelt >= 30) {
		// Let's use the BitMatrix horse-power
		long unsigned int ndup, nloci;
		static long unsigned int subtotal_NFC;
		clock_t time0;
		double dt;
		static double total_time;

		init_ppheadtail_row2key(key, dups, &(headtail->ppheadtail.row2key));
		ndup = (long unsigned int) headtail->ppheadtail.row2key.nelt;
		nloci = (long unsigned int) tb_end_buf->nelt;
		subtotal_NFC += ndup * nloci;
		Rprintf("START using BitMatrix horse-power: "
			"ndup=%lu nloci=%lu subtotal_NFC=%lu max_mm=%d\n",
			ndup, nloci, subtotal_NFC, max_mm);
		time0 = clock();
		if (headtail->max_Hwidth > 0)
			preprocess_head(&(headtail->head), &(headtail->ppheadtail));
		if (headtail->max_Twidth > 0)
			preprocess_tail(&(headtail->tail), &(headtail->ppheadtail));
		dt = (clock() - time0) * 1.00 / CLOCKS_PER_SEC;
		total_time += dt;
		Rprintf("END using BitMatrix horse-power: "
			"time=%g total_time=%g (in seconds)\n",
			dt, total_time);
		//Just preprocessing for now so we can benchmark it...
		//return;
	}
*/
	match_dup_headtail(key, headtail,
			S, tb_end_buf, max_mm,
			matchpdict_buf);
	if (dups == R_NilValue)
		return;
	for (i = 0, dup = INTEGER(dups);
	     i < LENGTH(dups);
	     i++, dup++)
	{
		dup_key = *dup - 1;
		match_dup_headtail(dup_key, headtail,
				S, tb_end_buf, max_mm,
				matchpdict_buf);
	}
	return;
}

static void match_dups_headtail2(int key, SEXP low2high,
		const HeadTail *headtail,
		const RoSeq *S, const IntAE *tb_end_buf, int max_mm,
		MatchPDictBuf *matchpdict_buf)
{
	int i;

	for (i = 0; i < tb_end_buf->nelt; i++) {
		_match_pdict_flanks(key, low2high,
				headtail,
				S, tb_end_buf->elts[i], max_mm, 1,
				matchpdict_buf);
	}
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
	int nkeys, i, key;

	SEXP dups;
	long unsigned int ndup, nloci;
	static long unsigned int total_NFC;  // Total Number of Flank Comparisons

#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] ENTERING _match_pdict_all_flanks()\n");
#endif
	tb_matching_keys = &(matchpdict_buf->tb_matches.matching_keys);
	nkeys = tb_matching_keys->nelt;
	for (i = 0; i < nkeys; i++) {
		key = tb_matching_keys->elts[i];
		dups = VECTOR_ELT(low2high, key);
		ndup = (long unsigned int) (1 + (dups == R_NilValue ? 0 : LENGTH(dups)));
		tb_end_buf = matchpdict_buf->tb_matches.match_ends.elts + key;
		nloci = (long unsigned int) tb_end_buf->nelt;
		total_NFC += ndup * nloci;
		match_dups_headtail(key, low2high,
				headtail,
				S, tb_end_buf, max_mm,
				matchpdict_buf);
	}
	//Rprintf("_match_pdict_all_flanks(): total_NFC=%lu\n", total_NFC);
#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] LEAVING _match_pdict_all_flanks()\n");
#endif
	return;
}

