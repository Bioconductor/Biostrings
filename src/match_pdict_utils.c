/****************************************************************************
 *                                                                          *
 *       Low-level utility functions used by the matchPDict() C code        *
 *                           Author: Herve Pages                            *
 *                                                                          *
 ****************************************************************************/
#include "Biostrings.h"
#include "IRanges_interface.h"
#include <S.h> /* for Salloc() */

#include <limits.h> /* for ULONG_MAX */
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


/****************************************************************************
 * Manipulation of the MatchPDictBuf buffer.
 *
 * The MatchPDictBuf struct is used for storing the matches found by the
 * matchPDict() function (and family).
 */

TBMatchBuf _new_TBMatchBuf(int tb_length, int tb_width,
		const int *head_widths, const int *tail_widths)
{
	static TBMatchBuf buf;

	buf.is_init = 1;
	buf.tb_width = tb_width;
	buf.head_widths = head_widths;
	buf.tail_widths = tail_widths;
	buf.PSlink_ids = new_IntAE(0, 0, 0);
	buf.match_ends = new_IntAEAE(tb_length, tb_length);
	return buf;
}

void _TBMatchBuf_report_match(TBMatchBuf *buf, int PSpair_id, int end)
{
	IntAE *end_buf;
	int nelt;

	if (!buf->is_init)
		return;
	end_buf = buf->match_ends.elts + PSpair_id;
	nelt = IntAE_get_nelt(end_buf);
	if (nelt == 0)
		IntAE_insert_at(&(buf->PSlink_ids),
				IntAE_get_nelt(&(buf->PSlink_ids)), PSpair_id);
	IntAE_insert_at(end_buf, nelt, end);
	return;
}

void _TBMatchBuf_flush(TBMatchBuf *buf)
{
	int nelt, i;
	const int *PSlink_id;

	if (!buf->is_init)
		return;
	nelt = IntAE_get_nelt(&(buf->PSlink_ids));
	for (i = 0, PSlink_id = buf->PSlink_ids.elts;
	     i < nelt;
	     i++, PSlink_id++)
	{
		IntAE_set_nelt(buf->match_ends.elts + *PSlink_id, 0);
	}
	IntAE_set_nelt(&(buf->PSlink_ids), 0);
	return;
}

MatchPDictBuf _new_MatchPDictBuf(SEXP matches_as, int tb_length, int tb_width,
		const int *head_widths, const int *tail_widths)
{
	const char *ms_mode;
	int ms_code;
	static MatchPDictBuf buf;

	ms_mode = CHAR(STRING_ELT(matches_as, 0));
	ms_code = _get_match_storing_code(ms_mode);
	if (ms_code == MATCHES_AS_NULL) {
		buf.tb_matches.is_init = 0;
	} else {
		buf.tb_matches = _new_TBMatchBuf(tb_length, tb_width,
					head_widths, tail_widths);
		buf.matches = _new_MatchBuf(ms_code, tb_length);
	}
	return buf;
}

void _MatchPDictBuf_report_match(MatchPDictBuf *buf, int PSpair_id, int tb_end)
{
	IntAE *PSlink_ids, *count_buf, *start_buf, *width_buf;
	int start, width;

	if (buf->tb_matches.is_init == 0)
		return;
	PSlink_ids = &(buf->matches.PSlink_ids);
	count_buf = &(buf->matches.match_counts);
	if (count_buf->elts[PSpair_id]++ == 0)
		IntAE_insert_at(PSlink_ids,
			IntAE_get_nelt(PSlink_ids), PSpair_id);
	width = buf->tb_matches.tb_width;
	start = tb_end - width + 1;
	if (buf->tb_matches.head_widths != NULL) {
		start -= buf->tb_matches.head_widths[PSpair_id];
		width += buf->tb_matches.head_widths[PSpair_id];
	}
	if (buf->tb_matches.tail_widths != NULL)
		width += buf->tb_matches.tail_widths[PSpair_id];
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _MatchPDictBuf_report_match():\n");
		Rprintf("  PSpair_id=%d  tb_end=%d  start=%d  width=%d\n",
			PSpair_id, tb_end, start, width);
	}
#endif
	if (buf->matches.match_starts.buflength != -1) {
		start_buf = buf->matches.match_starts.elts + PSpair_id;
		IntAE_insert_at(start_buf, IntAE_get_nelt(start_buf), start);
	}
	if (buf->matches.match_widths.buflength != -1) {
		width_buf = buf->matches.match_widths.elts + PSpair_id;
		IntAE_insert_at(width_buf, IntAE_get_nelt(width_buf), width);
	}
	return;
}

static void _MatchPDictBuf_report_match2(MatchPDictBuf *buf, int PSpair_id,
		int start, int width)
{
	if (buf->tb_matches.is_init == 0)
		return;
	_MatchBuf_report_match(&(buf->matches), PSpair_id, start, width);
	return;
}

void _MatchPDictBuf_flush(MatchPDictBuf *buf)
{
	if (buf->tb_matches.is_init == 0)
		return;
	_TBMatchBuf_flush(&(buf->tb_matches));
	_MatchBuf_flush(&(buf->matches));
	return;
}

void _MatchPDictBuf_append_and_flush(MatchBuf *buf1, MatchPDictBuf *buf2,
		int view_offset)
{
	if (buf2->tb_matches.is_init == 0)
		return;
	_MatchBuf_append_and_flush(buf1, &(buf2->matches), view_offset);
	_TBMatchBuf_flush(&(buf2->tb_matches));
	return;
}



/*****************************************************************************
 * Brute force flank comparison routines
 * -------------------------------------
 */

static void collect_grouped_keys(int key0, SEXP low2high, IntAE *grouped_keys)
{
	SEXP dups;
	int nelt, i, *key;

	nelt = 1;
	IntAE_set_nelt(grouped_keys, nelt);
	if (nelt > grouped_keys->buflength)
		error("Biostrings internal error in collect_grouped_keys(): "
		      "IntAE_get_nelt(grouped_keys) > grouped_keys->buflength");
	grouped_keys->elts[0] = key0;
	dups = VECTOR_ELT(low2high, key0);
	if (dups == R_NilValue)
		return;
	nelt += LENGTH(dups);
	IntAE_set_nelt(grouped_keys, nelt);
	if (nelt > grouped_keys->buflength)
		error("Biostrings internal error in collect_grouped_keys(): "
		      "IntAE_get_nelt(grouped_keys) > grouped_keys->buflength");
	memcpy(grouped_keys->elts + 1, INTEGER(dups),
			LENGTH(dups) * sizeof(int));
	for (i = 1, key = grouped_keys->elts + 1;
	     i < nelt;
	     i++, key++)
		*key -= 1;
	return;
}

static int nmismatch_in_HT(const Chars_holder *H, const Chars_holder *T,
		const Chars_holder *S, int Hshift, int Tshift, int max_nmis)
{
	int nmis;

	nmis = _nmismatch_at_Pshift(H, S, Hshift, max_nmis, NULL);
	if (nmis > max_nmis)
		return nmis;
	max_nmis -= nmis;
	nmis += _nmismatch_at_Pshift(T, S, Tshift, max_nmis, NULL);
	return nmis;
}

static void match_HT(const Chars_holder *H, const Chars_holder *T,
		const Chars_holder *S, int tb_end, int max_nmis, int min_nmis,
		MatchPDictBuf *matchpdict_buf, int key)
{
	int HTdeltashift, nmis;

	HTdeltashift = H->length + matchpdict_buf->tb_matches.tb_width;
	nmis = nmismatch_in_HT(H, T,
			S, tb_end - HTdeltashift, tb_end, max_nmis);
	if (nmis <= max_nmis && nmis >= min_nmis)
		_MatchPDictBuf_report_match(matchpdict_buf, key, tb_end);
	return;
}

static void match_headtail_for_loc(const HeadTail *headtail,
		const Chars_holder *S, int tb_end, int max_nmis, int min_nmis,
		MatchPDictBuf *matchpdict_buf)
{
	int nelt, i;
	const int *key;
	const Chars_holder *H, *T;

	nelt = IntAE_get_nelt(&(headtail->grouped_keys));
	for (i = 0, key = headtail->grouped_keys.elts;
	     i < nelt;
	     i++, key++)
	{
		H = headtail->head.elts + *key;
		T = headtail->tail.elts + *key;
		match_HT(H, T, S, tb_end, max_nmis, min_nmis,
				matchpdict_buf, *key);
	}
	return;
}

static void match_headtail_for_key(const HeadTail *headtail, int key,
		const Chars_holder *S, const IntAE *tb_end_buf,
		int max_nmis, int min_nmis,
		MatchPDictBuf *matchpdict_buf)
{
	const Chars_holder *H, *T;
	int nelt, j;
	const int *tb_end;

	H = headtail->head.elts + key;
	T = headtail->tail.elts + key;
	nelt = IntAE_get_nelt(tb_end_buf);
	for (j = 0, tb_end = tb_end_buf->elts;
	     j < nelt;
	     j++, tb_end++)
	{
		match_HT(H, T, S, *tb_end, max_nmis, min_nmis,
				matchpdict_buf, key);
	}
	return;
}

static void match_headtail_by_loc(const HeadTail *headtail,
		const Chars_holder *S, const IntAE *tb_end_buf,
		int max_nmis, int min_nmis,
		MatchPDictBuf *matchpdict_buf)
{
	int nelt, j;
	const int *tb_end;

	nelt = IntAE_get_nelt(tb_end_buf);
	for (j = 0, tb_end = tb_end_buf->elts;
	     j < nelt;
	     j++, tb_end++)
	{
		match_headtail_for_loc(headtail,
				S, *tb_end, max_nmis, min_nmis,
				matchpdict_buf);
	}
	return;
}

static void match_headtail_by_key(HeadTail *headtail,
		const Chars_holder *S, const IntAE *tb_end_buf,
		int max_nmis, int min_nmis,
		MatchPDictBuf *matchpdict_buf)
{
	int nelt, i;
	const int *key;

	nelt = IntAE_get_nelt(&(headtail->grouped_keys));
	for (i = 0, key = headtail->grouped_keys.elts;
	     i < nelt;
	     i++, key++)
	{
		match_headtail_for_key(headtail, *key,
				S, tb_end_buf, max_nmis, min_nmis,
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

#define MAX_REMAINING_KEYS 24  // >= 0 and < NBIT_PER_BITWORD
#define TMPMATCH_BMBUF_MAXNCOL 200

static PPHeadTail new_PPHeadTail(SEXP base_codes, int bmbuf_nrow,
		int max_Hwidth, int max_Twidth, int max_nmis)
{
	PPHeadTail ppheadtail;
	int i;

	ppheadtail.is_init = 1;
	if (LENGTH(base_codes) != 4)
		error("Biostrings internal error in _new_HeadTail(): "
			"LENGTH(base_codes) != 4");
	_init_byte2offset_with_INTEGER(&(ppheadtail.byte2offset),
				       base_codes, 1);
	if (max_Hwidth > 0)
		for (i = 0; i < 4; i++)
			ppheadtail.head_bmbuf[i] = _new_BitMatrix(bmbuf_nrow,
							max_Hwidth, 0UL);
	if (max_Twidth > 0)
		for (i = 0; i < 4; i++)
			ppheadtail.tail_bmbuf[i] = _new_BitMatrix(bmbuf_nrow,
							max_Twidth, 0UL);
	ppheadtail.nmis_bmbuf = _new_BitMatrix(bmbuf_nrow, max_nmis + 1, 0UL);
	ppheadtail.tmp_match_bmbuf = _new_BitMatrix(bmbuf_nrow, TMPMATCH_BMBUF_MAXNCOL, ULONG_MAX);
	ppheadtail.tmp_tb_end_buf = Salloc((long) TMPMATCH_BMBUF_MAXNCOL, int);
	//Rprintf("new_PPHeadTail():\n");
	//Rprintf("  nb of rows in each BitMatrix buffer=%d\n", bmbuf_nrow);
	return ppheadtail;
}

HeadTail _new_HeadTail(SEXP pdict_head, SEXP pdict_tail, SEXP pptb,
		SEXP max_mismatch, SEXP fixed,
		int with_ppheadtail)
{
	HeadTail headtail;
	int tb_length, max_nmis, fixedP, fixedS,
	    key, max_Hwidth, max_Twidth, max_HTwidth, HTwidth,
	    grouped_keys_buflength;
	SEXP low2high, dups, base_codes;
	RoSeqs head, tail;
	Chars_holder *H, *T;

	tb_length = _get_PreprocessedTB_length(pptb);
	low2high = _get_PreprocessedTB_low2high(pptb);
	max_nmis = INTEGER(max_mismatch)[0];
	fixedP = LOGICAL(fixed)[0];
	fixedS = LOGICAL(fixed)[1];
	if (pdict_head == R_NilValue) {
		head = _alloc_RoSeqs(tb_length);
		for (key = 0, H = head.elts; key < tb_length; key++, H++)
			H->length = 0;
	} else {
		head = _new_RoSeqs_from_XStringSet(tb_length, pdict_head);
	}
	if (pdict_tail == R_NilValue) {
		tail = _alloc_RoSeqs(tb_length);
		for (key = 0, T = tail.elts; key < tb_length; key++, T++)
			T->length = 0;
	} else {
		tail = _new_RoSeqs_from_XStringSet(tb_length, pdict_tail);
	}
	max_Hwidth = max_Twidth = max_HTwidth = grouped_keys_buflength = 0;
	for (key = 0, H = head.elts, T = tail.elts;
	     key < tb_length;
	     key++, H++, T++)
	{
		if (H->length > max_Hwidth)
			max_Hwidth = H->length;
		if (T->length > max_Twidth)
			max_Twidth = T->length;
		HTwidth = H->length + T->length;
		if (HTwidth > max_HTwidth)
			max_HTwidth = HTwidth;
		dups = VECTOR_ELT(low2high, key);
		if (dups != R_NilValue && LENGTH(dups) > grouped_keys_buflength)
			grouped_keys_buflength = LENGTH(dups);
	}
	grouped_keys_buflength++;
	headtail.head = head;
	headtail.tail = tail;
	headtail.max_Hwidth = max_Hwidth;
	headtail.max_Twidth = max_Twidth;
	headtail.max_HTwidth = max_HTwidth;
	headtail.grouped_keys = new_IntAE(grouped_keys_buflength, grouped_keys_buflength, 0);
	//Rprintf("_new_HeadTail():\n");
	//Rprintf("  tb_length=%d max_nmis=%d\n", tb_length, max_nmis);
	//Rprintf("  max_Hwidth=%d max_Twidth=%d max_HTwidth=%d\n",
	//	max_Hwidth, max_Twidth, max_HTwidth);
	//Rprintf("  grouped_keys_buflength=%d\n", grouped_keys_buflength);

	/* The (max_nmis <= 4) and (max_Hwidth + max_Twidth <= 10 + 4 * max_nmis)
	   criteria together with the MAX_REMAINING_KEYS value above (20)
	   are optimized for the Core 2 Duo arch (64bit) */
	if (with_ppheadtail
	 && (max_nmis < max_HTwidth)
	 && (max_nmis <= 4)
	 && (max_Hwidth + max_Twidth <= 10 + 4 * max_nmis)
	 && (fixedP && fixedS)) {
		/* The base codes for the head and tail are assumed to be the
		   same as for the Trusted Band */
		base_codes = _get_PreprocessedTB_base_codes(pptb);
		headtail.ppheadtail = new_PPHeadTail(base_codes,
					grouped_keys_buflength,
					max_Hwidth, max_Twidth, max_nmis);
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
		/* Set all bits to 1 */
		_BitMatrix_set_val(bmbuf + i, ULONG_MAX);
	}
	return;
}

static void init_nmis_bmbuf(BitMatrix *bmbuf, int nrow)
{
	if (nrow > bmbuf->nword_per_col * NBIT_PER_BITWORD)
		error("Biostrings internal error in init_nmis_bmbuf(): "
		      "not enough rows in 'bmbuf'");
	bmbuf->nrow = nrow;
	/* Set all bits to 0 */
	_BitMatrix_set_val(bmbuf, 0UL);
	return;
}

static void preprocess_H(const Chars_holder *H,
		const ByteTrTable *byte2offset, BitMatrix *bmbuf0, int i)
{
	int j, offset;
	const char *c;
	BitMatrix *bmbuf;

	for (j = 0, c = H->seq + H->length - 1; j < H->length; j++, c--) {
		offset = byte2offset->byte2code[(unsigned char) *c];
		if (offset == NA_INTEGER)
			error("preprocess_H(): don't know how to handle "
			      "non-base letters in the preprocessed head or "
			      "tail of a PDict object yet, sorry ==> FIXME");
		bmbuf = bmbuf0 + offset;
		_BitMatrix_set_bit(bmbuf, i, j, 0);
	}
	for (offset = 0; offset < 4; offset++) {
		bmbuf = bmbuf0 + offset;
		for (j = H->length; j < bmbuf->ncol; j++)
			_BitMatrix_set_bit(bmbuf, i, j, 0);
	}
	return;
}

static void preprocess_T(const Chars_holder *T,
		const ByteTrTable *byte2offset, BitMatrix *bmbuf0, int i)
{
	int j, offset;
	const char *c;
	BitMatrix *bmbuf;

	for (j = 0, c = T->seq; j < T->length; j++, c++) {
		offset = byte2offset->byte2code[(unsigned char) *c];
		if (offset == NA_INTEGER)
			error("preprocess_T(): don't know how to handle "
			      "non-base letters in the preprocessed head or "
			      "tail of a PDict object yet, sorry ==> FIXME");
		bmbuf = bmbuf0 + offset;
		_BitMatrix_set_bit(bmbuf, i, j, 0);
	}
	for (offset = 0; offset < 4; offset++) {
		bmbuf = bmbuf0 + offset;
		for (j = T->length; j < bmbuf->ncol; j++)
			_BitMatrix_set_bit(bmbuf, i, j, 0);
	}
	return;
}

static void preprocess_head(const RoSeqs *head, const IntAE *grouped_keys,
		const ByteTrTable *byte2offset, BitMatrix *bmbuf0)
{
	int nelt, i, *key;

	nelt = IntAE_get_nelt(grouped_keys);
	init_headortail_bmbuf(bmbuf0, nelt);
	for (i = 0, key = grouped_keys->elts;
	     i < nelt;
	     i++, key++)
		preprocess_H(head->elts + *key, byte2offset, bmbuf0, i);
	return;
}

static void preprocess_tail(const RoSeqs *tail, const IntAE *grouped_keys,
		const ByteTrTable *byte2offset, BitMatrix *bmbuf0)
{
	int nelt, i, *key;

	nelt = IntAE_get_nelt(grouped_keys);
	init_headortail_bmbuf(bmbuf0, nelt);
	for (i = 0, key = grouped_keys->elts;
	     i < nelt;
	     i++, key++)
		preprocess_T(tail->elts + *key, byte2offset, bmbuf0, i);
	return;
}

static BitCol match_ppheadtail_for_loc(HeadTail *headtail, int tb_width,
		const Chars_holder *S, int tb_end, int max_nmis, int min_nmis)
{
	BitMatrix *nmis_bmbuf;
	const BitMatrix *head_bmbuf, *tail_bmbuf;
	int j1, j2, offset;
	char s;
	BitCol bitcol, max_nmis_bitcol, min_nmis_bitcol;

	nmis_bmbuf = &(headtail->ppheadtail.nmis_bmbuf);
	// Match the heads
	head_bmbuf = headtail->ppheadtail.head_bmbuf;
	for (j1 = 0, j2 = tb_end - tb_width - 1;
	     j1 < headtail->max_Hwidth;
	     j1++, j2--)
	{
		// 'j2' should be a safe location in 'S' because we call
		// match_ppheadtail_for_loc() only when 'tb_end' is guaranteed
		// not to be too close to 'S' boundaries.
		s = S->seq[j2];
		offset = headtail->ppheadtail.byte2offset.byte2code[(unsigned char) s];
		if (offset == NA_INTEGER) {
			_BitMatrix_Rrot1(nmis_bmbuf);
			continue;
		}
		bitcol = _BitMatrix_get_col(head_bmbuf + offset, j1);
		_BitMatrix_grow1rows(nmis_bmbuf, &bitcol);
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
		s = S->seq[j2];
		offset = headtail->ppheadtail.byte2offset.byte2code[(unsigned char) s];
		if (offset == NA_INTEGER) {
			_BitMatrix_Rrot1(nmis_bmbuf);
			continue;
		}
		bitcol = _BitMatrix_get_col(tail_bmbuf + offset, j1);
		_BitMatrix_grow1rows(nmis_bmbuf, &bitcol);
	}
	max_nmis_bitcol = _BitMatrix_get_col(nmis_bmbuf, max_nmis);
	if (min_nmis >= 1) {
		min_nmis_bitcol = _BitMatrix_get_col(nmis_bmbuf, min_nmis - 1);
		_BitCol_A_gets_BimpliesA(&max_nmis_bitcol, &min_nmis_bitcol);
	}
	return max_nmis_bitcol;
}

static void report_matches_for_loc(const BitCol *bitcol, HeadTail *headtail,
		int tb_end, MatchPDictBuf *matchpdict_buf)
{
	// Note that using _BitCol_get_bit() for this would be easier but is
	// also twice slower!
	BitWord *bitword;
	int i, i2;

	bitword = bitcol->bitword0;
	for (i = i2 = 0; i < bitcol->nbit; i++, i2++) {
		if (i2 >= NBIT_PER_BITWORD) {
			i2 = 0;
			bitword++;
		}
		if (!(*bitword & 1UL)) {
			//_MatchPDictBuf_report_match(matchpdict_buf,
			//	headtail->grouped_keys.elts[i], tb_end);
			int key, start, width;
			key = headtail->grouped_keys.elts[i];
			width = headtail->head.elts[key].length
			      + matchpdict_buf->tb_matches.tb_width
			      + headtail->tail.elts[key].length;
			start = tb_end + headtail->tail.elts[key].length - width + 1;
			_MatchPDictBuf_report_match2(matchpdict_buf, key, start, width);
		}
		*bitword >>= 1;
	}
	return;
}

/*
static void flush_tmp_match_bmbuf(HeadTail *headtail, MatchPDictBuf *matchpdict_buf)
{
	BitMatrix *tmp_match_bmbuf;
	const int *tmp_tb_end_buf;
	BitWord *bitword2, *bitword, mask;
	int i, j;

	tmp_match_bmbuf = &(headtail->ppheadtail.tmp_match_bmbuf);
	tmp_tb_end_buf = headtail->ppheadtail.tmp_tb_end_buf;
	bitword2 = tmp_match_bmbuf->bitword00;
	mask = 1UL;
	for (i = 0; i < tmp_match_bmbuf->nrow; i++) {
		if (mask == 0UL) { // this means that i % NBIT_PER_BITWORD == 0
			bitword2++;
			mask = 1UL;
		}
		bitword = bitword2;
		for (j = 0; j < tmp_match_bmbuf->ncol; j++) {
			if (!(*bitword & mask)) {
				_MatchPDictBuf_report_match(matchpdict_buf,
					headtail->grouped_keys.elts[i], tmp_tb_end_buf[j]);
			}
			bitword += tmp_match_bmbuf->nword_per_col;
		}
		mask <<= 1;
	}
	tmp_match_bmbuf->ncol = 0;
	return;
}
*/

static void match_ppheadtail0(HeadTail *headtail,
		const Chars_holder *S, const IntAE *tb_end_buf,
		int max_nmis, int min_nmis,
		MatchPDictBuf *matchpdict_buf)
{
	BitMatrix *tmp_match_bmbuf;
	int nelt, min_safe_tb_end, max_safe_tb_end, j, ncol;
	const int *tb_end;
	BitCol bitcol;

	if (headtail->max_Hwidth > 0)
		preprocess_head(&(headtail->head), &(headtail->grouped_keys),
			&(headtail->ppheadtail.byte2offset),
			headtail->ppheadtail.head_bmbuf);
	if (headtail->max_Twidth > 0)
		preprocess_tail(&(headtail->tail), &(headtail->grouped_keys),
			&(headtail->ppheadtail.byte2offset),
			headtail->ppheadtail.tail_bmbuf);
	tmp_match_bmbuf = &(headtail->ppheadtail.tmp_match_bmbuf);
	tmp_match_bmbuf->nrow = IntAE_get_nelt(&(headtail->grouped_keys));
	tmp_match_bmbuf->ncol = 0;

	min_safe_tb_end = headtail->max_Hwidth
			+ matchpdict_buf->tb_matches.tb_width;
	max_safe_tb_end = S->length - headtail->max_Twidth;
	nelt = IntAE_get_nelt(tb_end_buf);
	for (j = 0, tb_end = tb_end_buf->elts;
	     j < nelt;
	     j++, tb_end++)
	{
		if (*tb_end < min_safe_tb_end || max_safe_tb_end < *tb_end) {
			match_headtail_for_loc(headtail,
					S, *tb_end, max_nmis, min_nmis,
					matchpdict_buf);
			continue;
		}
		// From now 'tb_end' is guaranteed to be "safe" i.e. not too
		// close to 'S' boundaries.
		init_nmis_bmbuf(&(headtail->ppheadtail.nmis_bmbuf),
				IntAE_get_nelt(&(headtail->grouped_keys)));
		bitcol = match_ppheadtail_for_loc(headtail,
				matchpdict_buf->tb_matches.tb_width,
				S, *tb_end, max_nmis, min_nmis);
		report_matches_for_loc(&bitcol, headtail, *tb_end, matchpdict_buf);
/*
		ncol = tmp_match_bmbuf->ncol;
		_BitMatrix_set_col(tmp_match_bmbuf, ncol, &bitcol);
		tmp_match_bmbuf->ncol++;
		headtail->ppheadtail.tmp_tb_end_buf[ncol] = *tb_end;
		if (tmp_match_bmbuf->ncol == TMPMATCH_BMBUF_MAXNCOL)
			flush_tmp_match_bmbuf(headtail, matchpdict_buf);
*/
	}
/*
	flush_tmp_match_bmbuf(headtail, matchpdict_buf);
*/
	return;
}

static void match_ppheadtail(HeadTail *headtail,
		const Chars_holder *S, const IntAE *tb_end_buf,
		int max_nmis, int min_nmis,
		MatchPDictBuf *matchpdict_buf)
{
	int nelt, nkey0, nkey1, nkey2, i;
	const int *key;

	nkey0 = IntAE_get_nelt(&(headtail->grouped_keys));
	nkey2 = nkey0 % NBIT_PER_BITWORD;
	if (nkey2 > MAX_REMAINING_KEYS) {
		match_ppheadtail0(headtail,
			S, tb_end_buf, max_nmis, min_nmis,
			matchpdict_buf);
		return;
	}
	nkey1 = nkey0 - nkey2;
	if (nkey1 != 0) {
		IntAE_set_nelt(&(headtail->grouped_keys), nkey1);
		match_ppheadtail0(headtail,
			S, tb_end_buf, max_nmis, min_nmis,
			matchpdict_buf);
		IntAE_set_nelt(&(headtail->grouped_keys), nkey0);
	}
	/* FIXME: If headtail->grouped_keys is guaranteed to not grow
	   during the loop below, then extract its nelt before entering
	   the loop. */
	for (i = nkey1, key = headtail->grouped_keys.elts + nkey1;
	     i < IntAE_get_nelt(&(headtail->grouped_keys));
	     i++, key++)
	{
		match_headtail_for_key(headtail, *key,
				S, tb_end_buf, max_nmis, min_nmis,
				matchpdict_buf);
	}
	return;
}

/*
static void BENCHMARK_match_ppheadtail(HeadTail *headtail,
		const Chars_holder *S, const IntAE *tb_end_buf,
		int max_nmis, int min_nmis,
		MatchPDictBuf *matchpdict_buf)
{
	clock_t time0;
	double dt1, dt2;
	int i;

	time0 = clock();
	for (i = 0; i < 100; i++) {
		match_ppheadtail0(headtail,
			S, tb_end_buf, max_nmis, min_nmis,
			matchpdict_buf);
	}
	dt1 = (double) (clock() - time0) / CLOCKS_PER_SEC;
	time0 = clock();
	for (i = 0; i < 100; i++) {
		match_headtail_by_key(headtail,
			S, tb_end_buf, max_nmis, min_nmis,
			matchpdict_buf);
	}
	dt2 = (double) (clock() - time0) / CLOCKS_PER_SEC;
	Rprintf("%.3f\t%.3f\n", dt1, dt2);
	return;
}
*/


/*****************************************************************************
 * _match_pdict_flanks_at() and _match_pdict_all_flanks()
 * ------------------------------------------------------
 */

void _match_pdict_flanks_at(int key0, SEXP low2high,
		HeadTail *headtail,
		const Chars_holder *S, int tb_end,
		int max_nmis, int min_nmis, int fixedP,
		MatchPDictBuf *matchpdict_buf)
{
/*
	static ncalls = 0;

	ncalls++;
	Rprintf("_match_pdict_flanks_at(): ncalls=%d key0=%d tb_end=%d\n",
		ncalls, key0, tb_end);
*/
	collect_grouped_keys(key0, low2high, &(headtail->grouped_keys));
	match_headtail_for_loc(headtail,
		S, tb_end, max_nmis, min_nmis,
		matchpdict_buf);
	return;
}

/* If 'headtail' is empty (i.e. headtail->max_HTwidth == 0) then
   _match_pdict_all_flanks() just propagates the matches to the duplicates */
void _match_pdict_all_flanks(SEXP low2high,
		HeadTail *headtail,
		const Chars_holder *S,
		int max_nmis, int min_nmis,
		MatchPDictBuf *matchpdict_buf)
{
	const IntAE *tb_PSlink_ids, *tb_end_buf;
	int nelt, i, key0;

	unsigned long int ndup, nloci, NFC; // NFC = Number of Flank Comparisons
	static unsigned long int total_NFC = 0UL, subtotal_NFC = 0UL;

#ifdef DEBUG_BIOSTRINGS
	if (debug)
		Rprintf("[DEBUG] ENTERING _match_pdict_all_flanks()\n");
#endif
	tb_PSlink_ids = &(matchpdict_buf->tb_matches.PSlink_ids);
	nelt = IntAE_get_nelt(tb_PSlink_ids);
	for (i = 0; i < nelt; i++) {
		key0 = tb_PSlink_ids->elts[i];
		collect_grouped_keys(key0, low2high, &(headtail->grouped_keys));
		tb_end_buf = matchpdict_buf->tb_matches.match_ends.elts + key0;
/*
		ndup = (unsigned long int) IntAE_get_nelt(&(headtail->grouped_keys));
		nloci = (unsigned long int) IntAE_get_nelt(tb_end_buf);
		NFC = ndup * nloci;
		total_NFC += NFC;
*/
		if (headtail->ppheadtail.is_init
		 && IntAE_get_nelt(tb_end_buf) >= 15) {
			// Use the BitMatrix horse-power
/*
			Rprintf("_match_pdict_all_flanks(): "
				"key0=%d "
				"IntAE_get_nelt(&(headtail->grouped_keys))=%d "
				"tb_end_buf->nelt=%d\n",
				key0,
				IntAE_get_nelt(&(headtail->grouped_keys)),
				tb_end_buf->nelt);
			subtotal_NFC += NFC;
*/
			match_ppheadtail(headtail, S, tb_end_buf,
					max_nmis, min_nmis,
					matchpdict_buf);
/*
			BENCHMARK_match_ppheadtail(headtail, S, tb_end_buf,
					max_nmis, min_nmis,
					matchpdict_buf);
*/
		} else {
			// Use brute force
			match_headtail_by_key(headtail, S, tb_end_buf,
					max_nmis, min_nmis,
					matchpdict_buf);
		}

/*
BENCHMARK
=========
Uncomment below (and comment out the use of brute force above), reinstall
and use with the following code:

library(BSgenome.Dmelanogaster.UCSC.dm3)
chr3R <- unmasked(Dmelanogaster$chr3R)
chr3R_50000 <- subseq(chr3R, end=50000)

library(drosophila2probe)
dict0 <- DNAStringSet(drosophila2probe$sequence, end=12)
pdict6 <- PDict(dict0, tb.start=4, tb.end=9)

mi6 <- matchPDict(pdict6, chr3R_50000, max.mismatch=1)
sum(countIndex(mi6))  # 17896

		if (headtail->ppheadtail.is_init
		 && IntAE_get_nelt(tb_end_buf) >= 40) {
			clock_t time0;
			double dt1, dt2;
			int j;

			time0 = clock();
			for (j = 0; j < 5000; j++)
				match_ppheadtail(headtail, S, tb_end_buf,
					max_nmis, min_nmis,
					matchpdict_buf);
			dt1 = (double) (clock() - time0) / CLOCKS_PER_SEC;
			time0 = clock();
			for (j = 0; j < 5000; j++)
				match_headtail_by_key(headtail, S, tb_end_buf,
					max_nmis, min_nmis,
					matchpdict_buf);
			dt2 = (double) (clock() - time0) / CLOCKS_PER_SEC;
			if (dt1 > dt2) {
				Rprintf("_match_pdict_all_flanks(): "
					"IntAE_get_nelt(&(headtail->grouped_keys))=%d "
					"tb_end_buf->nelt=%d\n",
					IntAE_get_nelt(&(headtail->grouped_keys)),
					IntAE_get_nelt(tb_end_buf));
				Rprintf("  --> dt1=%.3f dt2=%.3f\n", dt1, dt2);
			}
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

