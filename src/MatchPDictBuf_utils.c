/*****************************************************************************
 * Low-level manipulation of the MatchPDictBuf buffer
 * --------------------------------------------------
 * The MatchPDictBuf struct is used for storing the matches found by the
 * matchPDict() function.
 */
#include "Biostrings.h"
#include "IRanges_interface.h"


static int debug = 0;

SEXP debug_MatchPDictBuf_utils()
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

Seq2MatchBuf _new_Seq2MatchBuf(int nseq)
{
	static Seq2MatchBuf buf;

	buf.matching_keys = new_IntAE(0, 0, 0);
	buf.match_counts = new_IntAE(nseq, nseq, 0);
	buf.match_starts = new_IntAEAE(nseq, nseq);
	buf.match_widths = new_IntAEAE(nseq, nseq);
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
		buf->match_counts.elts[*key] =
			buf->match_starts.elts[*key].nelt =
			buf->match_widths.elts[*key].nelt = 0;
	}
	buf->matching_keys.nelt = 0;
	return;
}

SEXP _Seq2MatchBuf_which_asINTEGER(Seq2MatchBuf *buf)
{
	IntAE_shift(&(buf->matching_keys), 1);
	IntAE_qsort(&(buf->matching_keys));
	return IntAE_asINTEGER(&(buf->matching_keys));
}

SEXP _Seq2MatchBuf_counts_asINTEGER(Seq2MatchBuf *buf)
{
	return IntAE_asINTEGER(&(buf->match_counts));
}

SEXP _Seq2MatchBuf_starts_asLIST(Seq2MatchBuf *buf)
{
	return IntAEAE_asLIST(&(buf->match_starts), 1);
}

SEXP _Seq2MatchBuf_ends_asLIST(Seq2MatchBuf *buf)
{
	IntAEAE_sum_and_shift(&(buf->match_starts), &(buf->match_widths), -1);
	return IntAEAE_asLIST(&(buf->match_starts), 1);
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
			return IntAEAE_toEnvir(&(buf->match_starts), env, 1);
		return _Seq2MatchBuf_starts_asLIST(buf);
	    case MATCHES_AS_ENDS:
		if (env != R_NilValue) {
			IntAEAE_sum_and_shift(&(buf->match_starts), &(buf->match_widths), -1);
			return IntAEAE_toEnvir(&(buf->match_starts), env, 1);
		}
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
		buf.count_only = buf.matches_as == MATCHES_AS_WHICH ||
				 buf.matches_as == MATCHES_AS_COUNTS;
		buf.tb_matches = _new_TBMatchBuf(nseq, tb_width, head_widths, tail_widths);
		buf.matches = _new_Seq2MatchBuf(nseq);
	}
	return buf;
}

void _MatchPDictBuf_report_match(MatchPDictBuf *buf, int key, int tb_end)
{
	IntAE *matching_keys, *count_buf, *start_buf, *width_buf;
	int start, width;

	if (!buf->tb_matches.is_init)
		return;
	matching_keys = &(buf->matches.matching_keys);
	count_buf = &(buf->matches.match_counts);
	if (count_buf->elts[key]++ == 0)
		IntAE_insert_at(matching_keys, matching_keys->nelt, key);
	if (buf->count_only)
		return;
	start_buf = buf->matches.match_starts.elts + key;
	width_buf = buf->matches.match_widths.elts + key;
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
	IntAE_insert_at(start_buf, start_buf->nelt, start);
	IntAE_insert_at(width_buf, width_buf->nelt, width);
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
	IntAE *starts_buf1, *starts_buf2, *widths_buf1, *widths_buf2;

	if (buf2->matches_as == MATCHES_AS_NULL)
		return;
	buf2_matches = &(buf2->matches);
	for (i = 0, key = buf2_matches->matching_keys.elts;
	     i < buf2_matches->matching_keys.nelt;
	     i++, key++)
	{
		if (buf1->match_counts.elts[*key] == 0)
			IntAE_insert_at(&(buf1->matching_keys),
					buf1->matching_keys.nelt, *key);
		buf1->match_counts.elts[*key] += buf2_matches->match_counts.elts[*key];
		starts_buf1 = buf1->match_starts.elts + *key;
		starts_buf2 = buf2_matches->match_starts.elts + *key;
		IntAE_append_shifted_vals(starts_buf1,
			starts_buf2->elts, starts_buf2->nelt, view_offset);
		widths_buf1 = buf1->match_widths.elts + *key;
		widths_buf2 = buf2_matches->match_widths.elts + *key;
		IntAE_append_shifted_vals(widths_buf1,
			widths_buf2->elts, widths_buf2->nelt, view_offset);
	}
	_MatchPDictBuf_flush(buf2);
	return;
}

