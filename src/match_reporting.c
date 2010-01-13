/****************************************************************************
 *           MatchBuf manipulation and match reporting facilities           *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"
#include "IRanges_interface.h"

static int debug = 0;

SEXP debug_match_reporting()
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

int _get_match_storing_code(const char *ms_mode)
{
	if (strcmp(ms_mode, "MATCHES_AS_NULL") == 0)
		return MATCHES_AS_NULL;
	if (strcmp(ms_mode, "MATCHES_AS_WHICH") == 0)
		return MATCHES_AS_WHICH;
	if (strcmp(ms_mode, "MATCHES_AS_COUNTS") == 0)
		return MATCHES_AS_COUNTS;
	if (strcmp(ms_mode, "MATCHES_AS_STARTS") == 0)
		return MATCHES_AS_STARTS;
	if (strcmp(ms_mode, "MATCHES_AS_ENDS") == 0)
		return MATCHES_AS_ENDS;
	if (strcmp(ms_mode, "MATCHES_AS_RANGES") == 0)
		return MATCHES_AS_RANGES;
	if (strcmp(ms_mode, "MATCHES_AS_NORMALRANGES") == 0)
		return MATCHES_AS_NORMALRANGES;
	if (strcmp(ms_mode, "MATCHES_AS_COVERAGE") == 0)
		return MATCHES_AS_COVERAGE;
	error("Biostrings internal error in _get_match_storing_code(): "
	      "\"%s\": unknown match storing mode", ms_mode);
	return -1; /* keeps gcc -Wall happy */
}


/****************************************************************************
 * MatchBuf manipulation.
 */

MatchBuf _new_MatchBuf(int ms_code, int nPSpair)
{
	int count_only;
	static MatchBuf buf;

	if (ms_code != MATCHES_AS_NULL
	 && ms_code != MATCHES_AS_WHICH
	 && ms_code != MATCHES_AS_COUNTS
	 && ms_code != MATCHES_AS_STARTS
	 && ms_code != MATCHES_AS_ENDS
	 && ms_code != MATCHES_AS_RANGES)
		error("Biostrings internal error in _new_MatchBuf(): ",
		      "%d: unsupported match storing code", ms_code);
	count_only = ms_code == MATCHES_AS_WHICH ||
		     ms_code == MATCHES_AS_COUNTS;
	buf.ms_code = ms_code;
	buf.PSlink_ids = new_IntAE(0, 0, 0);
	buf.match_counts = new_IntAE(nPSpair, nPSpair, 0);
	if (count_only) {
		/* By setting 'buflength' to -1 we indicate that these
		   buffers must not be used */
		buf.match_starts.buflength = -1;
		buf.match_widths.buflength = -1;
	} else {
		buf.match_starts = new_IntAEAE(nPSpair, nPSpair);
		buf.match_widths = new_IntAEAE(nPSpair, nPSpair);
	}
	buf.current_PSpair_id = 0;
	return buf;
}

void _MatchBuf_flush(MatchBuf *buf)
{
	int i;
	const int *PSlink_id;

	for (i = 0, PSlink_id = buf->PSlink_ids.elts;
	     i < buf->PSlink_ids.nelt;
	     i++, PSlink_id++)
	{
		buf->match_counts.elts[*PSlink_id] = 0;
		if (buf->match_starts.buflength != -1)
			buf->match_starts.elts[*PSlink_id].nelt = 0;
		if (buf->match_widths.buflength != -1)
			buf->match_widths.elts[*PSlink_id].nelt = 0;
	}
	buf->PSlink_ids.nelt = 0;
	return;
}

void _MatchBuf_report_match(MatchBuf *buf,
		int PSpair_id, int start, int width)
{
	IntAE *PSlink_ids, *count_buf, *start_buf, *width_buf;

	PSlink_ids = &(buf->PSlink_ids);
	count_buf = &(buf->match_counts);
	if (count_buf->elts[PSpair_id]++ == 0)
		IntAE_insert_at(PSlink_ids, PSlink_ids->nelt, PSpair_id);
	if (buf->match_starts.buflength != -1) {
		start_buf = buf->match_starts.elts + PSpair_id;
		IntAE_insert_at(start_buf, start_buf->nelt, start);
	}
	if (buf->match_widths.buflength != -1) {
		width_buf = buf->match_widths.elts + PSpair_id;
		IntAE_insert_at(width_buf, width_buf->nelt, width);
	}
	return;
}

SEXP _MatchBuf_which_asINTEGER(const MatchBuf *buf)
{
	SEXP ans;
	int i;

	PROTECT(ans = IntAE_asINTEGER(&(buf->PSlink_ids)));
	sort_int_array(INTEGER(ans), LENGTH(ans), 0);
	for (i = 0; i < LENGTH(ans); i++)
		INTEGER(ans)[i]++;
	UNPROTECT(1);
	return ans;
}

SEXP _MatchBuf_counts_asINTEGER(const MatchBuf *buf)
{
	return IntAE_asINTEGER(&(buf->match_counts));
}

SEXP _MatchBuf_starts_asLIST(const MatchBuf *buf)
{
	if (buf->match_starts.buflength == -1)
		error("Biostrings internal error: _MatchBuf_starts_asLIST() "
		      "was called in the wrong context");
	return IntAEAE_asLIST(&(buf->match_starts), 1);
}

static SEXP _MatchBuf_starts_toEnvir(const MatchBuf *buf, SEXP env)
{
	if (buf->match_starts.buflength == -1)
		error("Biostrings internal error: _MatchBuf_starts_toEnvir() "
		      "was called in the wrong context");
	return IntAEAE_toEnvir(&(buf->match_starts), env, 1);
}

SEXP _MatchBuf_ends_asLIST(const MatchBuf *buf)
{
	if (buf->match_starts.buflength == -1
	 || buf->match_widths.buflength == -1)
		error("Biostrings internal error: _MatchBuf_ends_asLIST() "
		      "was called in the wrong context");
	IntAEAE_sum_and_shift(&(buf->match_starts), &(buf->match_widths), -1);
	return IntAEAE_asLIST(&(buf->match_starts), 1);
}

static SEXP _MatchBuf_ends_toEnvir(const MatchBuf *buf, SEXP env)
{
	if (buf->match_starts.buflength == -1
	 || buf->match_widths.buflength == -1)
		error("Biostrings internal error: _MatchBuf_ends_toEnvir() "
		      "was called in the wrong context");
	IntAEAE_sum_and_shift(&(buf->match_starts), &(buf->match_widths), -1);
	return IntAEAE_toEnvir(&(buf->match_starts), env, 1);
}

SEXP _MatchBuf_as_MIndex(const MatchBuf *buf)
{
	error("_MatchBuf_as_MIndex(): IMPLEMENT ME!");
	return R_NilValue;
}

SEXP _MatchBuf_as_SEXP(const MatchBuf *buf, SEXP env)
{
	switch (buf->ms_code) {
	    case MATCHES_AS_NULL:
		return R_NilValue;
	    case MATCHES_AS_WHICH:
		return _MatchBuf_which_asINTEGER(buf);
	    case MATCHES_AS_COUNTS:
		return _MatchBuf_counts_asINTEGER(buf);
	    case MATCHES_AS_STARTS:
		if (env != R_NilValue)
			return _MatchBuf_starts_toEnvir(buf, env);
		return _MatchBuf_starts_asLIST(buf);
	    case MATCHES_AS_ENDS:
		if (env != R_NilValue)
			return _MatchBuf_ends_toEnvir(buf, env);
		return _MatchBuf_ends_asLIST(buf);
	    case MATCHES_AS_RANGES:
		return _MatchBuf_as_MIndex(buf);
	}
	error("Biostrings internal error in _MatchBuf_as_SEXP(): "
	      "unsupported 'ms_code' value %d", buf->ms_code);
	return R_NilValue;
}


/****************************************************************************
 * Internal match buffer instance with a simple API.
 */

static int ms_code;
static MatchBuf internal_match_buf;
static int match_shift;

/* 'nPSpair' is ignored for now. */
void _init_match_reporting(const char *ms_mode, int nPSpair)
{
	ms_code = _get_match_storing_code(ms_mode);
	if (ms_code != MATCHES_AS_NULL
	 && ms_code != MATCHES_AS_COUNTS
	 && ms_code != MATCHES_AS_RANGES)
		error("Biostrings internal error in _init_match_reporting(): "
		      "\"%s\": unsupported match storing mode", ms_mode);
	internal_match_buf = _new_MatchBuf(ms_code, 1);
	match_shift = 0;
	return;
}

void _drop_reported_matches()
{
	_MatchBuf_flush(&internal_match_buf);
	return;
}

void _shift_match_on_reporting(int shift)
{
	match_shift = shift;
}

void _report_match(int start, int width)
{
	start += match_shift;
	_MatchBuf_report_match(&internal_match_buf, 0, start, width);
	return;
}

int _get_match_count()
{
	return internal_match_buf.match_counts.elts[0];
}

SEXP _reported_matches_asSEXP()
{
	SEXP start, width, ans;

	switch (ms_code) {
	    case MATCHES_AS_COUNTS:
		return ScalarInteger(_get_match_count());
	    case MATCHES_AS_RANGES:
		PROTECT(start = IntAE_asINTEGER(internal_match_buf.match_starts.elts));
		PROTECT(width = IntAE_asINTEGER(internal_match_buf.match_widths.elts));
		PROTECT(ans = new_IRanges("IRanges", start, width, R_NilValue));
		UNPROTECT(3);
		return ans;
	}
	return R_NilValue;
}

