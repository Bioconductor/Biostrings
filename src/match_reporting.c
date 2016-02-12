/****************************************************************************
 *           MatchBuf manipulation and match reporting facilities           *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "Biostrings.h"
#include "IRanges_interface.h"
#include "S4Vectors_interface.h"


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
	static MatchBuf match_buf;

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
	match_buf.ms_code = ms_code;
	match_buf.PSlink_ids = new_IntAE(0, 0, 0);
	match_buf.match_counts = new_IntAE(nPSpair, nPSpair, 0);
	if (count_only) {
		/* No match_starts and match_widths buffers in that case */
		match_buf.match_starts = NULL;
		match_buf.match_widths = NULL;
	} else {
		match_buf.match_starts = new_IntAEAE(nPSpair, nPSpair);
		match_buf.match_widths = new_IntAEAE(nPSpair, nPSpair);
	}
	return match_buf;
}

void _MatchBuf_report_match(MatchBuf *match_buf,
		int PSpair_id, int start, int width)
{
	IntAE *PSlink_ids, *count_buf, *start_buf, *width_buf;

	PSlink_ids = match_buf->PSlink_ids;
	count_buf = match_buf->match_counts;
	if (count_buf->elts[PSpair_id]++ == 0)
		IntAE_insert_at(PSlink_ids,
			IntAE_get_nelt(PSlink_ids), PSpair_id);
	if (match_buf->match_starts != NULL) {
		start_buf = match_buf->match_starts->elts[PSpair_id];
		IntAE_insert_at(start_buf, IntAE_get_nelt(start_buf), start);
	}
	if (match_buf->match_widths != NULL) {
		width_buf = match_buf->match_widths->elts[PSpair_id];
		IntAE_insert_at(width_buf, IntAE_get_nelt(width_buf), width);
	}
	return;
}

void _MatchBuf_flush(MatchBuf *match_buf)
{
	int nelt, i, PSlink_id;

	nelt = IntAE_get_nelt(match_buf->PSlink_ids);
	for (i = 0; i < nelt; i++) {
		PSlink_id = match_buf->PSlink_ids->elts[i];
		match_buf->match_counts->elts[PSlink_id] = 0;
		if (match_buf->match_starts != NULL)
			IntAE_set_nelt(match_buf->match_starts->elts[PSlink_id], 0);
		if (match_buf->match_widths != NULL)
			IntAE_set_nelt(match_buf->match_widths->elts[PSlink_id], 0);
	}
	IntAE_set_nelt(match_buf->PSlink_ids, 0);
	return;
}

void _MatchBuf_append_and_flush(MatchBuf *match_buf1,
		MatchBuf *match_buf2, int view_offset)
{
	int nelt, i, PSlink_id;
	IntAE *start_buf1, *start_buf2, *width_buf1, *width_buf2;

	if (match_buf1->ms_code == MATCHES_AS_NULL
	 || match_buf2->ms_code == MATCHES_AS_NULL)
		return;
	if (IntAE_get_nelt(match_buf1->match_counts) !=
	    IntAE_get_nelt(match_buf2->match_counts)
	 || match_buf1->ms_code != match_buf2->ms_code)
		error("Biostrings internal error in "
		      "_MatchBuf_append_and_flush(): "
		      "buffers are incompatible");
	nelt = IntAE_get_nelt(match_buf2->PSlink_ids);
	for (i = 0; i < nelt; i++) {
		PSlink_id = match_buf2->PSlink_ids->elts[i];
		if (match_buf1->match_counts->elts[PSlink_id] == 0)
			IntAE_insert_at(match_buf1->PSlink_ids,
				IntAE_get_nelt(match_buf1->PSlink_ids),
				PSlink_id);
		match_buf1->match_counts->elts[PSlink_id] +=
			match_buf2->match_counts->elts[PSlink_id];
		if (match_buf1->match_starts != NULL) {
			start_buf1 = match_buf1->match_starts->elts[PSlink_id];
			start_buf2 = match_buf2->match_starts->elts[PSlink_id];
			IntAE_append_shifted_vals(start_buf1,
				start_buf2->elts, IntAE_get_nelt(start_buf2),
				view_offset);
		}
		if (match_buf1->match_widths != NULL) {
			width_buf1 = match_buf1->match_widths->elts[PSlink_id];
			width_buf2 = match_buf2->match_widths->elts[PSlink_id];
			IntAE_append(width_buf1,
				width_buf2->elts, IntAE_get_nelt(width_buf2));
		}
	}
	_MatchBuf_flush(match_buf2);
	return;
}

SEXP _MatchBuf_which_asINTEGER(const MatchBuf *match_buf)
{
	SEXP ans;
	int i;

	PROTECT(ans = new_INTEGER_from_IntAE(match_buf->PSlink_ids));
	sort_int_array(INTEGER(ans), LENGTH(ans), 0);
	for (i = 0; i < LENGTH(ans); i++)
		INTEGER(ans)[i]++;
	UNPROTECT(1);
	return ans;
}

SEXP _MatchBuf_counts_asINTEGER(const MatchBuf *match_buf)
{
	return new_INTEGER_from_IntAE(match_buf->match_counts);
}

SEXP _MatchBuf_starts_asLIST(const MatchBuf *match_buf)
{
	if (match_buf->match_starts == NULL)
		error("Biostrings internal error: _MatchBuf_starts_asLIST() "
		      "was called in the wrong context");
	return new_LIST_from_IntAEAE(match_buf->match_starts, 1);
}

static SEXP _MatchBuf_starts_toEnvir(const MatchBuf *match_buf, SEXP env)
{
	if (match_buf->match_starts == NULL)
		error("Biostrings internal error: _MatchBuf_starts_toEnvir() "
		      "was called in the wrong context");
	return IntAEAE_toEnvir(match_buf->match_starts, env, 1);
}

static SEXP _MatchBuf_widths_asLIST(const MatchBuf *match_buf)
{
	if (match_buf->match_widths == NULL)
		error("Biostrings internal error: _MatchBuf_widths_asLIST() "
		      "was called in the wrong context");
	return new_LIST_from_IntAEAE(match_buf->match_widths, 1);
}

SEXP _MatchBuf_ends_asLIST(const MatchBuf *match_buf)
{
	if (match_buf->match_starts == NULL
	 || match_buf->match_widths == NULL)
		error("Biostrings internal error: _MatchBuf_ends_asLIST() "
		      "was called in the wrong context");
	IntAEAE_sum_and_shift(match_buf->match_starts,
			      match_buf->match_widths, -1);
	return new_LIST_from_IntAEAE(match_buf->match_starts, 1);
}

static SEXP _MatchBuf_ends_toEnvir(const MatchBuf *match_buf, SEXP env)
{
	if (match_buf->match_starts == NULL
	 || match_buf->match_widths == NULL)
		error("Biostrings internal error: _MatchBuf_ends_toEnvir() "
		      "was called in the wrong context");
	IntAEAE_sum_and_shift(match_buf->match_starts,
			      match_buf->match_widths, -1);
	return IntAEAE_toEnvir(match_buf->match_starts, env, 1);
}

/*
 * Returns the result of _MatchBuf_starts_asLIST(match_buf) and
 * _MatchBuf_widths_asLIST(match_buf) as the 2 components of an ordinary list.
 */
SEXP _MatchBuf_as_Ranges(const MatchBuf *match_buf)
{
	SEXP ans, ans_elt1, ans_elt2;

	PROTECT(ans = NEW_LIST(2));
	PROTECT(ans_elt1 = _MatchBuf_starts_asLIST(match_buf));
	SET_VECTOR_ELT(ans, 0, ans_elt1);
	UNPROTECT(1);
	PROTECT(ans_elt2 = _MatchBuf_widths_asLIST(match_buf));
	SET_VECTOR_ELT(ans, 1, ans_elt2);
	UNPROTECT(2);
	return ans;
}

SEXP _MatchBuf_as_SEXP(const MatchBuf *match_buf, SEXP env)
{
	switch (match_buf->ms_code) {
	    case MATCHES_AS_NULL:
		return R_NilValue;
	    case MATCHES_AS_WHICH:
		return _MatchBuf_which_asINTEGER(match_buf);
	    case MATCHES_AS_COUNTS:
		return _MatchBuf_counts_asINTEGER(match_buf);
	    case MATCHES_AS_STARTS:
		if (env != R_NilValue)
			return _MatchBuf_starts_toEnvir(match_buf, env);
		return _MatchBuf_starts_asLIST(match_buf);
	    case MATCHES_AS_ENDS:
		if (env != R_NilValue)
			return _MatchBuf_ends_toEnvir(match_buf, env);
		return _MatchBuf_ends_asLIST(match_buf);
	    case MATCHES_AS_RANGES:
		return _MatchBuf_as_Ranges(match_buf);
	}
	error("Biostrings internal error in _MatchBuf_as_SEXP(): "
	      "unknown 'match_buf->ms_code' value %d", match_buf->ms_code);
	return R_NilValue;
}


/****************************************************************************
 * Internal match buffer instance with a simple API.
 */

static MatchBuf internal_match_buf;
int active_PSpair_id;
static int match_shift;

void _init_match_reporting(const char *ms_mode, int nPSpair)
{
	int ms_code;

	ms_code = _get_match_storing_code(ms_mode);
	internal_match_buf = _new_MatchBuf(ms_code, nPSpair);
	active_PSpair_id = 0;
	match_shift = 0;
	return;
}

void _set_active_PSpair(int PSpair_id)
{
	active_PSpair_id = PSpair_id;
	return;
}

void _set_match_shift(int shift)
{
	match_shift = shift;
}

void _report_match(int start, int width)
{
	start += match_shift;
	_MatchBuf_report_match(&internal_match_buf,
			active_PSpair_id, start, width);
	return;
}

/* Drops reported matches for all PSpairs! */
void _drop_reported_matches()
{
	_MatchBuf_flush(&internal_match_buf);
	return;
}

int _get_match_count()
{
	return internal_match_buf.match_counts->elts[active_PSpair_id];
}

SEXP _reported_matches_asSEXP()
{
	SEXP start, width, ans;

	switch (internal_match_buf.ms_code) {
	    case MATCHES_AS_NULL:
		return R_NilValue;
	    case MATCHES_AS_COUNTS:
	    case MATCHES_AS_WHICH:
		return ScalarInteger(_get_match_count());
	    case MATCHES_AS_RANGES:
		PROTECT(start = new_INTEGER_from_IntAE(
		  internal_match_buf.match_starts->elts[active_PSpair_id]));
		PROTECT(width = new_INTEGER_from_IntAE(
		  internal_match_buf.match_widths->elts[active_PSpair_id]));
		PROTECT(ans = new_IRanges("IRanges", start, width, R_NilValue));
		UNPROTECT(3);
		return ans;
	}
	error("Biostrings internal error in _reported_matches_asSEXP(): "
	      "invalid 'internal_match_buf.ms_code' value %d",
	      internal_match_buf.ms_code);
	return R_NilValue;
}

MatchBuf *_get_internal_match_buf()
{
	return &internal_match_buf;
}

