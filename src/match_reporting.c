/****************************************************************************
 * Match reporting facilities
 * --------------------------
 */
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
	error("Biostrings internal error in _get_match_storing_code(): ",
	      "\"%s\": unknown match storing mode", ms_mode);
	return -1; /* keeps gcc -Wall happy */
}

static int ms_code;
static int match_count;
static RangeAE matchbuf;
static int match_shift;

void _init_match_reporting(const char *ms_mode)
{
	ms_code = _get_match_storing_code(ms_mode);
	if (ms_code != MATCHES_AS_NULL
	 && ms_code != MATCHES_AS_COUNTS
	 && ms_code != MATCHES_AS_RANGES)
		error("Biostrings internal error in _init_match_reporting(): ",
		      "\"%s\": unsupported match storing mode", ms_mode);
	match_count = match_shift = 0;
	matchbuf = new_RangeAE(0, 0);
	return;
}

void _drop_reported_matches()
{
	match_count = 0;
	matchbuf.start.nelt = matchbuf.width.nelt = 0;
	return;
}

void _shift_match_on_reporting(int shift)
{
	match_shift = shift;
}

void _report_match(int start, int width)
{
	start += match_shift;
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _report_match(): "
			"match found at start=%d width=%d\n", start, width);
	}
#endif
	switch (ms_code) {
		case MATCHES_AS_COUNTS:
			match_count++;
		break;
		case MATCHES_AS_RANGES:
			RangeAE_insert_at(&matchbuf, matchbuf.start.nelt, start, width);
		break;
	}
	return;
}

SEXP _reported_matches_asSEXP()
{
	switch (ms_code) {
		case MATCHES_AS_COUNTS:
			return ScalarInteger(match_count);
		case MATCHES_AS_RANGES:
			return RangeAE_asIRanges(&matchbuf);
	}
	return R_NilValue;
}

