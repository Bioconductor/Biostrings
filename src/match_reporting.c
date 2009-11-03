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

#define	DEVNULL		0  // Matches are not stored at all -> NULL is returned.
#define	COUNTONLY	1  // Matches are counted only -> a single integer is
			   // returned.
#define	ASIRANGES	2  // Matches are fully stored -> an IRanges object is
			   // returned.

/* Not supported yet */
#define	ASNORMALIRANGES	3  // Matches are stored in an IRanges object
			   // that is kept normal at each addition.
#define	ASCOVERAGE	4  // Only the coverage of the matches is stored
			   // in an XInteger of XRleInteger object.

static int mrmode;
static int match_count;
static RangeAE matchbuf;
static int match_shift;

void _init_match_reporting(const char *mode)
{
	if (strcmp(mode, "DEVNULL") == 0)
		mrmode = DEVNULL;
	else if (strcmp(mode, "COUNTONLY") == 0)
		mrmode = COUNTONLY;
	else if (strcmp(mode, "ASIRANGES") == 0)
		mrmode = ASIRANGES;
	else
		error("\"%s\": unsupported match reporting mode", mode);
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
	switch (mrmode) {
		case COUNTONLY:
			match_count++;
		break;
		case ASIRANGES:
			RangeAE_insert_at(&matchbuf, matchbuf.start.nelt, start, width);
		break;
	}
	return;
}

SEXP _reported_matches_asSEXP()
{
	switch (mrmode) {
		case COUNTONLY:
			return ScalarInteger(match_count);
		case ASIRANGES:
			return RangeAE_asIRanges(&matchbuf);
	}
	return R_NilValue;
}

