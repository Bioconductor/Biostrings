#include "Biostrings.h"

#include <stdlib.h>
#include <ctype.h> /* for isspace() */

static int debug = 0;

SEXP debug_utils()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'utils.c'\n", debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'utils.c'\n");
#endif
	return R_NilValue;
}


/* Like fgets() except that:
 *   - the string stored into the buffer pointed to by s is right-trimmed i.e.
 *     all the rightmost white-space characters were removed,
 *   - return the length of the string stored into the buffer pointed to by s
 *     on success and -1 on error or when end of file occurs while no
 *     characters have been read.
 */
int fgets_rtrimmed(char *s, int size, FILE *stream)
{
	char *s1;
	int line_len, i;
	long pos0;

	pos0 = ftell(stream);
	s1 = fgets(s, size, stream);
	if (s1 == NULL)
		return -1;
	/* 2 almost equivalent ways to get the length of the current line,
	   "almost" because of they will differ if a line contains embedded
	   NUL characters */
	line_len = ftell(stream) - pos0; /* should be faster than strlen() */
	/* line_len = strlen(s); */
	i = line_len - 1;
	while (i >= 0 && isspace(s[i])) i--;
	line_len = i + 1;
	s[line_len] = 0;
	return line_len;
}

/*
 * Values in 'codes' must represent byte values i.e. values >= 0 and < 256.
 * Output is written to 'chrtrtable' which must be a writable int array of
 * length CHRTRTABLE_LENGTH (256).
 */
void _init_chrtrtable(const int *codes, int len, int *chrtrtable)
{
	int code, *offset_p, offset;

	for (code = 0, offset_p = chrtrtable;
	     code < CHRTRTABLE_LENGTH;
	     code++, offset_p++)
		*offset_p = -1;
	for (offset = 0; offset < len; offset++, codes++) {
		code = *codes;
		if (code < 0 || code >= CHRTRTABLE_LENGTH)
			error("Biostrings internal error in _init_chrtrtable(): "
			      "invalid code %d", code);
		chrtrtable[(unsigned char) code] = offset;
	}
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		for (code = 0, offset_p = chrtrtable; code < CHRTRTABLE_LENGTH; code++, offset_p++)
			Rprintf("[DEBUG] _init_chrtrtable(): code=%d offset=%d\n", code, *offset_p);
	}
#endif
	return;
}

