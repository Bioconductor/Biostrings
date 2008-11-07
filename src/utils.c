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

#ifdef DEBUG_BIOSTRINGS
void print_ByteTrTable(const ByteTrTable byte2code)
{
	int byte;

	for (byte = 0; byte < BYTETRTABLE_LENGTH; byte++)
		Rprintf("[DEBUG]   byte=%d offset=%d\n", byte, byte2code[byte]);
	return;
}
#endif

void _init_ByteTrTable_with_lkup(ByteTrTable byte2code, SEXP lkup)
{
	int byte;

	if (LENGTH(lkup) > BYTETRTABLE_LENGTH)
		error("Biostrings internal error in _init_ByteTrTable_with_lkup(): "
		      "LENGTH(lkup) > BYTETRTABLE_LENGTH");
	for (byte = 0; byte < LENGTH(lkup); byte++)
		byte2code[byte] = INTEGER(lkup)[byte];
	for ( ; byte < BYTETRTABLE_LENGTH; byte++)
		byte2code[byte] = NA_INTEGER;
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _init_ByteTrTable_with_lkup():\n");
		print_ByteTrTable(byte2code);
	}
#endif
	return;
}

SEXP _new_lkup_from_ByteTrTable(const ByteTrTable *byte2code)
{
	SEXP ans;
	int byte;

	if (byte2code == NULL)
		return R_NilValue;
	PROTECT(ans = NEW_INTEGER(BYTETRTABLE_LENGTH));
	for (byte = 0; byte < BYTETRTABLE_LENGTH; byte++)
		INTEGER(ans)[byte] = (*byte2code)[byte];
	UNPROTECT(1);
	return ans;
}

/*
 * Values in 'bytes' must represent byte values i.e. values >= 0 and < 256.
 * The byte offsets are written to 'byte2offset'.
 */
void _init_ByteTrTable_with_offsets(ByteTrTable byte2offset, const int *bytes, int nbytes)
{
	int byte, offset;

	if (nbytes > BYTETRTABLE_LENGTH)
		error("Biostrings internal error in _init_ByteTrTable_with_offsets(): ",
		      "nbytes > BYTETRTABLE_LENGTH");
	for (byte = 0; byte < BYTETRTABLE_LENGTH; byte++)
		byte2offset[byte] = NA_INTEGER;
	for (offset = 0; offset < nbytes; offset++, bytes++) {
		byte = *bytes;
		if (byte < 0 || byte >= BYTETRTABLE_LENGTH)
			error("Biostrings internal error in _init_ByteTrTable_with_offsets(): "
			      "invalid byte %d", byte);
		byte2offset[(unsigned char) byte] = offset;
	}
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _init_ByteTrTable_with_offsets():\n");
		print_ByteTrTable(byte2offset);
	}
#endif
	return;
}

