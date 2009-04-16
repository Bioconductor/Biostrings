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
static void print_ByteTrTable(const ByteTrTable byte2code)
{
	int byte, code;

	Rprintf("[DEBUG]   Byte Translation Table:\n");
	for (byte = 0; byte < BYTETRTABLE_LENGTH; byte++) {
		Rprintf("[DEBUG]     byte=%d ", byte);
		if (32 <= byte && byte < 128)
			Rprintf("['%c']", byte);
		else
			Rprintf("     ");
		Rprintf(" -> code=");
		code = byte2code[byte];
		if (code == NA_INTEGER)
			Rprintf("NA\n");
		else
			Rprintf("%d\n", code);
	}
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

static void set_byte2offset_elt(ByteTrTable byte2offset,
		int byte, int offset, int error_on_dup)
{
	int *offset_p;

	if (byte < 0 || byte >= BYTETRTABLE_LENGTH)
		error("Biostrings internal error in set_byte2offset_elt(): "
		      "invalid byte value %d", byte);
	offset_p = byte2offset + (unsigned char) byte;
	if (*offset_p == NA_INTEGER) {
		*offset_p = offset;
		return;
	}
	if (error_on_dup)
		error("Biostrings internal error in set_byte2offset_elt(): "
		      "duplicated byte value %d", byte);
	return;
}

/*
 * Values in 'bytes' must represent byte values i.e. values >= 0 and < 256.
 * The byte offsets are written to 'byte2offset'.
*/
void _init_byte2offset_with_INTEGER(ByteTrTable byte2offset, SEXP bytes, int error_on_dup)
{
	int byte, offset;

	if (LENGTH(bytes) > BYTETRTABLE_LENGTH)
		error("Biostrings internal error in _init_byte2offset_with_INTEGER(): ",
		      "LENGTH(bytes) > BYTETRTABLE_LENGTH");
	for (byte = 0; byte < BYTETRTABLE_LENGTH; byte++)
		byte2offset[byte] = NA_INTEGER;
	for (offset = 0; offset < LENGTH(bytes); offset++) {
		byte = INTEGER(bytes)[offset];
		set_byte2offset_elt(byte2offset, byte, offset, error_on_dup);
	}
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _init_byte2offset_with_INTEGER():\n");
		print_ByteTrTable(byte2offset);
	}
#endif
	return;
}

void _init_byte2offset_with_RoSeq(ByteTrTable byte2offset, const RoSeq *seq, int error_on_dup)
{
	int byte, offset;

	if (seq->nelt > BYTETRTABLE_LENGTH)
		error("Biostrings internal error in _init_byte2offset_with_RoSeq(): ",
		      "seq->nelt > BYTETRTABLE_LENGTH");
	for (byte = 0; byte < BYTETRTABLE_LENGTH; byte++)
		byte2offset[byte] = NA_INTEGER;
	for (offset = 0; offset < seq->nelt; offset++) {
		byte = seq->elts[offset];
		set_byte2offset_elt(byte2offset, byte, offset, error_on_dup);
	}
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _init_byte2offset_with_RoSeq():\n");
		print_ByteTrTable(byte2offset);
	}
#endif
	return;
}

TwobitOligoMapper _new_TwobitOligoMapper(SEXP base_codes, int oligo_width, int endianness)
{
	TwobitOligoMapper tom;

	if (LENGTH(base_codes) != 4)
		error("_new_TwobitOligoMapper(): 'base_codes' must be of length 4");
	if (oligo_width < 1 || oligo_width > 15)
		error("_new_TwobitOligoMapper(): 'oligo_width' must be >=1 and <= 15");
	_init_byte2offset_with_INTEGER(tom.eightbit2twobit, base_codes, 1);
	tom.oligo_width = oligo_width;
	tom.endianness = endianness;
	tom.nbit_in_mask = (oligo_width - 1) * 2;
	tom.twobit_mask = (1 << tom.nbit_in_mask) - 1;
	if (endianness == 1)
		tom.twobit_mask <<= 2;
	tom.nb_valid_prev_char = 0;
	tom.current_signature = 0;
	return tom;
}

void _reset_twobit_signature(TwobitOligoMapper *tom)
{
	tom->nb_valid_prev_char = 0;
	tom->current_signature = 0;
	return;
}

int _next_twobit_signature(TwobitOligoMapper *tom, const char *c)
{
	int c_as_twobit;

	c_as_twobit = tom->eightbit2twobit[(unsigned char) *c];
	if (c_as_twobit == NA_INTEGER) {
		tom->nb_valid_prev_char = 0;
		return NA_INTEGER;
	}
	tom->nb_valid_prev_char++;
	tom->current_signature &= tom->twobit_mask;
	if (tom->endianness == 1) {
		tom->current_signature >>= 2;
		c_as_twobit <<= tom->nbit_in_mask;
	} else {
		tom->current_signature <<= 2;
	}
	tom->current_signature += c_as_twobit;
	if (tom->nb_valid_prev_char < tom->oligo_width)
		return NA_INTEGER;
	return tom->current_signature;
}

