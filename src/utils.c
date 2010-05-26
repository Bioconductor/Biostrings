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

void _init_byte2offset_with_cachedCharSeq(ByteTrTable byte2offset,
		const cachedCharSeq *seq, int error_on_dup)
{
	int byte, offset;

	if (seq->length > BYTETRTABLE_LENGTH)
		error("Biostrings internal error in _init_byte2offset_with_cachedCharSeq(): ",
		      "seq->length > BYTETRTABLE_LENGTH");
	for (byte = 0; byte < BYTETRTABLE_LENGTH; byte++)
		byte2offset[byte] = NA_INTEGER;
	for (offset = 0; offset < seq->length; offset++) {
		byte = seq->seq[offset];
		set_byte2offset_elt(byte2offset, byte, offset, error_on_dup);
	}
#ifdef DEBUG_BIOSTRINGS
	if (debug) {
		Rprintf("[DEBUG] _init_byte2offset_with_cachedCharSeq():\n");
		print_ByteTrTable(byte2offset);
	}
#endif
	return;
}

TwobitEncodingBuffer _new_TwobitEncodingBuffer(SEXP base_codes, int buflength, int endianness)
{
	TwobitEncodingBuffer teb;

	if (LENGTH(base_codes) != 4)
		error("_new_TwobitEncodingBuffer(): 'base_codes' must be of length 4");
	if (buflength < 1 || buflength > 15)
		error("_new_TwobitEncodingBuffer(): 'buflength' must be >= 1 and <= 15");
	_init_byte2offset_with_INTEGER(teb.eightbit2twobit, base_codes, 1);
	teb.buflength = buflength;
	teb.endianness = endianness;
	teb.nbit_in_mask = (buflength - 1) * 2;
	teb.twobit_mask = (1 << teb.nbit_in_mask) - 1;
	if (endianness == 1)
		teb.twobit_mask <<= 2;
	teb.lastin_twobit = NA_INTEGER;
	teb.nb_valid_prev_char = 0;
	teb.current_signature = 0;
	return teb;
}

void _reset_twobit_signature(TwobitEncodingBuffer *teb)
{
	teb->lastin_twobit = NA_INTEGER;
	teb->nb_valid_prev_char = 0;
	teb->current_signature = 0;
	return;
}

int _shift_twobit_signature(TwobitEncodingBuffer *teb, char c)
{
	int lastin_twobit;

	lastin_twobit = teb->lastin_twobit =
		teb->eightbit2twobit[(unsigned char) c];
	if (lastin_twobit == NA_INTEGER) {
		teb->nb_valid_prev_char = 0;
		return NA_INTEGER;
	}
	teb->nb_valid_prev_char++;
	teb->current_signature &= teb->twobit_mask;
	if (teb->endianness == 1) {
		teb->current_signature >>= 2;
		lastin_twobit <<= teb->nbit_in_mask;
	} else {
		teb->current_signature <<= 2;
	}
	teb->current_signature += lastin_twobit;
	if (teb->nb_valid_prev_char < teb->buflength)
		return NA_INTEGER;
	return teb->current_signature;
}

int _get_twobit_signature(TwobitEncodingBuffer *teb, const cachedCharSeq *seq)
{
	int i, twobit_sign;
	const char *c;

	if (seq->length != teb->buflength)
		error("_get_twobit_signature(): seq->length != teb->buflength");
	for (i = 0, c = seq->seq; i < seq->length; i++, c++)
		twobit_sign = _shift_twobit_signature(teb, *c);
	return twobit_sign;
}

/* 'at' must contain 1-based locations in 'seq'. */
int _get_twobit_signature_at(TwobitEncodingBuffer *teb, const cachedCharSeq *seq,
		const int *at, int at_length)
{
	int i, j, twobit_sign;

	if (at_length != teb->buflength)
		error("_get_twobit_signature_at(): at_length != teb->buflength");
	for (i = 0; i < at_length; i++) {
		j = at[i];
		if (j == NA_INTEGER || j < 1 || j > seq->length)
			return -1;
		twobit_sign = _shift_twobit_signature(teb, seq->seq[j - 1]);
	}
	return twobit_sign;
}

