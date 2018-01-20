/*****************************************************************************
 Biostrings C interface: typedefs and defines
 --------------------------------------------

   The Biostrings C interface is split in 2 files:
     1. Biostrings_defines.h (this file): contains the typedefs and defines
        of the interface.
     2. Biostrings_interface.h (in this directory): contains the prototypes
        of the Biostrings C routines that are part of the interface.

   Please consult Biostrings_interface.h for how to use this interface in your
   package.

 *****************************************************************************/
#ifndef BIOSTRINGS_DEFINES_H
#define BIOSTRINGS_DEFINES_H

#include "IRanges_defines.h"
#include "XVector_defines.h"

#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <limits.h> /* for CHAR_BIT */


/*
 * A simple typedef for representing a fast translation table from bytes to
 * integer values. Always put value NA_INTEGER in the table for bytes that are
 * not mapped.
 */
#define BYTETRTABLE_LENGTH 256
typedef struct byte_tr_table {
	int byte2code[BYTETRTABLE_LENGTH];
} ByteTrTable;

typedef struct twobit_encoding_buffer {
	ByteTrTable eightbit2twobit;
	int buflength;
	int endianness;  /* move bits to the left if 0, to the right if 1 */
	int nbit_in_mask;
	int twobit_mask;
	int lastin_twobit;
	int nb_valid_prev_char;
	int current_signature;
} TwobitEncodingBuffer;


/*
 * A simple typedef for representing a "bytewise op table".
 */
typedef struct bytewise_op_table {
	unsigned char xy2val[256][256];
} BytewiseOpTable;


/*
 * A structure for holding an array of pointers to read-only non
 * null-terminated sequences of chars.
 */
typedef struct roseqs {
	Chars_holder *elts;
	int nelt;
} RoSeqs;


/*
 * *_holder structs.
 */

typedef XVectorList_holder XStringSet_holder;

typedef struct xstringset_list_holder {
        const char *classname;
        int length;
        const int *end;
        XStringSet_holder unlistData_holder;
} XStringSetList_holder;

typedef struct mindex_holder {
	const char *classname;
	int length;
	SEXP width0;
	SEXP names;
	SEXP ends;
	SEXP dups0_high2low;
	SEXP dups0_low2high;
} MIndex_holder;


/*
 * The BitCol, BitMatrix and HeadTail structs are used for preprocessing
 * and fast matching of the head and tail of a PDict object.
 */
typedef unsigned long int BitWord;

#define NBIT_PER_BITWORD (sizeof(BitWord) * CHAR_BIT)

typedef struct bitcol {
	BitWord *bitword0;
	int nword;
	int nbit; /* <= nword * NBIT_PER_BITWORD */
} BitCol;

typedef struct bitmatrix {
	BitWord *bitword00;
	int nword_per_col;
	int nrow; /* <= nword_per_col * NBIT_PER_BITWORD */
	int ncol;
} BitMatrix;

typedef struct ppheadtail {
	int is_init;
	ByteTrTable byte2offset;
	BitMatrix head_bmbuf[4], tail_bmbuf[4];
	BitMatrix nmis_bmbuf;
	BitMatrix tmp_match_bmbuf;
	int *tmp_tb_end_buf;
} PPHeadTail;

typedef struct headtail {
	RoSeqs head, tail;
	int max_Hwidth, max_Twidth, max_HTwidth;
	IntAE *grouped_keys;
	PPHeadTail ppheadtail;
} HeadTail;


/*
 * Match storing modes.
 * np = nb of pattern sequences. ns = nb of subject sequences.
 * The "key" of a pattern (or subject) sequence is its 1-based position in the
 * pattern dictionary (or in the set of subject sequences).
 *                    |                           |       np = N, ns = 1      |
 *                    |        np = ns = 1        |    or np = 1, ns = N      |
 * -------------------|-------------------------------------------------------|
 *  MATCHES_AS_NULL   |             Matches are not stored at all.            |
 *                    |                   NULL is returned.                   |
 * -------------------|-------------------------------------------------------|
 *  MATCHES_AS_WHICH  |                           | Only the matching "keys"  |
 *                    |                           | are stored and returned.  |
 * -------------------|-------------------------------------------------------|
 *  MATCHES_AS_COUNTS |                Matches are counted only.              |
 *                    |    An integer vector of length np * ns is returned.   |
 * -------------------|-------------------------------------------------------|
 *  MATCHES_AS_STARTS | Only the starts of the matches are stored.            |
 *                    | An integer vector is      | A list of integer vectors |
 *                    | returned.                 | is returned (with         |
 *                    |                           | eventually NULL elements).|
 * -------------------|-------------------------------------------------------|
 *  MATCHES_AS_ENDS   | Only the ends of the matches are stored.              |
 *                    | An integer vector is      | A list of integer vectors |
 *                    | returned.                 | is returned (with         |
 *                    |                           | eventually NULL elements).|
 * -------------------|-------------------------------------------------------|
 *  MATCHES_AS_RANGES | The starts and ends of the matches are stored.        |
 *                    | An IntegerRanges object   | An MIndex object is       |
 *                    | is returned.              | returned.                 |
 */
#define MATCHES_AS_NULL		0
#define MATCHES_AS_WHICH	1
#define MATCHES_AS_COUNTS	2
#define MATCHES_AS_STARTS	3
#define MATCHES_AS_ENDS		4
#define MATCHES_AS_RANGES	5
#define MATCHES_AS_NORMALRANGES	6  // not supported yet
#define MATCHES_AS_COVERAGE	7  // supported yet

/* The 'PSlink_ids' field contains the ids of the pattern/subject pairs that
   are linked by at least 1 match. */
typedef struct match_buf {
	int ms_code;
	IntAE *PSlink_ids;
	IntAE *match_counts;
	IntAEAE *match_starts;  /* can be missing! (i.e. set to NULL) */
	IntAEAE *match_widths;  /* can be missing! (i.e. set to NULL) */
} MatchBuf;


/*
 * The MatchPDictBuf struct is used for storing the matches found by the
 * matchPDict() function (and family).
 */
typedef struct tbmatch_buf {
	int is_init;
	int tb_width;
	const int *head_widths;
	const int *tail_widths;
	IntAE *PSlink_ids;
	IntAEAE *match_ends;
} TBMatchBuf;

typedef struct matchpdict_buf {
	TBMatchBuf tb_matches;
	MatchBuf matches;
} MatchPDictBuf;

#endif
