/*****************************************************************************
 Biostrings C interface: typedefs and defines
 --------------------------------------------

   The Biostrings C interface is splitted in 2 files:
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

#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <limits.h> /* for CHAR_BIT */


/*
 * A simple typedef for representing a fast translation table from bytes to
 * integer values. Always put value NA_INTEGER in the table for bytes that are
 * not mapped.
 */
#define BYTETRTABLE_LENGTH 256
typedef int ByteTrTable[BYTETRTABLE_LENGTH];

typedef struct twobit_encoding_buffer {
	ByteTrTable eightbit2twobit;
	int buflength;
	int endianness;  /* move bits to the left if 0, to the right if 1 */
	int nbit_in_mask;
	int twobit_mask;
	int nb_valid_prev_char;
	int current_signature;
} TwobitEncodingBuffer;


/*
 * Two structures for holding pointers to read-only non null-terminated
 * sequences of chars:
 *
 *   o RoSeq:  array of const chars (think of this as a pointer to a non
 *             null-terminated sequence of chars);
 *   o RoSeqs: array of arrays of const chars;
 */
typedef struct roseq {
	const char *elts;
	int nelt;
} RoSeq;

typedef struct roseqs {
	RoSeq *elts;
	int nelt;
} RoSeqs;


/*
 * Use the CachedXStringSet struct for fast extraction of the elements of
 * an XStringSet object in a loop.
 */
typedef struct cachedxstringset {
	int *start;
	int *width;
	// cannot use const char * here because of function
	// _write_RoSeq_to_CachedXStringSet_elt()
	char *super_elts;
	int super_nelt;
	const char *baseClass;
	const ByteTrTable *enc_byte2code;
	const ByteTrTable *dec_byte2code;
} CachedXStringSet;


/*
 * The BitCol and BitMatrix structs are used for preprocessing and fast
 * matching of the head and tail of a PDict object.
 */
typedef long unsigned int BitWord;

#define NBIT_PER_BITWORD (sizeof(BitWord) * CHAR_BIT)

typedef struct bitcol {
	BitWord *words;
	int nword;
	int nbit; /* <= nword * NBIT_PER_BITWORD */
} BitCol;

typedef struct bitmatrix {
	BitWord *words;
	int nword_per_col;
	int nrow; /* <= nword_per_col * NBIT_PER_BITWORD */
	int ncol;
} BitMatrix;


/*
 * The MatchPDictBuf struct is used for storing the matches found by the
 * matchPDict() function.
 */
#define MATCHES_AS_NULL		0  // Matches are not stored at all -> NULL is
				   // returned.
#define MATCHES_AS_WHICH	1  // Only the matching keys are returned.
#define MATCHES_AS_COUNTS	2  // Matches are counted only -> an integer
				   // vector is returned.
#define MATCHES_AS_STARTS	3  // Only the ends of the matches are stored
				   // -> a list of integer vectors is returned
				   // (with eventually NULL elements).
#define MATCHES_AS_ENDS		4  // Only the ends of the matches are stored
				   // -> a list of integer vectors is returned
				   // (with eventually NULL elements).
/* Not supported yet */
#define MATCHES_AS_MINDEX	5  // Matches are fully stored -> an MIndex
				   // object is returned.

typedef struct tbmatch_buf {
	int is_init;
	int tb_width;
	const int *head_widths;
	const int *tail_widths;
	IntAE matching_keys;
	IntAEAE match_ends;
} TBMatchBuf;

typedef struct seq2match_buf {
	IntAE matching_keys;
	IntAE match_counts;
	IntAEAE match_starts;
	IntAEAE match_widths;
} Seq2MatchBuf;

typedef struct matchpdict_buf {
	int matches_as;
	TBMatchBuf tb_matches;
	Seq2MatchBuf matches;
} MatchPDictBuf;

#endif
