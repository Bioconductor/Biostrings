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
	char *super_elts;
	int super_nelt;
	const char *baseClass;
	const int *enc_chrtrtable;
	const int *dec_chrtrtable;
} CachedXStringSet;


/*
 * Match reporting modes (more modes will be added soon...)
 */
#define COUNT_MRMODE	1
#define START_MRMODE	2


#endif
