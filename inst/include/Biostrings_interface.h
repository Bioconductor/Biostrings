/*****************************************************************************
 Biostrings C interface: prototypes
 ----------------------------------

   The Biostrings C interface is splitted in 2 files:
     1. Biostrings_defines.h (in this directory): contains the typedefs and
        defines of the interface.
     2. Biostrings_interface.h (this file): contains the prototypes of the
        Biostrings C routines that are part of the interface.

 -----------------------------------------------------------------------------

   This file contains the prototypes of a subset of Biostrings C routines that
   can be called by C code in other packages.
   In order to use these routines in your package, you need to do the
   following:

   a. Add the Biostrings package to the 'Depends:', 'Imports:' and 'LinkingTo:'
      fields of your DESCRIPTION file. Note that most of the times you should
      already have an 'Imports:' field with at least the methods package
      listed in it (this is safe to have even if you are not sure that you
      need it).

   b. If it doesn't already have one, add a NAMESPACE file to your package and
      start it with the following lines:

        useDynLib(mypackage) # replace mypackage by the name of your package
        import(methods)      # only if you have it in your 'Imports:' field
        import(Biostrings)

   c. Add a Biostrings_stubs.c file to your src/ folder with the following
      line (and nothing else) in it:

        #include "_Biostrings_stubs.c" // note the underscore!

   d. Add the following line:

        #include "Biostrings_interface.h"

      in each .c file where you need to call a Biostrings C routine.
      Note that this include should come after any R or system include so the
      top of your .c file will typically look like this:

        // R includes
        #include <Rdefines.h>
        #include <R_ext/Rdynload.h>
        ... maybe more R includes ...

        // System includes
        #include <stdio.h>
        #include <stdlib.h>
        ... maybe more system includes ...

        #include "Biostrings_interface.h"

   e. Any Biostrings C routine defined in Biostrings_interface.h can now be
      called in your .c file. For example, you could write the following
      function:

        SEXP print_XString_bytes(SEXP xstring)
        {
            RoSeq x;
            int i;
            const char *x_char;

            x = get_XString_asRoSeq(xstring);
            for (i = 0, x_char = x.elts; i < x.nelt; i++, x_char++)
                Rprintf("%x ", *x_char);
            Rprintf("\n");
            return R_NilValue;
        }

      to display the sequence of an XString object (in hexadecimal format).
      Don't forget to register the print_XString_bytes() function
      if you want to make it a .Call entry point!

   f. 2 IMPORTANT THINGS TO REMEMBER ABOUT XString OBJECTS:
        o they are NOT null-terminated like standard strings in C: they can
          contain the null byte so you should never use the C standard string
          functions on them;
        o DNAString and RNAString objects have their data ENCODED: for
          example, if you know that the 'xstring' argument in the above code
          will always point to a DNAString instance, then you could replace
            Rprintf("%x ", *x_char);
          by
            Rprintf("%x(%c) ", *x_char, DNAdecode(*x_char));
          Note that this code will work properly only if 'xstring' is a
          DNAString instance!

   Please consult the "System and foreign language interfaces" section in the
   Writing R Extensions manual for more information:

     http://cran.r-project.org/doc/manuals/R-exts.html

 *****************************************************************************/
#include "Biostrings_defines.h"


/*
 * Low-level manipulation of the extendable buffers.
 */

CharBuf new_CharBuf_from_string(
	const char *string
);

CharBBuf new_CharBBuf(
	int buflength,
	int nelt
);

void append_string_to_CharBBuf(
	CharBBuf *cbbuf,
	const char *string
);


/*
 * Creates an IRanges object from a set of sequences (only the lengths of the
 * sequences are used).
 */
SEXP new_IRanges_from_RoSeqs(const char *class, RoSeqs seqs);


/*
 * Low-level manipulation of XString and XStringSet objects.
 */
char DNAencode(char c);

char DNAdecode(char code);

char RNAencode(char c);

char RNAdecode(char code);

RoSeq get_XString_asRoSeq(SEXP x);

const char *get_XStringSet_baseClass(SEXP x);

int get_XStringSet_length(SEXP x);

CachedXStringSet new_CachedXStringSet(SEXP x);

RoSeq get_CachedXStringSet_elt_asRoSeq(
	CachedXStringSet *x,
	int i
);

RoSeq get_XStringSet_elt_asRoSeq(
	SEXP x,
	int i
);

SEXP new_XStringSet_from_RoSeqs(
	const char *baseClass,
	RoSeqs seqs
);

void set_XStringSet_names(
	SEXP x,
	SEXP names
);

SEXP alloc_XStringSet(
	const char *baseClass,
	int length,
	int super_length
);

void write_RoSeq_to_CachedXStringSet_elt(
	CachedXStringSet *x,
	int i,
	const RoSeq *seq,
	int encode
);

void write_RoSeq_to_XStringSet_elt(
	SEXP x,
	int i,
	const RoSeq *seq,
	int encode
);


/*
 * Converting a set of sequences from one internal representation into another.
 */

RoSeqs new_RoSeqs_from_BBuf(CharBBuf cbbuf);

SEXP new_STRSXP_from_RoSeqs(RoSeqs seqs, SEXP lkup);


/*
 * Match reporting facilities.
 *
 * init_match_reporting():
 *   Look at Biostrings_defines.h in this folder for the valid values of the
 *   'mrmode' arg (the match reporting mode).
 *
 * report_match():
 *   In mode COUNT_MRMODE, it ignores the values of its 'start' and 'end' args.
 *   In mode START_MRMODE, it ignores the value of its 'end' arg.
 *   So yes, for now, the 'end' arg is always ignored, but modes that look at
 *   this argument will be added in the future.
 *
 * reported_matches_asSEXP():
 *   The SEXP returned by reported_matches_asSEXP() is an integer vector
 *   of length 1 in mode COUNT_MRMODE and of length the number of matches
 *   in mode START_MRMODE. IMPORTANT: It is returned UNPROTECTED!
 */
void init_match_reporting(int mrmode);

int report_match(int start, int end);

SEXP reported_matches_asSEXP();

