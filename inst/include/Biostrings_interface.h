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
            cachedCharSeq x;
            int i;
            const char *x_char;

            x = cache_XRaw(xstring);
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
 * Low-level manipulation of XString and XStringSet objects.
 */

char DNAencode(char c);

char DNAdecode(char code);

char RNAencode(char c);

char RNAdecode(char code);

int get_XStringSet_length(SEXP x);

const char *get_XStringSet_xsbaseclassname(SEXP x);

cachedXStringSet cache_XStringSet(SEXP x);

int get_cachedXStringSet_length(const cachedXStringSet *cached_x);

cachedCharSeq get_cachedXStringSet_elt(
	const cachedXStringSet *cached_x,
	int i
);

void set_XStringSet_names(
	SEXP x,
	SEXP names
);


/*
 * Match reporting facilities.
 */

void init_match_reporting(const char *ms_mode, int nPSpair);

void set_active_PSpair(int PSpair_id);

void set_match_shift(int shift);

void report_match(int start, int width);

void drop_reported_matches();

int get_match_count();

SEXP reported_matches_asSEXP();


/*
 * MIndex abstract accessor functions.
 */

cachedMIndex cache_MIndex(SEXP x);

int get_cachedMIndex_length(const cachedMIndex *cached_x);

int get_cachedMIndex_elt_width0(const cachedMIndex *cached_x, int i);

cachedIRanges get_cachedMIndex_elt(const cachedMIndex *cached_x, int i);


/*
 * A BOYER-MOORE-LIKE MATCHING ALGO
 */

int match_pattern_boyermoore(
	const cachedCharSeq *P,
	const cachedCharSeq *S,
	int nfirstmatches,
	int walk_backward
);

