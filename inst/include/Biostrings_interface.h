/*****************************************************************************
 Biostrings C interface
 ----------------------

 This file contains the prototypes of a subset of Biostrings C routines that
 can be called by C code in other packages.
 In order to use these routines in your package, you need to do the following:

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

        SEXP print_BString_charseq_as_bytes(SEXP x_BString)
        {
            const char *S, *current_Schar;
            int nS, i;

            S = get_BString_charseq(x_BString, &nS);
            for (i = 0, current_Schar = S; i < nS; i++, current_Schar++)
                Rprintf("%x ", *current_Schar);
            Rprintf("\n");
            return R_NilValue;
        }

      to display the sequence of a BString object (in hexadecimal format).
      Don't forget to register the print_BString_charseq_as_bytes() function
      if you want to make it a .Call entry point!

   f. 3 IMPORTANT THINGS TO REMEMBER ABOUT BString OBJECTS:
        o their data are IMMUTABLE: this is enforced by having
          get_BString_charseq() returning a const char *
        o they are NOT null-terminated like standard strings in C: they can
          contain the null byte so you should never use the C standard string
          functions on them
        o DNAString and RNAString objects have their data ENCODED: for
          example, if you know that the x_BString SEXP in the code above will
          always be containing a DNAString instance, then you could replace
            Rprintf("%x ", *current_Schar);
          by
            Rprintf("%x(%c) ", *current_Schar, DNAdecode(*current_Schar));
          Note that this code will only work if x_BString is a DNAString
          instance!

 Please consult the "System and foreign language interfaces" section in the
 Writing R Extensions manual for more information:

   http://cran.r-project.org/doc/manuals/R-exts.html#System-and-foreign-language-interfaces

 *****************************************************************************/
#include "Biostrings_defines.h"

char DNAencode(char c);

char DNAdecode(char code);

char RNAencode(char c);

char RNAdecode(char code);

const char *get_BString_charseq(SEXP x, int *length);

/*
 * Look at Biostrings_defines.h in this folder for the valid values of the
 * 'mrmode' arg (the match reporting mode).
 */
void init_match_reporting(int mrmode);

/*
 * In mode COUNT_MRMODE, the values of the 'start' and 'end' args are ignored.
 * In mode START_MRMODE, the value of the 'end' arg is ignored.
 * So yes, for now, the 'end' arg is always ignored, but modes that look at
 * this argument will be added in the future.
 */
int report_match(int start, int end);

/*
 * The SEXP returned by reported_matches_asSEXP() is unprotected! It's an
 * integer vector of length 1 in mode COUNT_MRMODE and of length the number of
 * matches in mode START_MRMODE.
 */
SEXP reported_matches_asSEXP();

// more functions to come soon...

