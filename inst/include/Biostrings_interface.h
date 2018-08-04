/*****************************************************************************
 Biostrings C interface: prototypes
 ----------------------------------

   The Biostrings C interface is split in 2 files:
     1. Biostrings_defines.h (in this directory): contains the typedefs and
        defines of the interface.
     2. Biostrings_interface.h (this file): contains the prototypes of the
        Biostrings C routines that are part of the interface.

 -----------------------------------------------------------------------------

   This file contains the prototypes of a subset of Biostrings C routines that
   can be called by C code from other packages. In order to use these routines
   in your package, you need to do the following:

   a. Add the IRanges, XVector, and Biostrings package to the 'LinkingTo:'
      field of your DESCRIPTION file. Each of them should also be listed
      either in the 'Depends:' or in the 'Imports:' field, or in both (despite
      'R CMD check' NOTE discouraging this). Finally it's recommended that you
      have the methods package in 'Imports:', preferably in first position.

   b. If it doesn't already have one, add a NAMESPACE file to your package and
      start it with the following lines:

        useDynLib(mypackage) # replace mypackage by the name of your package
        import(methods)      # if you have it in your 'Imports:' field
        import(IRanges)      # if you need to import stuff from IRanges
        import(XVector)      # if you need to import stuff from XVector
        import(Biostrings)   # if you need to import stuff from Biostrings

   c. Add a Biostrings_stubs.c file to your src/ folder with the following
      line (and nothing else) in it:

        #include "_Biostrings_stubs.c" // note the underscore!

      If you also need to call C routines defined in XVector (e.g. hold_XRaw),
      IRanges (e.g. hold_IRanges), and/or S4Vectors (e.g. copy_vector_block),
      then add an XVector_stubs.c, IRanges_stubs.c, and/or S4Vector_stubs.c
      file in a similar fashion.

   d. Add the following line:

        #include "Biostrings_interface.h"

      in your .c files where you need to call Biostrings C routines.
      Add similar lines for XVector_interface.h and/or IRanges_interface.h
      in case you also need to call XVector and/or IRanges C routines.

      Note that those includes should come after any R or system include so
      the top of your .c file will typically look like this:

        // R includes
        #include <Rdefines.h>
        #include <R_ext/Rdynload.h>
        ... maybe more R includes ...

        // System includes
        #include <stdio.h>
        #include <stdlib.h>
        ... maybe more system includes ...

        #include "XVector_interface.h"
        #include "Biostrings_interface.h"

   e. Any C routine defined in Biostrings_interface.h, XVector_interface.h, or
      IRanges_interface.h can now be called in your .c file. For example, you
      could write the following function:

        SEXP print_XString_bytes(SEXP xstring)
        {
            Chars_holder x;
            int i;
            const char *x_char;

            x = hold_XRaw(xstring);
            for (i = 0, x_char = x.ptr; i < x.length; i++, x_char++)
                Rprintf("%x ", *x_char);
            Rprintf("\n");
            return R_NilValue;
        }

      to display the sequence of an XString object (in hexadecimal format).
      Don't forget to register the print_XString_bytes() function if you want
      to make it a .Call entry point!

   f. 2 IMPORTANT THINGS TO REMEMBER ABOUT XString OBJECTS:

      o They are NOT null-terminated like standard strings in C: they can
        contain the null byte so you should never use the C standard string
        functions on them;

      o DNAString and RNAString objects have their data ENCODED: for
        example, if you can assume that the 'xstring' argument in the above
        code will always point to a DNAString instance, then you could replace
            Rprintf("%x ", *x_char);
        by
            Rprintf("%x(%c) ", *x_char, DNAdecode(*x_char));
        Of course this code would work properly only if 'xstring' is actually
        a DNAString instance!

   Please consult the "System and foreign language interfaces" section in the
   Writing R Extensions manual for more information:

     http://cran.r-project.org/doc/manuals/R-exts.html

   Don't hesitate to ask on the bioc-devel mailing list for questions or
   suggestions about this:

     http://bioconductor.org/help/mailing-list/

   Thanks for using the Biostrings package!

 *****************************************************************************/
#include "Biostrings_defines.h"


/*
 * Low-level manipulation of XString objects.
 * (see XString_class.c)
 */

char DNAencode(char c);

char DNAdecode(char code);

char RNAencode(char c);

char RNAdecode(char code);

/*
 * Low-level manipulation of XStringSet objects.
 * (see XStringSet_class.c)
 */

int get_XStringSet_length(SEXP x);

XStringSet_holder hold_XStringSet(SEXP x);

int get_length_from_XStringSet_holder(const XStringSet_holder *x_holder);

Chars_holder get_elt_from_XStringSet_holder(
	const XStringSet_holder *x_holder,
	int i
);

XStringSet_holder get_linear_subset_from_XStringSet_holder(
	const XStringSet_holder *x_holder,
	int offset,
	int length
);

void set_XStringSet_names(
	SEXP x,
	SEXP names
);

/*
 * Low-level manipulation of XStringSetList objects.
 * (see XStringSetList_class.c)
 */

XStringSetList_holder hold_XStringSetList(SEXP x);

int get_length_from_XStringSetList_holder(
	const XStringSetList_holder *x_holder
);

XStringSet_holder get_elt_from_XStringSetList_holder(
	const XStringSetList_holder *x_holder,
	int i
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

MIndex_holder hold_MIndex(SEXP x);

int get_length_from_MIndex_holder(const MIndex_holder *x_holder);

int get_width0_elt_from_MIndex_holder(const MIndex_holder *x_holder, int i);

IRanges_holder get_elt_from_MIndex_holder(const MIndex_holder *x_holder, int i);


/*
 * A BOYER-MOORE-LIKE MATCHING ALGO
 */

int match_pattern_boyermoore(
	const Chars_holder *P,
	const Chars_holder *S,
	int nfirstmatches,
	int walk_backward
);

