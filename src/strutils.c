/*
 * Copyright R. Gentleman, W. Huber 2003-2007, all rights reserved
 *
 */

#include "Biostrings.h"

#include <R.h>
#include <Rinternals.h>
#include "R_ext/Arith.h"
#include "R_ext/Error.h"
#include "R_ext/Applic.h" /* machar */

#include <string.h>
#include <stdlib.h>

/*--------------------------*/
/*  the complementary base  */
/*--------------------------*/
char compbase(char c) {
  char bases[] = "TACGtacgn";
  char compl[] = "ATGCatgcn";
  char* p;

  p = strchr(bases, (int) c);
  if(p==NULL) {
    error("Character %c does not code for a nucleic acid.", c);
  }
  return(compl[p-bases]);
}

/*------------------------------------------------------
 get longest consecutive stretch of consecutive letters 
 ------------------------------------------------------*/
SEXP MP_longestConsecutive(SEXP x, SEXP letter)
{
  int i, j, ncur, nmax;
  const char *pc, *seq;
  char c;
  SEXP rv;

  /* Check and preprocess function arguments */
  if (!isString(x))
    error("'x' must be a string.");

  if (!isString(letter) || length(letter)!=1)
    error("'letter' must be a character variable of length 1.");
  pc = CHAR(STRING_ELT(letter, 0));
  if(strlen(pc)!=1) {
      error("'letter' must contain exactly one character but contains %d.", 
          strlen(pc));
  }
  c = *pc;

  PROTECT(rv = allocVector(INTSXP, length(x)));

  for(i=0; i<length(x); i++) {    
    if(STRING_ELT(x, i) == NA_STRING){
      nmax = NA_INTEGER;
    } else {
      seq = CHAR(STRING_ELT(x, i));
      nmax = ncur = 0;
      for(j=0; j<strlen(seq); j++) {
	if(seq[j]==c) {
	  ncur++;
	  if(ncur>nmax)
	    nmax=ncur;
	} else {
	  ncur=0;
	}
      }
    }
    INTEGER(rv)[i] = nmax;
  }

  UNPROTECT(1);
  return(rv);
}
