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

static char dna_base_comp(char c) {
  char bases[] = "TACGNtacgn";
  char compl[] = "ATGCNatgcn";
  char* p;

  p = strchr(bases, (int) c);
  if (p == NULL) {
    error("'%c' does not code for DNA.", c);
  }
  return compl[p - bases];
}

static char rna_base_comp(char c) {
  char bases[] = "UACGNuacgn";
  char compl[] = "AUGCNaugcn";
  char* p;

  p = strchr(bases, (int) c);
  if (p == NULL) {
    error("'%c' does not code for RNA.", c);
  }
  return compl[p - bases];
}

static void rev_comp(const char *in, int len, char *out, char (*base_comp)(char))
{
    int i, j;

    for(i = 0, j = len - 1; j >= 0; i++, j--) {
	out[i] = (*base_comp)(in[j]);
    }
}

static SEXP do_revcomp(SEXP x, char (*base_comp)(char))
{
    SEXP r, xk;
    char *rev;
    const char *orig;
    int k, n, m;

    if( !isString(x) )
        error("argument must be a string");

    n = length(x);
    PROTECT(r = allocVector(STRSXP, n));
    for (k = 0; k < n; k++ ) {
        xk = STRING_ELT(x, k);
        if (xk == NA_STRING) {
            SET_STRING_ELT(r, k, NA_STRING);
        } else {
            m = length(xk);
            /* CallocCharBuf itself takes care of the +1 for the \0 */
            rev = CallocCharBuf(m);
            orig = CHAR(xk);
            rev_comp(orig, m, rev, base_comp);
            SET_STRING_ELT(r, k, mkChar(rev));
            Free(rev);
        }
    }
    UNPROTECT(1);
    return r;
}

/* R interface for reverse RNA complement */
SEXP MP_rna_revcomp(SEXP x)
{
    return do_revcomp(x, &rna_base_comp);
}

/* R interface for reverse DNA complement */
SEXP MP_dna_revcomp(SEXP x)
{
    return do_revcomp(x, &dna_base_comp);
}

/*------------------------------------------*/
/* R interface                              */
/* reverse all elements of the input string */
/*------------------------------------------*/
SEXP MP_revstring(SEXP x)
{
  SEXP r, xk;
  char *rev;
  const char *orig;
  int i, j, k, n, m;

  if( !isString(x) )
    error("argument must be a string");

  n = length(x);   
  PROTECT(r = allocVector(STRSXP, n));
  for(k=0; k<n; k++ ) {
    xk = STRING_ELT(x, k);
    if(xk == NA_STRING){
      SET_STRING_ELT(r, k, NA_STRING);
    } else {
      m = length(xk);
      /* CallocCharBuf itself takes care of the +1 for the \0 */
      rev = CallocCharBuf(m);  
      orig = CHAR(xk);
      for(i=0, j = m-1; j>=0; i++, j--) {
	rev[i] = orig[j];
      }
      SET_STRING_ELT(r, k, mkChar(rev));
      Free(rev);
    }
  }
  UNPROTECT(1);
  return(r);
} 

/*------------------------------------------------*/
/* replace an element with its complementary base */
/*------------------------------------------------*/
SEXP MP_complementSeq(SEXP x, SEXP start, SEXP stop)
{
  SEXP r, xk;
  const char *orig;
  char *rev;
  int i, i0, i1, j0, j1, k, n, m;

  /* Bureaucracy */
  if (!isString(x))
    error("'x' must be a string.");
  if (!isInteger(start) || length(start)!=1)
    error("'start' must be an integer variable of length 1.");
  if (!isInteger(stop) || length(stop)!=1)
    error("'stop' must be an integer variable of length 1.");

  i0 = *INTEGER(start)-1;
  i1 = *INTEGER(stop);
  if (i0 < 0)
    error("'start' must be >=1.");
  if (i1 < 0)
    error("'stop' must be >=0.");

  n = length(x);  
  PROTECT(r = allocVector(STRSXP, n));
  for(k=0; k<n; k++) {   
    xk = STRING_ELT(x, k);
    if(xk == NA_STRING){
      SET_STRING_ELT(r, k, NA_STRING);
    } else {
      m = length(xk);
      rev = CallocCharBuf(m);
      orig = CHAR(xk);
      j0 = (i0<m) ? i0 : m;
      j1 = (i1==0) ? m : ((i1<m) ? i1 : m);
      for(i=0; i<j0; i++) 
	rev[i] = orig[i];
      for(i=j0; i<j1; i++) 
	rev[i] = compbase(orig[i]);
      for(i=j1; i<m; i++) 
	rev[i] = orig[i];
      SET_STRING_ELT(r, k, mkChar(rev));
      Free(rev);
    }
  }
  UNPROTECT(1);
  return(r);
}

/*------------------------------------------------*/
/* get the CGAT content of a sequence             */
/*------------------------------------------------*/
SEXP MP_basecontent(SEXP x, SEXP dna)
{
  SEXP rv, rownames, colnames, dimnames, dim;
  const char *seq;
  int i, j, n, ia, ic, ig, it;
  int using_dna;

  if( !isString(x) )
    error("argument must be a string");

  if (!isLogical(dna))
      error("argument 'dna' must be TRUE/FALSE");

  using_dna = LOGICAL(dna)[0];
  if (using_dna == NA_LOGICAL)
      using_dna = 1;

  n = length(x);
  PROTECT(rv = allocVector(INTSXP, n*4));

  for(i=0; i<n; i++) {
    if(STRING_ELT(x, i) == NA_STRING){
      ia = ic = ig = it = NA_INTEGER;
    } else {
      ia = ic = ig = it = 0;
      seq = CHAR(STRING_ELT(x, i));
      for(j=0; j<strlen(seq); j++) {
	switch(seq[j]) {
	case 'a': 
	case 'A':
	  ia++;
	  break;
	case 'c':
	case 'C':
	  ic++;
	  break;
	case 'g':
	case 'G':
	  ig++;
	  break;
	case 't':
	case 'T':
          if (!using_dna)
              error("unknown base '%c' in row %d, col %d", seq[j], i+1, j+1);
	  it++;
	  break;
	case 'u':
	case 'U':
          if (using_dna)
              error("unknown base '%c' in row %d, col %d", seq[j], i+1, j+1);
	  it++;
	  break;
	default:
	  error("Unknown base %c in row %d, column %d.", seq[j], i+1, j+1);
	}
      }
    }
    INTEGER(rv)[i    ] = ia; 
    INTEGER(rv)[i+n  ] = it; 
    INTEGER(rv)[i+n*2] = ic; 
    INTEGER(rv)[i+n*3] = ig; 
  }

  /* dim */
  PROTECT(dim = allocVector(INTSXP, 2));
  INTEGER(dim)[0] = n;
  INTEGER(dim)[1] = 4;
  setAttrib(rv, R_DimSymbol, dim);

  /* dim names */
  PROTECT(colnames = allocVector(STRSXP, 4));
  SET_STRING_ELT(colnames, 0, mkChar("A"));
  if (using_dna)
      SET_STRING_ELT(colnames, 1, mkChar("T"));
  else
      SET_STRING_ELT(colnames, 1, mkChar("U"));
  SET_STRING_ELT(colnames, 2, mkChar("C"));
  SET_STRING_ELT(colnames, 3, mkChar("G"));

  /* dim names */
  PROTECT(rownames = allocVector(STRSXP, n));
  PROTECT(dimnames = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(dimnames, 0, rownames);
  SET_VECTOR_ELT(dimnames, 1, colnames);
  setAttrib(rv, R_DimNamesSymbol, dimnames);

  UNPROTECT(5);
  return(rv);
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
