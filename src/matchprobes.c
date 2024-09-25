/*
 * Copyright R. Gentleman and W. Huber 2002-2003, all rights reserved
 * 
 * A program to match probes on a DNA chip
 * to a given sequence
 */

#include "Biostrings.h"

#include <R.h>
#include <Rinternals.h>
#include "R_ext/Arith.h"
#include "R_ext/Error.h"
#include "R_ext/Applic.h" /* machar */
#include "R_ext/RS.h" 		/* CallocCharBuf */

#include <string.h>
#include <stdlib.h>

/* The position (counting starts from 0) within a short 
   oligonucleotide probe on an Affymetrix genechip where
   the base is flipped between "match" and "mismatch" probes */
#define MISMATCHPOSITION 12

/*------------------------------------------------------------*/
/* A data structure for storing a match result between a pair */
/* of sequences                                               */
/* pos1: the position of the match in sequence 1              */
/* pos2: the position of the match in sequence 2              */
/* len:  the length of the match                              */
/* type: type of match (e.g. PM or MM )                       */
/*------------------------------------------------------------*/
typedef struct {
  int pos1;
  int pos2;
  int len;
  int type;
}  Match; 

typedef struct {
  int rec;
  Match match;
} MatchWithRec;


/* The general signature of sequence match functions       */
/* and one particular instance:                            */
/* typedef Match* (*R_MatchFun)(char*, char*, struct Matchres*); */

/*--------------------------------------------------------------------*/
/* Test whether a sequence x contains a substring y or a substring    */
/* that is obtained from y by flipping its 12th position              */
/*                                                                    */
/* Arguments:                                                         */
/* INPUT                                                              */
/*  x, y: pointers to the "long" (haystack) and the short (needle)    */
/*        sequence, respectively                                      */
/* OUTPUT                                                             */
/*  rv:   pointer to a struct Match (see above)                       */
/*--------------------------------------------------------------------*/
void strstr_with_pmormm(const char* x, const char* y, Match* rv) {
  char *m;
  char *scratch;
  int len;

  /* this the return value if nothing is found */  
  rv->pos1 = 0;
  rv->pos2 = 0;
  rv->len  = 0;
  rv->type = 0;

  if ((m = strstr(x, y))) {
    rv->pos1 = (m - x)/sizeof(char) + 1;
    rv->pos2 = 1;
    rv->len  = strlen(y);
    rv->type = 1;
    return;
  }

  len = strlen(y);
  if (len < MISMATCHPOSITION) {
      error("Sequence y is too short: must at least have length %d.", 
            MISMATCHPOSITION);
  }

  scratch = CallocCharBuf(len);
  strcpy(scratch, y);
  scratch[MISMATCHPOSITION] = compbase(scratch[MISMATCHPOSITION]);

  if ((m = strstr(x, scratch))) {
    rv->pos1 = (m - x)/sizeof(char) + 1;
    rv->pos2 = 1;
    rv->len  = len;
    rv->type = 2;
  }

  R_Free(scratch);
  return;
}

/*---------------------------------------------------------------*/
/* R interface:                                                  */
/* query:    a vector of mode 'character'                        */
/* probepos: a logical. If TRUE, return a list with the position */
/*           of the match too                                    */
/*---------------------------------------------------------------*/
SEXP MP_matchprobes(SEXP query, SEXP records, 
		    SEXP probepos)
     /*		      SEXP probepos, SEXP method)  */
{
  int i, k, ct, rvlen;

  /* R vectors that are used in the construction of the return value */
  SEXP matchesrec, matchespos, rv, rvnames, rvmatch, rvpos;

  /* query and record sequences */
  const char *queryseq, *recseq;
  /* their length */
  int nrqu, nrrec;

  /* logical: should the match positions also be returned?  */
  int lprobepos;

  /* intermediate buffer to store the match results before  */
  /* bringing them into matchespos and matchesres  */
  MatchWithRec *mbuf;
  Match m;

  /*
  if(Rf_isFunction(method)) {
    matchfun = (R_DistFun) myRDist;
    isFun = 1;
  } else
    matchfun = (R_DistFun) R_ExternalPtrAddr(method);
  */

  if(!isString(query) )
    error("Argument query must be a string");

  if(!isLogical(probepos))
    error("Argument probepos must be logical.");
  lprobepos = asLogical(probepos);
  
  nrqu  = length(query);
  nrrec = length(records);

  if(lprobepos) {
    rvlen = 2;
    PROTECT(rvpos = allocVector(VECSXP, nrqu));
  } else {
    rvlen = 1;
    rvpos = NULL;
  }
  PROTECT(rv      = allocVector(VECSXP, rvlen));
  PROTECT(rvnames = allocVector(VECSXP, rvlen));
  PROTECT(rvmatch = allocVector(VECSXP, nrqu));

  /* allocate memory for the match buffer */
  mbuf = R_Calloc(nrrec, MatchWithRec);

  /* loop over the query sequences in 'query' */
  for(k=0; k<nrqu; k++) {  

    R_CheckUserInterrupt();

    ct = 0;
    if(STRING_ELT(query, k) != NA_STRING) {
      queryseq = CHAR(STRING_ELT(query, k));

      /* loop over the sequences on the array */    
      for (i=0; i<nrrec; i++) {
	if(STRING_ELT(records, i) != NA_STRING) {
	  recseq = CHAR(STRING_ELT(records, i));
	
	  strstr_with_pmormm(queryseq, recseq, &m);
	  if (m.type!=0) {
	    mbuf[ct].rec = i+1;
	    if (m.type==2) 
	      mbuf[ct].rec *= -1;
	    mbuf[ct].match.pos1 =  m.pos1;
	    ct++;
	  } /* if */
	} /* if */
      } /* for i */
    } /* if */
    matchesrec = allocVector(INTSXP, ct); 
    SET_VECTOR_ELT(rvmatch, k, matchesrec);
    for(i=0; i<ct; i++)
      INTEGER(matchesrec)[i] = mbuf[i].rec;

    if(lprobepos){
      matchespos = allocVector(INTSXP, ct);
      SET_VECTOR_ELT(rvpos, k, matchespos);
      for(i=0; i<ct; i++)
	INTEGER(matchespos)[i] = mbuf[i].match.pos1;
    }
  } /* for k */

  SET_VECTOR_ELT(rv,      0, rvmatch);
  SET_VECTOR_ELT(rvnames, 0, mkChar("match"));
  if(lprobepos) {
    SET_VECTOR_ELT(rv,      1, rvpos);
    SET_VECTOR_ELT(rvnames, 1, mkChar("pos"));
  }
  setAttrib(rv, R_NamesSymbol, rvnames);

  R_Free(mbuf);
  UNPROTECT(2+rvlen);
  return(rv);
}
