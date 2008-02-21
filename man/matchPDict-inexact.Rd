\name{matchPDict-inexact}

\alias{matchPDict-inexact}

% Classes
\alias{class:TBdna_PDict}
\alias{TBdna_PDict-class}

% Generics and methods:
\alias{initialize,TBdna_PDict-method}
\alias{head,TBdna_PDict-method}
\alias{tail,TBdna_PDict-method}
\alias{show,TBdna_PDict-method}
\alias{duplicated,TBdna_PDict-method}
\alias{patternFrequency,TBdna_PDict-method}


\title{Inexact matching with matchPDict()}

\description{
  The \code{matchPDict} function can be used for inexact matching.

  For the moment, only one kind of inexact matching is provided where
  the user must specify the length of the Trusted Prefix for each
  pattern stored in the dictionary and the maximum number of mismatching
  letters allowed in each pattern suffix (or "tail"). Please refer
  to the \code{\link{PDict}} documentation for how to create a
  Trusted Prefix DNA dictionary (preprocessing).

  The inexact matches are searched by using a 2-stage approach:
  (1) use the Aho-Corasick algorithm for fast lookup of all the prefixes
  in the subject, followed by
  (2) comparison of each tail with the short substrings in the subject
  located next to the positions found in (1).
  Because performance will mostly depend on the number of positions found in
  (1) it's important to minimize this number by keeping the length of the
  Trusted Prefix as large as possible.
  A simple empiric rule is: if TP is the length of the Trusted Prefix
  and L the length of the subject, then L / (4**TP) should be as small as
  possible, typically < 10 or 20.

  Note that, when \code{matchPDict} is used for inexact matching,
  dictionaries with patterns of variable length are accepted.
  In addition, when \code{fixed=FALSE}, then IUPAC extended letters placed
  in the tail of the dictionary are interpreted as ambiguities.

  However, like in exact matching mode, it only works with a dictionary
  of DNA patterns.
}

\author{H. Pages}

\seealso{
  \link{PDict},
  \link{matchPDict-exact}
}

\examples{
  ## WARNING! THIS IS STILL WORK IN PROGRESS!!!
  ## A known issue is that the "startIndex", "[[" and "unlist" methods are
  ## BROKEN on the ViewsIndex object returned by matchPDict() called on a
  ## TBdna_PDict object!

  ## -------------------------------------------------------------------------
  ## A. WITH UNNAMED PATTERNS
  ## -------------------------------------------------------------------------

  library(drosophila2probe)
  dict0 <- drosophila2probe$sequence  # The input dictionary.

  ## Specify the length of the Trusted Prefix via the 'tb.width' argument
  pdict9 <- PDict(dict0, tb.width=9)
  pdict9
  tail(pdict9)
  #sum(duplicated(pdict9)) # not ready yet!
  #table(patternFrequency(pdict9)) # not ready yet!

  library(BSgenome.Dmelanogaster.FlyBase.r51)
  chr3R <- Dmelanogaster[["3R"]]
  chr3R
  table(countPDict(pdict9, chr3R, max.mismatch=1))
  table(countPDict(pdict9, chr3R, max.mismatch=3))
  table(countPDict(pdict9, chr3R, max.mismatch=5))

  ## -------------------------------------------------------------------------
  ## B. COMPARISON WITH EXACT MATCHING
  ## -------------------------------------------------------------------------

  ## When the input dictionary is of constant width, using 'max.mismatch=0'
  ## on the Trusted Prefix dictionary is equivalent to (but much slower than)
  ## using the full dictionary.
  pdict0 <- PDict(dict0)
  count0 <- countPDict(pdict0, chr3R)
  count_tp9_mm0 <- countPDict(pdict9, chr3R, max.mismatch=0)
  identical(count_tp9_mm0, count0)    # TRUE
  
  ## -------------------------------------------------------------------------
  ## C. VARIABLE LENGTH DICTIONARY
  ## -------------------------------------------------------------------------

  pdict <- PDict(c("ACNG", "GT", "CGT", "AC"), tb.width=2)
  pdict
  tail(pdict) 
  endIndex(matchPDict(pdict, DNAString("ACGGACCG"), max.mismatch=0))
  endIndex(matchPDict(pdict, DNAString("ACGGACCG"), max.mismatch=1))
  vindex <- matchPDict(pdict, DNAString("ACGGACCG"), max.mismatch=2)
  endIndex(vindex)

  startIndex(vindex)                  # CURRENTLY BROKEN!
  vindex[[1]]                         # CURRENTLY BROKEN!
  unlist(vindex)                      # CURRENTLY BROKEN!

  ## -------------------------------------------------------------------------
  ## D. WITH AMBIGUITIES IN THE TAIL
  ## -------------------------------------------------------------------------

  ## Ambiguities (IUPAC extended letters) are supported in the tail only:
  endIndex(matchPDict(pdict, DNAString("ACGGACCG"), fixed=FALSE))
}

\keyword{methods}

