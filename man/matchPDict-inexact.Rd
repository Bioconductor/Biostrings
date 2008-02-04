\name{matchPDict-inexact}

\alias{matchPDict-inexact}

% Classes
\alias{class:TPdna_PDict}
\alias{TPdna_PDict-class}

% Generics and methods:
\alias{initialize,TPdna_PDict-method}
\alias{tail,TPdna_PDict-method}
\alias{show,TPdna_PDict-method}


\title{Inexact matching with matchPDict()}

\description{
  The \code{matchPDict} function can be used for inexact matching.

  For the moment, only one kind of inexact matching is provided where
  the user must specify the length of the "trusted prefix" for each
  pattern stored in the dictionary and the maximum number of mismatching
  letters allowed in each pattern suffix (or "tail").
  The inexact matches are searched by using a 2-stage approach:
  (1) use the Aho-Corasick algorithm for fast lookup of all the prefixes
  in the subject, followed by
  (2) comparison of each tail with the short substrings in the subject
  located next to the positions found in (1).
  Performance will mostly depend on the number of positions found in (1)
  so it's important to keep the length of the "trusted prefix" as large
  as possible in order to minimize this number.
  A simple empiric rule is: if W is the length of the "trusted prefix"
  and L the length of the subject, then L / (4**W) should be as small as
  possible, typically < 10 or 20.

  Note that, when used for inexact matching, \code{matchPDict} supports
  dictionaries with patterns of variable length.
  However, like in exact matching mode, it only works with a dictionary
  of DNA patterns.
}

\author{H. Pages}

\seealso{
  \link{matchPDict-exact}
}

\examples{
  ## WARNING! THIS IS STILL WORK IN PROGRESS!!!
  ## A known issue is that the "startIndex", "[[" and "unlist" methods are
  ## BROKEN on the ViewsIndex object returned by matchPDict() called on a
  ## TPdna_PDict object!

  ## Use the 'width' argument to specify the length of the "trusted prefix"
  pdict <- new("TPdna_PDict", c("ACNG", "GT", "CGT", "AC"), width=2)
  pdict
  tail(pdict) 
  endIndex(matchPDict(pdict, DNAString("ACGGACCG"), max.mismatch=0))
  endIndex(matchPDict(pdict, DNAString("ACGGACCG"), max.mismatch=1))
  vindex <- matchPDict(pdict, DNAString("ACGGACCG"), max.mismatch=2)
  endIndex(vindex)

  startIndex(vindex)   # CURRENTLY BROKEN!
  vindex[[1]]          # CURRENTLY BROKEN!
  unlist(vindex)       # CURRENTLY BROKEN!

  ## Ambiguities (IUPAC extended letters) are supported in the tail only:
  endIndex(matchPDict(pdict, DNAString("ACGGACCG"), fixed=FALSE))
}

\keyword{methods}

