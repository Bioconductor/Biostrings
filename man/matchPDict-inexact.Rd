\name{matchPDict-inexact}

\alias{matchPDict-inexact}


\title{Inexact matching with matchPDict()/countPDict()/whichPDict()}

\description{
  The \code{matchPDict}, \code{countPDict} and \code{whichPDict} functions
  efficiently find the occurrences in a text (the subject) of all patterns
  stored in a preprocessed dictionary.

  This man page shows how to use these functions for inexact matching
  or when the original dictionary has a variable width.

  See \code{?\link{matchPDict}} for how to use these functions for exact
  matching of a constant width dictionary i.e. a dictionary where all the
  patterns have the same length (same number of nucleotides).
}

\details{
  In this man page, we assume that you know how to preprocess
  a dictionary of DNA patterns that can then be used with
  \code{matchPDict}, \code{countPDict} or \code{\code{whichPDict}}.
  Please see \code{?\link{PDict}} if you don't.

  When using \code{matchPDict}, \code{countPDict} or \code{whichPDict}
  for inexact matching or when the original dictionary has a variable width,
  a Trusted Band must be defined during the preprocessing step.
  This is done thru the \code{tb.start}, \code{tb.end} and \code{tb.width}
  arguments of the \code{\link{PDict}} constructor (see \code{?\link{PDict}}
  for the details).

  Then \code{matchPDict}/\code{countPDict}/\code{whichPDict} can be called
  with a null or non-null \code{max.mismatch} value and the search for exact
  or inexact matches happens in 2 steps:
  (1) find all the exact matches of all the elements in the Trusted Band;
      then
  (2) for each element in the Trusted Band that has at least one exact
      match, compare the head and the tail of this element
      with the flanking sequences of the matches found in (1).

  Note that the number of exact matches found in (1) will decrease
  exponentially with the width of the Trusted Band.
  Here is a simple guideline in order to get reasonably good
  performance: if TBW is the width of the Trusted Band
  (\code{TBW <- tb.width(pdict)}) and L the number of letters in the
  subject (\code{L <- nchar(subject)}), then \code{L / (4^TBW)} should
  be kept as small as possible, typically < 10 or 20.

  In addition, when a Trusted Band has been defined during preprocessing,
  then \code{matchPDict}/\code{countPDict}/\code{whichPDict} can be called
  with \code{fixed=FALSE}. In this case, IUPAC extended letters in the head
  or the tail of the \link{PDict} object are treated as ambiguities.
}

\author{H. Pages}

\seealso{
  \link{PDict-class},
  \link{MIndex-class},
  \link{matchPDict}
}

\examples{
  ## WARNING! THIS IS STILL WORK IN PROGRESS!!!
  ## A known issue is that the "startIndex", "[[" and "unlist" methods are
  ## BROKEN on the MIndex object returned by matchPDict() called on a
  ## PDict object with a head or a tail!

  ## ---------------------------------------------------------------------
  ## A. WITH UNNAMED PATTERNS
  ## ---------------------------------------------------------------------

  library(drosophila2probe)
  dict0 <- DNAStringSet(drosophila2probe$sequence)
  dict0  # the original dictionary

  ## Preprocess the original dictionary and define a Trusted Band that
  ## goes from the 1st to the 9h nucleotide of each of its elements.
  pdict9 <- PDict(dict0, tb.end=9)
  pdict9
  tail(pdict9)
  sum(duplicated(pdict9))
  table(patternFrequency(pdict9))

  library(BSgenome.Dmelanogaster.UCSC.dm3)
  chr3R <- Dmelanogaster$chr3R
  chr3R
  table(countPDict(pdict9, chr3R, max.mismatch=1))
  table(countPDict(pdict9, chr3R, max.mismatch=3))
  table(countPDict(pdict9, chr3R, max.mismatch=5))

  ## ---------------------------------------------------------------------
  ## B. COMPARISON WITH EXACT MATCHING
  ## ---------------------------------------------------------------------

  ## When the original dictionary is of constant width, exact matching
  ## (i.e. 'max.mismatch=0' and 'fixed=TRUE) will be more efficient with
  ## a full-width Trusted Band (i.e. a Trusted Band that covers the entire
  ## dictionary) than with a Trusted Band of width < width(dict0).
  pdict0 <- PDict(dict0)
  count0 <- countPDict(pdict0, chr3R)
  count_tp9_mm0 <- countPDict(pdict9, chr3R, max.mismatch=0)
  identical(count_tp9_mm0, count0)    # TRUE
  
  ## ---------------------------------------------------------------------
  ## C. VARIABLE WIDTH DICTIONARY
  ## ---------------------------------------------------------------------

  pdict <- PDict(c("ACNG", "GT", "CGT", "AC"), tb.end=2)
  pdict
  tail(pdict) 
  endIndex(matchPDict(pdict, DNAString("ACGGACCG"), max.mismatch=0))
  endIndex(matchPDict(pdict, DNAString("ACGGACCG"), max.mismatch=1))
  mindex <- matchPDict(pdict, DNAString("ACGGACCG"), max.mismatch=2)
  endIndex(mindex)

  startIndex(mindex)                  # CURRENTLY BROKEN!
  mindex[[1]]                         # CURRENTLY BROKEN!
  unlist(mindex)                      # CURRENTLY BROKEN!

  ## ---------------------------------------------------------------------
  ## D. WITH AMBIGUITIES IN THE TAIL
  ## ---------------------------------------------------------------------

  ## Ambiguities (IUPAC extended letters) are supported in the tail only:
  endIndex(matchPDict(pdict, DNAString("ACGGACCG"), fixed=FALSE))
}

\keyword{methods}
