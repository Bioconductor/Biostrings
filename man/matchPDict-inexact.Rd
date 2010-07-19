\name{matchPDict-inexact}

\alias{matchPDict-inexact}


\title{Inexact matching with matchPDict()/countPDict()/whichPDict()}

\description{
  The \code{matchPDict}, \code{countPDict} and \code{whichPDict} functions
  efficiently find the occurrences in a text (the subject) of all patterns
  stored in a preprocessed dictionary.

  This man page shows how to use these functions for inexact (or fuzzy)
  matching or when the original dictionary has a variable width.

  See \code{?\link{matchPDict}} for how to use these functions for exact
  matching of a constant width dictionary i.e. a dictionary where all the
  patterns have the same length (same number of nucleotides).
}

\details{
  In this man page, we assume that you know how to preprocess
  a dictionary of DNA patterns that can then be used with
  \code{matchPDict}, \code{countPDict} or \code{\link{whichPDict}}.
  Please see \code{?\link{PDict}} if you don't.

  \code{matchPDict} and family support different kinds of inexact
  matching but with some restrictions. Inexact matching is controlled
  via the definition of a Trusted Band during the preprocessing step
  and/or via the \code{max.mismatch}, \code{min.mismatch} and \code{fixed}
  arguments.
  Defining a Trusted Band is also required when the original dictionary
  is not rectangular (variable width), even for exact matching.
  See \code{?\link{PDict}} for how to define a Trusted Band.

  Here is how \code{matchPDict} and family handle the Trusted Band
  defined on \code{pdict}:
  \itemize{
  \item (1) Find all the exact matches of all the elements in
        the Trusted Band.
  \item (2) For each element in the Trusted Band that has at
        least one exact match, compare the head and the tail of
        this element with the flanking sequences of the matches
        found in (1).
  }

  Note that the number of exact matches found in (1) will decrease
  exponentially with the width of the Trusted Band.
  Here is a simple guideline in order to get reasonably good
  performance: if TBW is the width of the Trusted Band
  (\code{TBW <- tb.width(pdict)}) and L the number of letters in the
  subject (\code{L <- nchar(subject)}), then \code{L / (4^TBW)} should
  be kept as small as possible, typically < 10 or 20.

  In addition, when a Trusted Band has been defined during preprocessing,
  then \code{matchPDict} and family can be called with \code{fixed=FALSE}.
  In this case, IUPAC ambiguity codes in the head or the tail of the
  \link{PDict} object are treated as ambiguities.

  Finally, \code{fixed="pattern"} can be used to indicate that IUPAC
  ambiguity codes in the subject should be treated as ambiguities.
  It only works if the density of codes is not too high.
  It works whether or not a Trusted Band has been defined on \code{pdict}.
}

\author{H. Pages}

\references{
  Aho, Alfred V.; Margaret J. Corasick (June 1975). "Efficient string
  matching: An aid to bibliographic search".
  Communications of the ACM 18 (6): 333-340.
}

\seealso{
  \link{PDict-class},
  \link{MIndex-class},
  \link{matchPDict}
}

\examples{
  ## ---------------------------------------------------------------------
  ## A. USING AN EXPLICIT TRUSTED BAND
  ## ---------------------------------------------------------------------

  library(drosophila2probe)
  dict0 <- DNAStringSet(drosophila2probe)
  dict0  # the original dictionary

  ## Preprocess the original dictionary by defining a Trusted Band that
  ## spans nucleotides 1 to 9 of each pattern.
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
  count0b <- countPDict(pdict9, chr3R, max.mismatch=0)
  identical(count0b, count0)  # TRUE
  
  ## ---------------------------------------------------------------------
  ## C. USING AN EXPLICIT TRUSTED BAND ON A VARIABLE WIDTH DICTIONARY
  ## ---------------------------------------------------------------------

  ## Here is a small variable width dictionary that contains IUPAC
  ## ambiguities (pattern 1 and 3 contain an N):
  dict0 <- DNAStringSet(c("TACCNG", "TAGT", "CGGNT", "AGTAG", "TAGT"))
  ## (Note that pattern 2 and 5 are identical.)

  ## If we only want to do exact matching, then it is recommended to use
  ## the widest possible Trusted Band i.e. to set its width to
  ## 'min(width(dict0))' because this is what will give the best
  ## performance. However, when 'dict0' contains IUPAC ambiguities (like
  ## in our case), it could be that one of them is falling into the
  ## Trusted Band so we get an error (only base letters can go in the
  ## Trusted Band for now):
  \dontrun{
    PDict(dict0, tb.end=min(width(dict0)))  # Error!
  }

  ## In our case, the Trusted Band cannot be wider than 3:
  pdict <- PDict(dict0, tb.end=3)
  tail(pdict)

  subject <- DNAString("TAGTACCAGTTTCGGG")

  m <- matchPDict(pdict, subject)
  countIndex(m)  # pattern 2 and 5 have 1 exact match
  m[[2]]

  ## We can take advantage of the fact that our Trusted Band doesn't cover
  ## the entire dictionary to allow inexact matching on the uncovered parts
  ## (the tail in our case):

  m <- matchPDict(pdict, subject, fixed=FALSE)
  countIndex(m)  # now pattern 1 has 1 match too
  m[[1]]

  m <- matchPDict(pdict, subject, max.mismatch=1)
  countIndex(m)  # now pattern 4 has 1 match too
  m[[4]]

  m <- matchPDict(pdict, subject, max.mismatch=1, fixed=FALSE)
  countIndex(m)  # now pattern 3 has 1 match too
  m[[3]]  # note that this match is "out of limit"
  Views(subject, m[[3]])

  m <- matchPDict(pdict, subject, max.mismatch=2)
  countIndex(m)  # pattern 4 gets 1 additional match
  m[[4]]

  ## Unlist all matches:
  unlist(m)

  ## ---------------------------------------------------------------------
  ## D. WITH IUPAC AMBIGUITY CODES IN THE SUBJECT
  ## ---------------------------------------------------------------------
  pdict <- PDict(c("ACAC", "TCCG"))
  as.list(matchPDict(pdict, DNAString("ACNCCGT")))
  as.list(matchPDict(pdict, DNAString("ACNCCGT"), fixed="pattern"))
  as.list(matchPDict(pdict, DNAString("ACWCCGT"), fixed="pattern"))
  as.list(matchPDict(pdict, DNAString("ACRCCGT"), fixed="pattern"))
  as.list(matchPDict(pdict, DNAString("ACKCCGT"), fixed="pattern"))

  dict <- DNAStringSet(c("TTC", "CTT"))
  pdict <- PDict(dict)
  subject <- DNAString("CYTCACTTC")
  mi1 <- matchPDict(pdict, subject, fixed="pattern")
  mi2 <- matchPDict(dict, subject, fixed="pattern")
  stopifnot(identical(as.list(mi1), as.list(mi2)))
}

\keyword{methods}
