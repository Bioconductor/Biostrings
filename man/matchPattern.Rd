\name{matchPattern}

\alias{matchPattern}
\alias{matchPattern,character-method}
\alias{matchPattern,XString-method}
\alias{matchPattern,XStringViews-method}
\alias{matchPattern,MaskedXString-method}
\alias{countPattern}
\alias{countPattern,character-method}
\alias{countPattern,XString-method}
\alias{countPattern,XStringViews-method}
\alias{countPattern,MaskedXString-method}

\alias{vmatchPattern}
\alias{vmatchPattern,character-method}
\alias{vmatchPattern,XStringSet-method}
\alias{vcountPattern}
\alias{vcountPattern,character-method}
\alias{vcountPattern,XStringSet-method}

\alias{matchDNAPattern}


\title{String searching functions}

\description{
  Generic that finds all matches of a pattern in a sequence (an XString object).
}

\usage{
  matchPattern(pattern, subject, algorithm="auto", max.mismatch=0, with.indels=FALSE, fixed=TRUE)
  countPattern(pattern, subject, algorithm="auto", max.mismatch=0, with.indels=FALSE, fixed=TRUE)
  vmatchPattern(pattern, subject, algorithm="auto", max.mismatch=0, with.indels=FALSE, fixed=TRUE)
  vcountPattern(pattern, subject, algorithm="auto", max.mismatch=0, with.indels=FALSE, fixed=TRUE)
}

\arguments{
  \item{pattern}{
    The pattern string.
  }
  \item{subject}{
    For \code{matchPattern} and \code{countPattern}:
    an \link{XString}, \link{XStringViews} or \link{MaskedXString}
    object.

    For \code{vmatchPattern} and \code{vcountPattern}:
    an \link{XStringSet} or \link{XStringViews}
    object.
  }
  \item{algorithm}{
    One of the following: \code{"auto"}, \code{"naive-exact"},
    \code{"naive-inexact"}, \code{"boyer-moore"}, \code{"shift-or"}
    or \code{"indels"}.
  }
  \item{max.mismatch}{
    The maximum number of mismatching letters allowed (see
    \code{\link{isMatching}} for the details).
    If non-zero, an inexact matching algorithm is used.
  }
  \item{with.indels}{
    [DOCUMENT ME]
  }
  \item{fixed}{
    If \code{FALSE} then IUPAC extended letters are interpreted as ambiguities
    (see \code{\link{isMatching}} for the details).
  }
}

\details{
  Available algorithms are: ``naive exact'', ``naive inexact'',
  ``Boyer-Moore-like'' and ``shift-or''. Not all of them can be
  used in all situations: restrictions depend on the length of
  the pattern, the class of the subject and the values of
  \code{max.mismatch} and \code{fixed}.

  When 2 different algorithms can be used for a given task,
  then choosing one or the other only affects the performance,
  not the result, so there is no "wrong choice" (strictly speaking).
  In short, it is better to just use \code{algorithm="auto"} (the default):
  this way \code{matchPattern} will choose the algo that is best suited
  for the task.
}

\value{
  An \link{XStringViews} object for \code{matchPattern}.

  A single integer for \code{countPattern}.

  An \link{MIndex} object for \code{vmatchPattern}.

  An integer vector for \code{vcountPattern}, with each element in
  the vector corresponding to the number of matches in the corresponding
  element of \code{subject}.
}

\seealso{
  \code{\link{isMatching}},
  \code{\link{mismatch}},
  \code{\link{matchPDict}},
  \code{\link{matchLRPatterns}},
  \code{\link{matchProbePair}},
  \code{\link{maskMotif}},
  \code{\link{alphabetFrequency}},
  \link{XStringViews-class},
  \link{MIndex-class}
}

\examples{
  ## ---------------------------------------------------------------------
  ## A. matchPattern()/countPattern()
  ## ---------------------------------------------------------------------

  ## A simple inexact matching example with a short subject
  x <- DNAString("AAGCGCGATATG")
  m1 <- matchPattern("GCNNNAT", x)
  m1
  m2 <- matchPattern("GCNNNAT", x, fixed=FALSE)
  m2
  as.matrix(m2)

  ## With DNA sequence of yeast chromosome number 1
  data(yeastSEQCHR1)
  yeast1 <- DNAString(yeastSEQCHR1)
  PpiI <- "GAACNNNNNCTC" # a restriction enzyme pattern
  match1.PpiI <- matchPattern(PpiI, yeast1, fixed=FALSE)
  match2.PpiI <- matchPattern(PpiI, yeast1, max.mismatch=1, fixed=FALSE)

  ## With a genome containing isolated Ns
  library(BSgenome.Celegans.UCSC.ce2)
  chrII <- Celegans[["chrII"]]
  alphabetFrequency(chrII)
  matchPattern("N", chrII)
  matchPattern("TGGGTGTCTTT", chrII) # no match
  matchPattern("TGGGTGTCTTT", chrII, fixed=FALSE) # 1 match

  ## Using wildcards ("N") in the pattern on a genome containing N-blocks
  library(BSgenome.Dmelanogaster.UCSC.dm3)
  chrX <- maskMotif(Dmelanogaster$chrX, "N")
  as(chrX, "XStringViews") # 4 non masked regions
  matchPattern("TTTATGNTTGGTA", chrX, fixed=FALSE)
  ## Can also be achieved with no mask
  masks(chrX) <- NULL
  matchPattern("TTTATGNTTGGTA", chrX, fixed="subject")

  ## ---------------------------------------------------------------------
  ## B. vmatchPattern()/vcountPattern()
  ## ---------------------------------------------------------------------

  Ebox <- DNAString("CANNTG")
  subject <- Celegans$upstream5000
  mindex <- vmatchPattern(Ebox, subject, fixed=FALSE)
  count_index <- countIndex(mindex)    # Get the number of matches per
                                       # subject element.
  sum(count_index)                     # Total number of matches.
  table(count_index)
  i0 <- which(count_index == max(count_index))
  subject[i0]                          # The subject element with most
                                       # occurrences.
  mindex[[i0]]                         # Its matches as an IRanges object.
  Views(subject[[i0]], start=start(mindex[[i0]]), end=end(mindex[[i0]])) # And
                                       # as an XStringViews object.

  ## ---------------------------------------------------------------------
  ## C. With indels (WORK-IN-PROGRESS, DO NOT USE YET!!!)
  ## ---------------------------------------------------------------------
  library(BSgenome.Celegans.UCSC.ce2)
  pattern <- DNAString("ACGGACCTAATGTTATC")
  subject <- Celegans$chrI

  ## Allowing up to 2 mismatching letters doesn't give any match
  matchPattern(pattern, subject, max.mismatch=2)

  ## But allowing up to 2 edit operations gives 3 matches (matches 3 and 4
  ## can be considered the same and need to be merged)
  system.time(m <- matchPattern(pattern, subject, max.mismatch=2, with.indels=TRUE))
  m

  ## pairwiseAlignment() returns the (first) best match only
  mat <- nucleotideSubstitutionMatrix(match=1, mismatch=0, baseOnly=TRUE)
  system.time(pwa <- pairwiseAlignment(pattern, subject, type="local",
                       substitutionMatrix=mat, gapOpening=0, gapExtension=1))
  pwa
}

\keyword{methods}
