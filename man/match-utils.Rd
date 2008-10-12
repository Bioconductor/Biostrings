\name{match-utils}

\alias{match-utils}

\alias{nmismatchStartingAt}
\alias{nmismatchEndingAt}
\alias{nmismatchStartingAt,character-method}
\alias{nmismatchEndingAt,character-method}
\alias{nmismatchStartingAt,XString-method}
\alias{nmismatchEndingAt,XString-method}

\alias{isMatching}
\alias{isMatching,character-method}
\alias{isMatching,XString-method}

\alias{mismatch}
\alias{mismatch,ANY,XStringViews-method}
\alias{mismatch,AlignedXStringSet,missing-method}
\alias{nmatch}
\alias{nmatch,ANY,XStringViews-method}
\alias{nmatch,PairwiseAlignment,missing-method}
\alias{nmatch,PairwiseAlignmentSummary,missing-method}
\alias{nmismatch}
\alias{nmismatch,ANY,XStringViews-method}
\alias{nmismatch,AlignedXStringSet,missing-method}
\alias{nmismatch,PairwiseAlignment,missing-method}
\alias{nmismatch,PairwiseAlignmentSummary,missing-method}

\alias{mismatchTable}
\alias{mismatchTable,AlignedXStringSet-method}
\alias{mismatchTable,QualityAlignedXStringSet-method}
\alias{mismatchTable,PairwiseAlignment-method}
\alias{mismatchSummary}
\alias{mismatchSummary,AlignedXStringSet-method}
\alias{mismatchSummary,QualityAlignedXStringSet-method}
\alias{mismatchSummary,PairwiseAlignment-method}
\alias{mismatchSummary,PairwiseAlignmentSummary-method}

\alias{coverage,XStringViews-method}
\alias{coverage,MaskedXString-method}
\alias{coverage,MIndex-method}
\alias{coverage,AlignedXStringSet-method}
\alias{coverage,PairwiseAlignment-method}
\alias{coverage,PairwiseAlignmentSummary-method}


\title{Utility functions related to pattern matching}

\description{
  This man page gives some background information about the concept of "match"
  ("exact match" or "inexact match") as understood by the various pattern
  matching functions available in the Biostrings package.

  The \code{nmismatchStartingAt}, \code{nmismatchEndingAt} and \code{isMatching}
  functions implement this concept.

  Other utility functions related to pattern matching are described here:
  the \code{mismatch} function for getting the positions of the mismatching
  letters of a given pattern relatively to its matches in a given subject,
  the \code{nmatch} and \code{nmismatch} functions for getting the number of
  matching and mismatching letters produced by the \code{mismatch} function,
  the \code{mismatchTable} function for getting the more than just the mismatching
  positions, and the \code{mismatchSummary} function for getting summaries of
  information produced by the \code{mismatchTable} function, and the \code{coverage}
  function that can be used to get the "coverage" of a subject by a given pattern or
  set of patterns.
}

\usage{
  nmismatchStartingAt(pattern, subject, starting.at=1, fixed=TRUE)
  nmismatchEndingAt(pattern, subject, ending.at=1, fixed=TRUE)
  isMatching(pattern, subject, start=1, max.mismatch=0, fixed=TRUE)
  mismatch(pattern, x, fixed=TRUE)
  nmatch(pattern, x, fixed=TRUE)
  nmismatch(pattern, x, fixed=TRUE)
  mismatchTable(x, shiftLeft=0L, shiftRight=0L, \dots)
  mismatchSummary(x, \dots)
  \S4method{mismatchSummary}{PairwiseAlignment}(x, weight=1L)
  \S4method{coverage}{MIndex}(x, start=NA, end=NA)
  \S4method{coverage}{XStringViews}(x, start=NA, end=NA, weight=1L)
  \S4method{coverage}{MaskedXString}(x, start=NA, end=NA, weight=1L)
  \S4method{coverage}{AlignedXStringSet}(x, start=NA, end=NA, weight=1L)
  \S4method{coverage}{PairwiseAlignment}(x, start=NA, end=NA, weight=1L)
}

\arguments{
  \item{pattern}{
    The pattern string.
  }
  \item{subject}{
    An \link{XString} object (or character vector) containing the subject
    sequence,
  }
  \item{starting.at}{
    An integer vector specifying the starting positions of the pattern
    relatively to the subject.
  }
  \item{ending.at}{
    An integer vector specifying the ending positions of the pattern
    relatively to the subject.
  }
  \item{start}{
    For \code{isMatching}: an integer vector specifying the starting positions
    of the pattern relatively to the subject.
    For \code{coverage}: a single integer specifying the position in \code{x}
    where to start the extraction of the coverage.
  }
  \item{max.mismatch}{
    The maximum number of mismatching letters allowed.
    Note that \code{isMatching} doesn't support the kind of inexact matching
    where a given number of insertions or deletions are allowed.
    Therefore all the "matches" (i.e. the substrings in the subject that match
    the pattern) have the length of the pattern.
  }
  \item{fixed}{
    Only with a \link{DNAString} or \link{RNAString} subject can a \code{fixed}
    value other than the default (\code{TRUE}) be used.

    With \code{fixed=FALSE}, ambiguities (i.e. letters from the IUPAC Extended
    Genetic Alphabet (see \code{\link{IUPAC_CODE_MAP}}) that are not from the
    base alphabet) in the pattern \_and\_ in the subject are interpreted as
    wildcards i.e. they match any letter that they stand for.

    \code{fixed} can also be a character vector, a subset
    of \code{c("pattern", "subject")}.
    \code{fixed=c("pattern", "subject")} is equivalent to \code{fixed=TRUE}
    (the default).
    An empty vector is equivalent to \code{fixed=FALSE}.
    With \code{fixed="subject"}, ambiguities in the pattern only
    are interpreted as wildcards.
    With \code{fixed="pattern"}, ambiguities in the subject only
    are interpreted as wildcards.
  }
  \item{x}{
    An \link{XStringViews} object for \code{mismatch} (typically, one returned
    by \code{matchPattern(pattern, subject)}).

    Typically an \link{XStringViews} or \link{MIndex} object for \code{coverage}
    but \link[IRanges]{IRanges}, \link[IRanges]{MaskCollection} and \link{MaskedXString}
    objects are accepted too.
  }
  \item{shiftLeft, shiftRight}{
    Non-positive and non-negative integers respectively that specify how many
    preceding and succeeding characters to and from the mismatch position to
    include in the mismatch substrings.
  }
  \item{weight}{
    An integer vector specifying how much each element in \code{x} counts.
  }
  \item{end}{
    A single integer specifying the position in \code{x} where to end the
    extraction of the coverage.
  }
  \item{\dots}{
    Further arguments to be passed to or from other methods.
  }
}

\value{
  \code{nmismatchStartingAt}, \code{nmismatchEndingAt}:  an integer vector of
  the same length as \code{starting.at} (or \code{ending.at}) reporting the
  number of mismatching letters for each starting (or ending) position.

  \code{isMatching}:  a logical vector of the same length as \code{start}.

  \code{mismatch}:  a list of integer vectors.
  
  \code{nmismatch}:  an integer vector containing the length of the vectors
  produced by \code{mismatch}.

  \code{mismatchTable}:  a data.frame containing the positions and substrings
  of the mismatches for the \code{AlignedXStringSet} or \code{PairwiseAlignment}
  object.

  \code{mismatchSummary}:  a list of data.frame objects containing counts and
  frequencies of the mismatches for the \code{AlignedXStringSet} or
  \code{PairwiseAlignment} object.

  \code{coverage}:  an \link[IRanges]{XRleInteger} object indicating the coverage
  of \code{x} in the interval specified by the \code{start} and \code{end}
  arguments.
  An integer value called the "coverage" can be associated to each position
  in \code{x}, indicating how many times this position is covered by the views
  or matches stored in \code{x}. For example, if \code{x} is an
  \link{XStringViews} object, the coverage of a given position in \code{x} is
  the number of views it belongs to.
  If \code{x} is an \link{MIndex} object, the coverage of a given position
  in \code{x} is the number of matches (or hits) it belongs to.
  Note that the positions in the returned \link[IRanges]{XRleInteger} object are
  to be interpreted as relative to the interval specified by the \code{start}
  and \code{end} arguments.
}

\seealso{
  \code{\link{matchPattern}},
  \code{\link{matchPDict}},
  \code{\link{IUPAC_CODE_MAP}},
  \link{XString-class},
  \link{XStringViews-class},
  \link{MIndex-class},
  \link[IRanges]{coverage},
  \link[IRanges]{IRanges-class},
  \link[IRanges]{MaskCollection-class},
  \link{MaskedXString-class},
  \link{AlignedXStringSet-class},
  \link{PairwiseAlignment-class}
}

\examples{
  ## ---------------------------------------------------------------------
  ## nmismatchStartingAt() / isMatching()
  ## ---------------------------------------------------------------------
  subject <- DNAString("GTATA")

  ## Pattern "AT" matches subject "GTATA" at position 3 (exact match)
  nmismatchStartingAt("AT", subject, starting.at=3)
  isMatching("AT", subject, start=3)

  ## ... but not at position 1
  nmismatchStartingAt("AT", subject)
  isMatching("AT", subject)

  ## ... unless we allow 1 mismatching letter (inexact match)
  isMatching("AT", subject, max.mismatch=1)

  ## Here we look at 6 different starting positions and find 3 matches if
  ## we allow 1 mismatching letter
  isMatching("AT", subject, start=0:5, max.mismatch=1)

  ## No match
  nmismatchStartingAt("NT", subject, starting.at=1:4)
  isMatching("NT", subject, start=1:4)

  ## 2 matches if N is interpreted as an ambiguity (fixed=FALSE)
  nmismatchStartingAt("NT", subject, starting.at=1:4, fixed=FALSE)
  isMatching("NT", subject, start=1:4, fixed=FALSE)

  ## max.mismatch != 0 and fixed=FALSE can be used together
  nmismatchStartingAt("NCA", subject, starting.at=0:5, fixed=FALSE)
  isMatching("NCA", subject, start=0:5, max.mismatch=1, fixed=FALSE)

  some_starts <- c(10:-10, NA, 6)
  subject <- DNAString("ACGTGCA")
  is_matching <- isMatching("CAT", subject, start=some_starts, max.mismatch=1)
  some_starts[is_matching]

  ## ---------------------------------------------------------------------
  ## mismatch() / nmismatch()
  ## ---------------------------------------------------------------------
  m <- matchPattern("NCA", subject, max.mismatch=1, fixed=FALSE)
  mismatch("NCA", m)
  nmismatch("NCA", m)

  ## ---------------------------------------------------------------------
  ## coverage()
  ## ---------------------------------------------------------------------
  coverage(m)

  ## See ?matchPDict for examples of using coverage() on an MIndex object...
}

\keyword{methods}
