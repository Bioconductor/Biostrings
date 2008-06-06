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
\alias{mismatch,XStringViews-method}

\alias{coverage}
\alias{coverage,IRanges-method}
\alias{coverage,MaskCollection-method}
\alias{coverage,XStringViews-method}
\alias{coverage,MaskedXString-method}
\alias{coverage,MIndex-method}


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
  and the \code{coverage} function that can be used to get the "coverage"
  of a subject by a given pattern or set of patterns.
}

\usage{
  nmismatchStartingAt(pattern, subject, starting.at=1, fixed=TRUE)
  nmismatchEndingAt(pattern, subject, ending.at=1, fixed=TRUE)
  isMatching(pattern, subject, start=1, max.mismatch=0, fixed=TRUE)
  mismatch(pattern, x, fixed=TRUE)
  coverage(x, start=NA, end=NA)
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
    but \link{IRanges}, \link{MaskCollection} and \link{MaskedXString} objects
    are accepted too.
  }
  \item{end}{
    A single integer specifying the position in \code{x}
    where to end the extraction of the coverage.
  }
}

\value{
  \code{nmismatchStartingAt}, \code{nmismatchEndingAt}: an integer vector of
  the same length as \code{starting.at} (or \code{ending.at}) reporting the
  number of mismatching letters for each starting (or ending) position.

  A logical vector of the same length as \code{start} for \code{isMatching}.

  A list of integer vectors for \code{mismatch}.

  An integer vector indicating the coverage of \code{x} in the interval
  specified by the \code{start} and \code{end} arguments.
  An integer value called the "coverage" can be associated to each position
  in \code{x}, indicating how many times this position is covered by the views
  or matches stored in \code{x}. For example, if \code{x} is an
  \link{XStringViews} object, the coverage of a given position in \code{x} is
  the number of views it belongs to.
  If \code{x} is an \link{MIndex} object, the coverage of a given position
  in \code{x} is the number of matches (or hits) it belongs to.
  Note that the positions in the returned vector are to be interpreted as
  relative to the interval specified by the \code{start} and \code{end}
  arguments.

}

\seealso{
  \code{\link{matchPattern}},
  \code{\link{matchPDict}},
  \code{\link{IUPAC_CODE_MAP}},
  \link{XString-class},
  \link{XStringViews-class},
  \link{MIndex-class},
  \link{IRanges-class},
  \link{MaskCollection-class},
  \link{MaskedXString-class}
}

\examples{
  ## ---------------------------------------------------------------------
  ## isMatching()
  ## ---------------------------------------------------------------------
  subject <- DNAString("GTATA")
  ## Pattern "AT" matches subject "GTATA" at position 3
  isMatching("AT", subject, start=3)
  ## No match at position 1
  isMatching("AT", subject)
  ## Match at position 1 if we allow 1 mismatching letter
  isMatching("AT", subject, max.mismatch=1)
  ## 3 matches in total if we allow 1 mismatching letter
  isMatching("AT", subject, start=0:5, max.mismatch=1)
  ## No match
  isMatching("NT", subject, start=1:4)
  ## 2 matches if N is interpreted as an ambiguity (fixed=FALSE)
  isMatching("NT", subject, start=1:4, fixed=FALSE)
  ## max.mismatch != 0 and fixed=FALSE can be used together
  isMatching("NCA", subject, start=0:5, max.mismatch=1, fixed=FALSE)

  some_starts <- c(10:-10, NA, 6)
  is_matching <- isMatching("CAT", DNAString("ACGTGCA"), start=some_starts, max.mismatch=1)
  some_starts[is_matching]

  ## ---------------------------------------------------------------------
  ## mismatch()
  ## ---------------------------------------------------------------------
  m <- matchPattern("NCA", subject, max.mismatch=1, fixed=FALSE)
  mismatch("NCA", m)

  ## ---------------------------------------------------------------------
  ## coverage()
  ## ---------------------------------------------------------------------
  coverage(m)

  x <- IRanges(start=c(-2L, 6L, 9L, -4L, 1L, 0L, -6L, 10L),
               width=c( 5L, 0L, 6L,  1L, 4L, 3L,  2L,  3L))
  coverage(x, start=-6, end=20)  # 'start' and 'end' must be specified for
                                 # an IRanges object.
  coverage(shift(x, 2), start=-6, end=20)
  coverage(restrict(x, 1, 10), start=-6, end=20)
  coverage(reduce(x), start=-6, end=20)
  coverage(gaps(x, start=-6, end=20), start=-6, end=20)

  mask1 <- Mask(mask.width=29, start=c(11, 25, 28), width=c(5, 2, 2))
  mask2 <- Mask(mask.width=29, start=c(3, 10, 27), width=c(5, 8, 1))
  mask3 <- Mask(mask.width=29, start=c(7, 12), width=c(2, 4))
  mymasks <- append(append(mask1, mask2), mask3)
  coverage(mymasks)

  ## See ?matchPDict for examples of using coverage() on an MIndex object...
}

\keyword{methods}

