\name{match-utils}

\alias{match-utils}

\alias{neditStartingAt}
\alias{neditEndingAt}
\alias{neditAt}
\alias{isMatchingStartingAt}
\alias{isMatchingEndingAt}
\alias{isMatchingAt}
\alias{neditStartingAt,character-method}
\alias{neditStartingAt,XString-method}
\alias{neditStartingAt,XStringSet-method}
\alias{neditEndingAt,character-method}
\alias{neditEndingAt,XString-method}
\alias{neditEndingAt,XStringSet-method}
\alias{isMatchingStartingAt,character-method}
\alias{isMatchingStartingAt,XString-method}
\alias{isMatchingStartingAt,XStringSet-method}
\alias{isMatchingEndingAt,character-method}
\alias{isMatchingEndingAt,XString-method}
\alias{isMatchingEndingAt,XStringSet-method}

\alias{mismatch}
\alias{mismatch,ANY,XStringViews-method}
\alias{nmatch}
\alias{nmatch,ANY,XStringViews-method}
\alias{nmismatch}
\alias{nmismatch,ANY,XStringViews-method}

\alias{coverage,XStringViews-method}
\alias{coverage,MaskedXString-method}
\alias{coverage,MIndex-method}

% Old stuff:
\alias{nmismatchStartingAt}
\alias{nmismatchEndingAt}
\alias{isMatching}


\title{Utility functions related to pattern matching}

\description{
  In this man page we define precisely and illustrate what a "match" of a
  pattern P in a subject S is in the context of the Biostrings package.
  This definition of a "match" is central to most pattern matching functions
  available in this package: unless specified otherwise, most of them will
  adhere to the definition provided here.

  \code{neditStartingAt}, \code{neditEndingAt}, \code{isMatchingStartingAt}
  and \code{isMatchingEndingAt} are low-level functions that implement
  some basic concepts. Once these concepts are understood, we can use them
  to provide a simple and concise definition of a "match".

  Other utility functions related to pattern matching are described here:
  the \code{mismatch} function for getting the positions of the mismatching
  letters of a given pattern relatively to its matches in a given subject,
  the \code{nmatch} and \code{nmismatch} functions for getting the number of
  matching and mismatching letters produced by the \code{mismatch} function,
  and the \code{coverage} function that can be used to get the "coverage" of
  a subject by a given pattern or set of patterns.
}

\usage{
  neditStartingAt(pattern, subject, starting.at=1, with.indels=FALSE, fixed=TRUE)
  neditEndingAt(pattern, subject, ending.at=1, with.indels=FALSE, fixed=TRUE)
  neditAt(pattern, subject, at=1, with.indels=FALSE, fixed=TRUE)

  isMatchingStartingAt(pattern, subject, starting.at=1,
                  max.mismatch=0, with.indels=FALSE, fixed=TRUE)
  isMatchingEndingAt(pattern, subject, ending.at=1,
                  max.mismatch=0, with.indels=FALSE, fixed=TRUE)
  isMatchingAt(pattern, subject, at=1,
                  max.mismatch=0, with.indels=FALSE, fixed=TRUE)

  mismatch(pattern, x, fixed=TRUE)
  nmatch(pattern, x, fixed=TRUE)
  nmismatch(pattern, x, fixed=TRUE)
  \S4method{coverage}{MIndex}(x, start=NA, end=NA)
  \S4method{coverage}{XStringViews}(x, start=NA, end=NA, weight=1L)
  \S4method{coverage}{MaskedXString}(x, start=NA, end=NA, weight=1L)
}

\arguments{
  \item{pattern}{
    The pattern string.
  }
  \item{subject}{
    An \link{XString}, \link{XStringSet} object, or character vector
    containing the subject sequence.
  }
  \item{starting.at, ending.at, at}{
    An integer vector specifying the starting (for \code{starting.at}
    and \code{at}) or ending (for \code{ending.at}) positions of the
    pattern relatively to the subject.
  }
  \item{max.mismatch}{
    See details below.
  }
  \item{with.indels}{
    See details below.
  }
  \item{fixed}{
    Only with a \link{DNAString} or \link{RNAString}-based subject can a
    \code{fixed} value other than the default (\code{TRUE}) be used.

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
    but \link[IRanges]{IRanges}, \link[IRanges]{MaskCollection} and
    \link{MaskedXString} objects are accepted too.
  }
  \item{start, end}{
    Two single integers specifying where to start and end the extraction of the
    coverage in \code{x}.
  }
  \item{weight}{
    An integer vector specifying how much each element in \code{x} counts.
  }
}

\details{
  A "match" of pattern P in subject S is a substring S' of S that is considered
  similar enough to P according to some distance (or metric) specified by the
  user. 2 distances are supported by most pattern matching functions in the
  Biostrings package. The first (and simplest) one is the "number of mismatching
  letters". It is defined only when the 2 strings to compare have the same
  length, so when this distance is used, only matches that have the same number
  of letters as P are considered.
  The second one is the "edit distance" (aka Levenshtein distance): it's
  the minimum number of operations needed to transform P into S', where an
  operation is an insertion, deletion, or substitution of a single letter.
  When this metric is used, matches can have a different number of letters
  than P.

  The \code{neditStartingAt} (and \code{neditEndingAt}) function implements
  these 2 distances.
  If \code{with.indels} is \code{FALSE} (the default), then the first distance
  is used i.e. \code{neditStartingAt} returns the "number of mismatching
  letters" between the pattern P and the substring S' of S starting at the
  positions specified in \code{starting.at} (note that \code{neditStartingAt}
  and \code{neditEndingAt} are vectorized so long vectors of integers can be
  passed thru the \code{starting.at} or \code{ending.at} arguments).
  If \code{with.indels} is \code{TRUE}, then the "edit distance" distance is
  used: for each position specified in \code{starting.at}, P is compared to
  all the substrings S' of S starting at this position and the smallest
  distance is returned. Note that this distance is guaranteed to be reached
  for a substrings of length < 2*length(P) so, of course, in practise,
  P only needs to be compared to a small number of substrings for every
  starting position.
}

\value{
  \code{neditStartingAt} and \code{neditEndingAt}: If \code{subject} is an
  \link{XString}, then returns an integer vector of the same length as
  \code{starting.at} (or \code{ending.at}). If \code{subject} is an
  \link{XStringSet}, then returns an integer matrix with
  \code{length(subject)} rows and \code{length(starting.at)}
  (or \code{length(ending.at)}) columns.

  \code{isMatchingStartingAt(...)} and \code{isMatchingEndingAt(...)}:
  If \code{subject} is an \link{XString}, then returns a logical vector
  defined by \code{neditStartingAt(...) <= max.mismatch} or
  \code{neditEndingAt(...) <= max.mismatch}, respectively.
  If \code{subject} is an \link{XStringSet}, then returns a logical matrix
  with \code{length(subject)} rows and \code{length(starting.at)}
  (or \code{length(ending.at)}) columns.
  \code{neditAt} and \code{isMatchingAt} are conveniency wrappers for
  \code{neditStartingAt} and \code{isMatchingStartingAt}, respectively.

  \code{mismatch}:  a list of integer vectors.
  
  \code{nmismatch}:  an integer vector containing the length of the vectors
  produced by \code{mismatch}.

  \code{coverage}:  an \link[IRanges]{Rle} object indicating the coverage of
  \code{x} in the interval specified by the \code{start} and \code{end}
  arguments.
  An integer value called the "coverage" can be associated to each position
  in \code{x}, indicating how many times this position is covered by the views
  or matches stored in \code{x}. For example, if \code{x} is an
  \link{XStringViews} object, the coverage of a given position in \code{x} is
  the number of views it belongs to.
  If \code{x} is an \link{MIndex} object, the coverage of a given position
  in \code{x} is the number of matches (or hits) it belongs to.
  Note that the positions in the returned \link[IRanges]{Rle} object are to be
  interpreted as relative to the interval specified by the \code{start} and
  \code{end} arguments.
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
  \link{align-utils}
}

\examples{
  ## ---------------------------------------------------------------------
  ## neditAt() / isMatchingAt()
  ## ---------------------------------------------------------------------
  subject <- DNAString("GTATA")

  ## Pattern "AT" matches subject "GTATA" at position 3 (exact match)
  neditAt("AT", subject, at=3)
  isMatchingAt("AT", subject, at=3)

  ## ... but not at position 1
  neditAt("AT", subject)
  isMatchingAt("AT", subject)

  ## ... unless we allow 1 mismatching letter (inexact match)
  isMatchingAt("AT", subject, max.mismatch=1)

  ## Here we look at 6 different starting positions and find 3 matches if
  ## we allow 1 mismatching letter
  isMatchingAt("AT", subject, at=0:5, max.mismatch=1)

  ## No match
  neditAt("NT", subject, at=1:4)
  isMatchingAt("NT", subject, at=1:4)

  ## 2 matches if N is interpreted as an ambiguity (fixed=FALSE)
  neditAt("NT", subject, at=1:4, fixed=FALSE)
  isMatchingAt("NT", subject, at=1:4, fixed=FALSE)

  ## max.mismatch != 0 and fixed=FALSE can be used together
  neditAt("NCA", subject, at=0:5, fixed=FALSE)
  isMatchingAt("NCA", subject, at=0:5, max.mismatch=1, fixed=FALSE)

  some_starts <- c(10:-10, NA, 6)
  subject <- DNAString("ACGTGCA")
  is_matching <- isMatchingAt("CAT", subject, at=some_starts, max.mismatch=1)
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
