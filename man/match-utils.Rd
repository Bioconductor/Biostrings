\name{match-utils}

\alias{isMatching}
\alias{isMatching,character-method}
\alias{isMatching,BString-method}
\alias{mismatch}
\alias{mismatch,BStringViews-method}


\title{Utility functions related to pattern matching}

\description{
  Some Utility functions related to pattern matching.
}

\usage{
  isMatching(pattern, subject, start=1, max.mismatch=0, fixed=TRUE)
  mismatch(pattern, x, fixed=TRUE)
}

\arguments{
  \item{pattern}{
    The pattern string.
  }
  \item{subject}{
    A \link{BString} (or derived) object containing the subject string,
  }
  \item{start}{
    Integer vector specifying where to place the pattern relatively to the subject.
  }
  \item{max.mismatch}{
    The maximum number of mismatches allowed. If non-zero, an inexact matching
    algorithm is used.
    Note that \code{matchPattern} (and \code{countPattern}) don't support the
    kind of inexact matching where a given number of insertions or deletions
    are allowed. Hence the matches found have always the same length (called
    the "width" in this particular context) as the pattern.
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
    A \link{BStringViews} object (typically, one returned
    by \code{matchPattern(pattern, subject)}).
  }
}

\value{
  A logical vector of the same length as \code{start} for \code{isMatching}.

  A list of integer vectors for \code{mismatch}.
}

\seealso{
  \code{\link{matchPattern}},
  \code{\link{matchPDict}},
  \code{\link{IUPAC_CODE_MAP}},
  \link{BStringViews-class},
  \link{DNAString-class}
}

\examples{
  isMatching("GC", "GCTGT")
  isMatching("GC", "GCTGT", start=4)
  isMatching("GC", "GCTGT", start=4, max.mismatch=1)
  isMatching("GN", "GCTGT", start=1:4, fixed=TRUE)
  isMatching("GN", "GCTGT", start=0, max.mismatch=1, fixed=TRUE)
}

\keyword{methods}

