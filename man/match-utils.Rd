\name{match-utils}

\alias{isMatching}
\alias{isMatching,character-method}
\alias{isMatching,XString-method}
\alias{mismatch}
\alias{mismatch,XStringViews-method}


\title{Utility functions related to pattern matching}

\description{
  This man page gives some background information about the concept of "match"
  ("exact match" or "inexact match") as understood by the various pattern matching
  functions available in Biostrings.

  The \code{isMatching} function implements this concept.

  Some other utility functions related to pattern matching might be added here
  in the future.
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
    An \link{XString} object containing the subject sequence,
  }
  \item{start}{
    Integer vector specifying where to place the pattern relatively to the subject.
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
    An \link{XStringViews} object (typically, one returned
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
  \link{XString-class},
  \link{XStringViews-class}
}

\examples{
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
  isMatching("NGA", subject, start=0:5, max.mismatch=1, fixed=FALSE)

  some_starts <- c(10:-10, NA, 6)
  is_matching <- isMatching("CAT", DNAString("ACGTGCA"), start=some_starts, max.mismatch=1)
  some_starts[is_matching]
}

\keyword{methods}

