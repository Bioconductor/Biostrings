\name{match-utils}

\alias{match-utils}

\alias{mismatch}
\alias{mismatch,ANY,XStringViews-method}
\alias{nmatch}
\alias{nmatch,ANY,XStringViews-method}
\alias{nmismatch}
\alias{nmismatch,ANY,XStringViews-method}

\alias{coverage,MaskedXString-method}
\alias{coverage,MIndex-method}


\title{Utility functions operating on the matches returned by a
  high-level matching function}

\description{
  Miscellaneous utility functions operating on the matches returned by a
  high-level matching function like \code{\link{matchPattern}},
  \code{\link{matchPDict}}, etc...
}

\usage{
  mismatch(pattern, x, fixed=TRUE)
  nmatch(pattern, x, fixed=TRUE)
  nmismatch(pattern, x, fixed=TRUE)
  \S4method{coverage}{MIndex}(x, shift=0L, width=NULL, weight=1L)
  \S4method{coverage}{MaskedXString}(x, shift=0L, width=NULL, weight=1L)
}

\arguments{
  \item{pattern}{
    The pattern string.
  }
  \item{x}{
    An \link{XStringViews} object for \code{mismatch} (typically, one returned
    by \code{matchPattern(pattern, subject)}).

    An \link{MIndex} object for \code{coverage}, or any object for
    which a \code{coverage} method is defined. See \code{?\link[IRanges]{coverage}}.
  }
  \item{fixed}{
    See \code{?`\link{lowlevel-matching}`}.
  }
  \item{shift, width}{
    See \code{?\link[IRanges]{coverage}}.
  }
  \item{weight}{
    An integer vector specifying how much each element in \code{x} counts.
  }
}

\details{
  The \code{mismatch} function gives the positions of the mismatching
  letters of a given pattern relatively to its matches in a given subject.

  The \code{nmatch} and \code{nmismatch} functions give the number of
  matching and mismatching letters produced by the \code{mismatch} function.

  The \code{coverage} function computes the "coverage" of a subject
  by a given pattern or set of patterns.
}

\value{
  \code{mismatch}:  a list of integer vectors.
  
  \code{nmismatch}:  an integer vector containing the length of the vectors
  produced by \code{mismatch}.

  \code{coverage}:  an \link[IRanges:Rle-class]{Rle} object indicating
  the coverage of \code{x}.
  See \code{?\link[IRanges]{coverage}} for the details.
  If \code{x} is an \link{MIndex} object, the coverage of a given position
  in the underlying sequence (typically the subject used during the search
  that returned \code{x}) is the number of matches (or hits) it belongs to.
}

\seealso{
  \link{lowlevel-matching},
  \code{\link{matchPattern}},
  \code{\link{matchPDict}},
  \link{XString-class},
  \link{XStringViews-class},
  \link{MIndex-class},
  \link[IRanges]{coverage},
  \link{align-utils}
}

\examples{
  ## ---------------------------------------------------------------------
  ## mismatch() / nmismatch()
  ## ---------------------------------------------------------------------
  subject <- DNAString("ACGTGCA")
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
