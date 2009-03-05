\name{subXString}

\alias{subXString}
\alias{subBString}

\alias{substr,XString-method}
\alias{substr,MaskedXString-method}
\alias{substring,XString-method}
\alias{substring,MaskedXString-method}


\title{Fast substring extraction}

\description{
  Functions for fast substring extraction.
}

\usage{
  subXString(x, start=NA, end=NA, length=NA)
  \S4method{substr}{XString}(x, start=NA, stop=NA)
  \S4method{substring}{XString}(text, first=NA, last=NA)
}

\arguments{
  \item{x}{
    An \link{XString} object for \code{subXString}.
    A character vector, an \link{XStringViews}, \link{XString},
    or \link{MaskedXString} object for \code{substr} or \code{substring}.
  }
  \item{start}{
    A numeric vector.
  }
  \item{end}{
    A numeric vector.
  }
  \item{length}{
    A numeric vector.
  }
  \item{stop}{
    A numeric vector.
  }
  \item{text}{
    A character vector, an \link{XStringViews} or an \link{XString}
    object.
  }
  \item{first}{
    A numeric vector.
  }
  \item{last}{
    A numeric vector.
  }
}

\details{
  \code{subXString} is deprecated in favor of \code{\link[IRanges]{subseq}}.
}

\value{
  An \link{XString} object of the same base type as \code{x} for
  \code{subXString}.

  A character vector for \code{substr} and \code{substring}.
}

\seealso{
  \code{\link[IRanges]{subseq}},
  \code{\link{letter}},
  \link{XString-class},
  \link{XStringViews-class}
}

\keyword{methods}
\keyword{manip}

