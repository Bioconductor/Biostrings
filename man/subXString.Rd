\name{subXString}

\alias{subXString}
\alias{subBString}

\alias{substr}
\alias{substr,XString-method}
\alias{substring}
\alias{substring,XString-method}


\title{Fast substring extraction}

\description{
  Functions for fast substring extraction.
}

\usage{
  subXString(x, start=NA, end=NA, length=NA)
  substr(x, start=NA, stop=NA)
  substring(text, first=NA, last=NA)
}

\arguments{
  \item{x}{
    An \link{XString} object for \code{subXString}.
    A character vector, an \link{XStringViews} or an \link{XString}
    object for \code{substr} or \code{substring}.
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
  \code{subXString} provides a very efficient way to extract a substring
  from an \link{XString} object.
  For example, extracting a 100Mb substring from Human chromosome 1 (250Mb)
  with \code{subXString} is almost instantaneous and has almost no memory
  footprint. In fact, the cost in time and memory of a call to \code{subXString}
  doesn't depend on the length of the original object or on the length of
  the extracted object.
  The "trick" behind this "feature" is very simple: \code{subXString}
  does NOT copy the sequence data.
}

\value{
  An \link{XString} object of the same subtype as \code{x} for
  \code{subXString}.

  A character vector for \code{substr} and \code{substring}.
}

\author{H. Pages}

\seealso{
  \code{\link{subseq}},
  \code{\link{letter}},
  \code{\link{views,XString-method}},
  \link{XString-class},
  \link{XStringViews-class}
}

\examples{
  subXString("AxyxyxBC", 7)
  s <- BString("AxyxyxBC")
  subXString(s, 7)  # same as subXString("AxyxyxBC", 7)
  subXString(s, , 7)
  subXString(s, , 7, 3)
  identical(subXString(s), s) # TRUE

  v <- views(s, c(6, 4,-1, NA), c(7, 6, 1, 1))
  ## 2 equivalent ways of keeping the last letter of each view
  views(subject(v), end(v), end(v))

  ## 2 equivalent ways of making the views wider
  views(subject(v), start(v)-3, end(v)+3)
}

\keyword{methods}
\keyword{manip}

