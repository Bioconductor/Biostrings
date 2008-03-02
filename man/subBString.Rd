\name{subBString}

\alias{subBString}

\alias{substr}
\alias{substr,BString-method}
\alias{substring}
\alias{substring,BString-method}


\title{Fast substring extraction}

\description{
  Functions for fast substring extraction.
}

\usage{
  subBString(x, start=NA, end=NA, length=NA)
  substr(x, start=NA, stop=NA)
  substring(text, first=NA, last=NA)
}

\arguments{
  \item{x}{
    A \code{\link{BString}} (or derived) object for \code{subBString}.
    A character vector, a \code{\link{BStringViews}} or a \code{\link{BString}}
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
    A character vector, a \code{\link{BStringViews}} or a \code{\link{BString}}
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
  \code{subBString} provides a very efficient way to extract a substring
  from a \link{BString} (or derived) object.
  For example, extracting a 100Mb substring from Human chromosome 1 (250Mb)
  with \code{subBString} is almost instantaneous and has almost no memory
  footprint. In fact, the cost in time and memory of a call to \code{subBString}
  doesn't depend on the length of the original object or on the length of
  the extracted object.
  The "trick" behind this "feature" is very simple: \code{subBString}
  does NOT copy the sequence data.
}

\value{
  [TODO]
}

\author{H. Pages}

\seealso{
  \code{\link{letter}},
  \code{\link{views}},
  \link{BString-class},
  \link{DNAString-class},
  \link{RNAString-class},
  \link{AAString-class},
  \link{BStringViews-class}
}

\examples{
  subBString("AxyxyxBC", 7)
  s <- BString("AxyxyxBC")
  subBString(s, 7)  # same as subBString("AxyxyxBC", 7)
  subBString(s, , 7)
  subBString(s, , 7, 3)
  identical(subBString(s), s) # TRUE

  v <- views(s, c(6, 4,-1, NA), c(7, 6, 1, 1))
  ## 2 equivalent ways of keeping the last letter of each view
  views(subject(v), end(v), end(v))

  ## 2 equivalent ways of making the views wider
  views(subject(v), start(v)-3, end(v)+3)
}

\keyword{methods}
\keyword{manip}

