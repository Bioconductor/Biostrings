\name{subBString}

\alias{subBString}
\alias{subviews}

\alias{substr}
\alias{substr,BString-method}
\alias{substr,BStringViews-method}
\alias{substring}
\alias{substring,BString-method}
\alias{substring,BStringViews-method}


\title{Substring extraction}

\description{
  Functions for fast substring extraction.
}

\usage{
  subBString(x, start=NA, end=NA, length=NA)
  subviews(x, start=NA, end=NA, width=NA, check.limits=TRUE)
  substr(x, start=NA, stop=NA)
  substring(text, first=NA, last=NA)
}

\arguments{
  \item{x}{
    A \code{\link{BString}} (or derived) object for \code{subBString}.
    A \code{\link{BStringViews}} object for \code{subviews}.
    Any object for which a method is defined for \code{substr}.
  }
  \item{start}{
    A numeric vector.
  }
  \item{end}{
    A numeric vector.
  }
  \item{length}{
  }
  \item{width}{
  }
  \item{check.limits}{
  }
  \item{stop}{
  }
  \item{text}{
  }
  \item{first}{
  }
  \item{last}{
  }
}

\details{
  For \code{subviews}, the start/end values must be given relatively to
  \code{x} views.
}

\value{
  [TODO]
}

\seealso{
  \code{\link{letter}},
  \code{\link{views}},
  \link{BString-class},
  \link{BStringViews-class}
}

\examples{
  subBString("AxyxyxBC", 7)
  subBString("AxyxyxBC", , 7)
  subBString("AxyxyxBC", , 7, 3)

  v <- views("AxyxyxBC", c(6, 4,-1, NA), c(7, 6, 1, 1))
  ## 2 equivalent ways of keeping the last letter of each view
  views(subject(v), end(v), end(v))
  subviews(v, end=width(v), width=1)
  
  ## 2 equivalent ways of making the views wider
  views(subject(v), start(v)-3, end(v)+3)
  subviews(v, -2, end=width(v)+3, check.limits=FALSE)
}

\keyword{methods}
\keyword{manip}

