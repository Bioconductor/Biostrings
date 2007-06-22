\name{subBString}

\alias{subBString}
\alias{subBStringViews}

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
  subBStringViews(x, start=NA, end=NA, width=NA, check.limits=TRUE)
  substr(x, start=NA, stop=NA)
  substring(text, first=NA, last=NA)
}

\arguments{
  \item{x}{
    A \code{\link{BString}} (or derived) object for \code{subBString}.
    A \code{\link{BStringViews}} object for \code{subBStringViews}.
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
  [TODO]
}

\value{
  [TODO]
}

\seealso{
  \link{BString-class},
  \link{BStringViews-class},
  \code{\link{views}}
}

\examples{
  subBString("AxyxyxBC", 7)
  subBString("AxyxyxBC", , 7)
  subBString("AxyxyxBC", , 7, 3)

  v <- views("AxyxyxBC", c(6, 4,-1, NA), c(7, 6, 1, 1))
  ## keep only the last letter of each view
  subBStringViews(v, end=width(v), width=1)
  ## enlarge the views  
  subBStringViews(v, -2, end=width(v)+3, check.limits=FALSE)
}

\keyword{methods}
\keyword{manip}

