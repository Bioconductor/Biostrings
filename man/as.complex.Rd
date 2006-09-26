\name{as.complex}

% Standard generic methods:
\alias{as.complex}
\alias{as.complex,DNAString-method}

\title{Converting a DNAString to a vector of complexes}

\description{
}

\usage{
  as.complex(x, ...)
}

\arguments{
  \item{x}{
    A \code{\link{DNAString}} object.
  }
  \item{...}{
    Accepts an extra argument \code{baseValues} that must be
    a named complex vector containing the values associated
    to each base e.g.
      \code{c("A"=1+0i, "G"=0+1i, "T"=-1+0i, "C"=0-1i)}
  }
}

\details{
}

\author{H. Pages}

\seealso{
  \link{DNAString}
}

\examples{
}

\keyword{methods}

