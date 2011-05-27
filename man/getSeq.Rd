\name{getSeq}

\alias{getSeq}

\title{getSeq}

\description{
   A generic function for extracting a set of sequences (or subsequences)
   from a sequence container like a \link[BSgenome]{BSgenome} object or
   others.
}

\usage{
getSeq(x, ...)
}

\arguments{
  \item{x}{
    A \link[BSgenome]{BSgenome} object or any other supported object.
    Do \code{showMethods("getSeq")} to get the list of all supported
    types for \code{x}. 
  }
  \item{...}{
    Any additional arguments needed by the specialized methods.
  }
}

\value{
  An \link{XStringSet} object or a character vector containing the
  extracted sequences. See individual methods for the details.
}

\seealso{
  \link[BSgenome]{getSeq,BSgenome-method},
  \link{XStringSet-class}
}

\examples{
## Note that you need to load the package(s) defining the specialized
## method(s) to have showMethods() display them:
library(BSgenome)
showMethods("getSeq")
}

\keyword{manip}
