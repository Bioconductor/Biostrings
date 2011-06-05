\name{getSeq}

\alias{getSeq}

\title{getSeq}

\description{
   A generic function for extracting a set of sequences (or subsequences)
   from a sequence container like a \link[BSgenome]{BSgenome} object or
   other.
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
  An \link{XString} object or an \link{XStringSet} object or a character
  vector containing the extracted sequence(s).

  See man pages of individual methods for the details e.g. with
  \code{?`\link[BSgenome]{getSeq,BSgenome-method}`} to access the man page
  of the method for \link[BSgenome]{BSgenome} objects (make sure the
  BSgenome package is loaded first). 
}

\seealso{
  \link[BSgenome]{getSeq,BSgenome-method},
  \link{XString-class},
  \link{XStringSet-class}
}

\examples{
## Note that you need to load the package(s) defining the specialized
## methods to have showMethods() display them and to be able to access
## their man pages:
library(BSgenome)
showMethods("getSeq")
}

\keyword{manip}
