\name{matchPattern}

\alias{matchPattern}
\alias{matchPattern,ANY,character-method}
\alias{matchPattern,ANY,BString-method}
\alias{matchPattern,ANY,BStringViews-method}
\alias{countPattern}
\alias{countPattern,ANY,character-method}
\alias{countPattern,ANY,BString-method}
\alias{countPattern,ANY,BStringViews-method}
\alias{mismatch}
\alias{mismatch,ANY,BStringViews-method}

\alias{matchDNAPattern}


\title{String searching functions}

\description{
  Generic that finds all matches of a pattern in a BString.
  Currently only the ``shift-or'' algorithm is available. 
}

\usage{
  matchPattern(pattern, subject, algorithm="shift-or", mismatch=0, fixed=NA)
  countPattern(pattern, subject, algorithm="shift-or", mismatch=0, fixed=NA)
  mismatch(pattern, x, fixed=NA)
}

\arguments{
  \item{pattern}{
    The pattern string.
  }
  \item{subject}{
    A \link{BString} (or \link{DNAString} or \link{RNAString}) object
    containing the subject string.
  }
  \item{algorithm}{
    One of the following: \code{"boyer-moore"}, \code{"forward-search"}
    or \code{"shift-or"} (currently, only \code{"shift-or"} is available).
  }
  \item{mismatch}{
    The number of mismatches allowed. If non-zero a fuzzy string searching
    algorithm is used for matching.
  }
  \item{fixed}{
    A logical.  If \code{TRUE}, \code{pattern} is a string to be matched as is.
  }
  \item{x}{
    A \link{BStringViews} object (typically, one returned
    by code{matchPattern(pattern, subject)}).
  }
}

\details{
}

\value{
}

\author{H. Pages}

\seealso{}

\examples{
}

\keyword{classes}

