\name{matchPattern}

\alias{matchPattern}
\alias{matchPattern,character-method}
\alias{matchPattern,BString-method}
\alias{matchPattern,BStringViews-method}
\alias{countPattern}
\alias{countPattern,character-method}
\alias{countPattern,BString-method}
\alias{countPattern,BStringViews-method}
\alias{mismatch}
\alias{mismatch,BStringViews-method}

\alias{matchDNAPattern}


\title{String searching functions}

\description{
  Generic that finds all matches of a pattern in a BString.
  Currently only a ``Boyer-Moore-like'' algorithm and the ``shift-or''
  algorithm are available. 
}

\usage{
  matchPattern(pattern, subject, algorithm="auto", mismatch=0, fixed=TRUE)
  countPattern(pattern, subject, algorithm="auto", mismatch=0, fixed=TRUE)
  mismatch(pattern, x, fixed=TRUE)
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
    One of the following: \code{"auto"}, \code{"naive"}, \code{"boyer-moore"},
    \code{"forward-search"} or \code{"shift-or"} (note that
    \code{"forward-search"} is not available yet).
  }
  \item{mismatch}{
    The number of mismatches allowed. If non-zero a fuzzy string searching
    algorithm is used for matching.
  }
  \item{fixed}{
    A logical. Only \link{DNAString} or \link{RNAString} support \code{fixed=FALSE}.
    In this case, ambiguities (i.e. letters from the IUPAC Extended Genetic Alphabet
    (see \code{\link{IUPAC_CODE_MAP}}) that are not from the base alphabet)
    match any letter that they stand for.
  }
  \item{x}{
    A \link{BStringViews} object (typically, one returned
    by \code{matchPattern(pattern, subject)}).
  }
}

\details{
}

\value{
}

\author{H. Pages}

\seealso{
  
}

\examples{
}

\keyword{methods}

