\name{AAString-class}
\docType{class}

% Classes
\alias{class:AAString}
\alias{AAString-class}

% Constructor-like functions and generics:
\alias{initialize,AAString-method}
\alias{AAString}


\title{The AAString class}

\description{
  A \code{AAString} object allows efficient storage
  and manipulation of a long amino acid sequence.
}

\details{
  The \code{AAString} class derives directly from the
  \code{\link{BString}} class (with no additional slot).
  All functions and methods described in the
  \code{\link{BString}} man page also work
  with a \code{AAString} object (inheritance).

  Unlike a \code{\link{BString}} object that allows storage
  of any non-empty string (based on a single-byte character set)
  a \code{AAString} object can only store a non-empty string
  based on a limited alphabet (the Amino Acid
  alphabet + the gap letter).
}

\section{The Amino Acid alphabet}{
  [TODO]
}

\section{Constructor-like functions and generics}{
  In the code snippet below,
  \code{src} can be a character vector
  or a \code{\link{BString}} (or derived) object.

  \describe{
    \item{}{
      \code{AAString(src)}:
      [TODO: Document me]
    }
  }
}

\author{H. Pages}

\seealso{
  \link{BString-class}, \link{BStringViews-class}
}

\examples{
  ## [TODO: Provide examples]
}

\keyword{methods}
\keyword{classes}
