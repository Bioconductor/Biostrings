\name{RNAString-class}
\docType{class}

% Classes
\alias{class:RNAString}
\alias{RNAString-class}

% Constants:
\alias{RNA_ALPHABET}

% Constructor-like functions and generics:
\alias{initialize,RNAString-method}
\alias{RNAString}

% Other functions and generics:
\alias{alphabet,RNAString-method}


\title{The RNAString class}

\description{
  A \code{RNAString} object allows efficient storage
  and manipulation of a long RNA sequence.
}

\details{
  The \code{RNAString} class derives directly from the
  \code{\link{BString}} class (with no additional slot).
  All functions and methods described in the
  \code{\link{BString}} man page also work
  with a \code{RNAString} object (inheritance).

  Unlike a \code{\link{BString}} object that allows storage
  of any non-empty string (based on a single-byte character set)
  a \code{RNAString} object can only store a non-empty string
  based on a limited alphabet (the IUPAC Extended Genetic
  alphabet + the gap letter).
  In addition, the letters stored in a \code{RNAString} object
  are encoded in a way that optimizes fast search algorithms.
}

\section{The IUPAC Extended Genetic alphabet}{

  [TODO]

  \describe{
    \item{}{
      \code{RNA_ALPHABET}:
      [TODO: Document me]
    }
  }
}

\section{Constructor-like functions and generics}{
  In the code snippet below,
  \code{src} can be a character vector
  or a \code{\link{BString}} (or derived) object.

  \describe{
    \item{}{
      \code{RNAString(src)}:
      [TODO: Document me]
    }
  }
}

\section{Other functions and generics}{
  In the code snippet below,
  \code{x} is a \code{RNAString} object.

  \describe{
    \item{}{
      \code{alphabet(x)}:
      [TODO: Document me]
    }
  }
}

\author{H. Pages}

\seealso{
  \link{BString-class}, \link{DNAString-class}, \link{BStringViews-class}
}

\examples{
  d <- DNAString("TTGAAAA-CTC-N")
  length(d)
  r <- RNAString(d)
}

\keyword{methods}
\keyword{classes}
