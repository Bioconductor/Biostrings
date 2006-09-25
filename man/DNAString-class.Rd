\name{DNAString-class}
\docType{class}

% Classes
\alias{class:DNAString}
\alias{DNAString-class}

% Constants:
\alias{DNA_ALPHABET}

% Constructor-like functions and generics:
\alias{initialize,DNAString-method}
\alias{DNAString}

% Other functions and generics:
\alias{alphabet,DNAString-method}


\title{The DNAString class}

\description{
  A \code{DNAString} object allows efficient storage
  and manipulation of a long DNA sequence.
}

\details{
  The \code{DNAString} class derives directly from the
  \code{\link{BString}} class (with no additional slot).
  All functions and methods described in the
  \code{\link{BString}} man page also work
  with a \code{DNAString} object (inheritance).

  Unlike a \code{\link{BString}} object that allows storage
  of any non-empty string (based on a single-byte character set)
  a \code{DNAString} object can only store a non-empty string
  based on a limited alphabet (the IUPAC Extended Genetic
  alphabet + the gap letter).
  In addition, the letters stored in a \code{DNAString} object
  are encoded in a way that optimizes fast search algorithms.
}

\section{The IUPAC Extended Genetic alphabet}{

  [TODO]

  \describe{
    \item{}{
      \code{DNA_ALPHABET}:
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
      \code{DNAString(src)}:
      [TODO: Document me]
    }
  }
}

\section{Other functions and generics}{
  In the code snippet below,
  \code{x} is a \code{DNAString} object.

  \describe{
    \item{}{
      \code{alphabet(x)}:
      [TODO: Document me]
    }
  }
}

\author{H. Pages}

\seealso{
  \link{BString-class}, \link{RNAString-class}, \link{BStringViews-class}
}

\examples{
  d <- DNAString("TTGAAAA-CTC-N")
  length(d)
}

\keyword{methods}
\keyword{classes}
