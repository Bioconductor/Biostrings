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
  based on the RNAString alphabet (see below).
  In addition, the letters stored in a \code{RNAString} object
  are encoded in a way that optimizes fast search algorithms.
}

\section{The RNAString alphabet}{
  This alphabet contains all letters from the
  IUPAC Extended Genetic Alphabet (see \code{\link{IUPAC_CODE_MAP}})
  where \code{"T"} is replaced by \code{"U"}
  + the gap letter \code{"-"}.
  It is stored in the \code{RNA_ALPHABET} constant (character vector).
  The \code{alphabet} method also returns
  \code{RNA_ALPHABET} when applied to a \code{RNAString}
  object and is provided for convenience only.
}

\section{Constructor-like functions and generics}{
  In the code snippet below,
  \code{src} can be a character vector
  or a \code{\link{BString}} (or derived) object.

  \describe{
    \item{}{
      \code{RNAString(src)}:
      Tries to convert \code{src} into an \code{RNAString} object.
      \code{src} can be a single string (character vector of length 1),
      a \link{BString} object or a \link{DNAString} object.
    }
  }
}

\section{Other functions and generics}{
  In the code snippets below,
  \code{x} must be a \code{\link{BString}} (or derived) object.

  \describe{
    \item{}{
      \code{alphabet(x)}:
      If \code{x} is of class \code{RNAString}, return the RNAString
      alphabet (see above).
      See the corrsponding man pages when \code{x} is of class
      \code{\link{BString}}, \code{\link{DNAString}} or \code{\link{AAString}}.
    }
  }
}

\author{H. Pages}

\seealso{
  \code{\link{IUPAC_CODE_MAP}},
  \code{\link{letter}},
  \link{BString-class},
  \link{DNAString-class},
  \link{BStringViews-class}
}

\examples{
  d <- DNAString("TTGAAAA-CTC-N")
  r <- RNAString(d)
  alphabet(r)
  RNA_ALPHABET
}

\keyword{methods}
\keyword{classes}
