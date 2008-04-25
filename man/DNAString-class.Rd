\name{DNAString-class}
\docType{class}

% Classes
\alias{class:DNAString}
\alias{DNAString-class}

% Constants:
\alias{DNA_ALPHABET}

% Constructor-like functions and generics:
\alias{DNAString}

% Other functions and generics:
\alias{alphabet,DNAString-method}


\title{DNAString objects}

\description{
  A DNAString object allows efficient storage
  and manipulation of a long DNA sequence.
}

\details{
  The DNAString class is a direct \link{XString} subtype
  (with no additional slot).
  Therefore all functions and methods described in the
  \link{XString} man page also work
  with a DNAString object (inheritance).

  Unlike the \link{BString} container that allows storage
  of any single string (based on a single-byte character set)
  the DNAString container can only store a string
  based on the DNA alphabet (see below).
  In addition, the letters stored in a DNAString object
  are encoded in a way that optimizes fast search algorithms.
}

\section{The DNA alphabet}{
  This alphabet contains all letters from the
  IUPAC Extended Genetic Alphabet (see \code{\link{IUPAC_CODE_MAP}})
  + the gap letter \code{"-"}.
  It is stored in the \code{DNA_ALPHABET} constant (character vector).
  The \code{alphabet} method also returns
  \code{DNA_ALPHABET} when applied to a DNAString
  object and is provided for convenience only.
}

\section{Constructor-like functions and generics}{
  In the code snippet below,
  \code{src} can be a single string (character vector of length 1),
  a \link{BString} object or an \link{RNAString} object.

  \describe{
    \item{}{
      \code{DNAString(src)}:
      Tries to convert \code{src} into a DNAString object by reading
      \code{nchar} letters starting at position \code{start} in \code{src}.
    }
  }
}

\section{Other functions and generics}{
  In the code snippets below,
  \code{x} must be an \link{XString} object.

  \describe{
    \item{}{
      \code{alphabet(x)}:
      If \code{x} is a DNAString object, then return the DNA
      alphabet (see above).
      See the corresponding man pages when \code{x} is of class
      \link{BString}, \link{RNAString} or \link{AAString}.
    }
  }
}

\author{H. Pages}

\seealso{
  \code{\link{IUPAC_CODE_MAP}},
  \code{\link{letter}},
  \code{\link{subXString}},
  \link{XString-class},
  \link{RNAString-class},
  \link{XStringViews-class}
}

\examples{
  DNA_ALPHABET
  d <- DNAString("TTGAAAA-CTC-N")
  length(d)
  alphabet(d)
}

\keyword{methods}
\keyword{classes}
