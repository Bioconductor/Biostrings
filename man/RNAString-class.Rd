\name{RNAString-class}
\docType{class}

% Classes
\alias{class:RNAString}
\alias{RNAString-class}

% Constants:
\alias{RNA_ALPHABET}

% Constructor-like functions and generics:
\alias{RNAString}


\title{RNAString objects}

\description{
  An RNAString object allows efficient storage
  and manipulation of a long RNA sequence.
}

\details{
  The RNAString class is a direct \link{XString} subclass
  (with no additional slot).
  Therefore all functions and methods described in the
  \link{XString} man page also work
  with an RNAString object (inheritance).

  Unlike the \link{BString} container that allows storage
  of any single string (based on a single-byte character set)
  the RNAString container can only store a string
  based on the RNA alphabet (see below).
  In addition, the letters stored in an RNAString object
  are encoded in a way that optimizes fast search algorithms.
}

\section{The RNA alphabet}{
  This alphabet contains all letters from the
  IUPAC Extended Genetic Alphabet (see \code{?\link{IUPAC_CODE_MAP}})
  where \code{"T"} is replaced by \code{"U"}
  + the gap (\code{"-"}) and the hard masking (\code{"+"}) letters.
  It is stored in the \code{RNA_ALPHABET} constant (character vector).
  The \code{alphabet} method also returns
  \code{RNA_ALPHABET} when applied to an RNAString
  object and is provided for convenience only.
}

\section{Constructor-like functions and generics}{
  In the code snippet below,
  \code{x} can be a single string (character vector of length 1),
  a \link{BString} object or a \link{DNAString} object.

  \describe{
    \item{}{
      \code{RNAString(x, start=1, nchar=NA, check=TRUE)}:
      Tries to convert \code{x} into an RNAString object by reading
      \code{nchar} letters starting at position \code{start} in \code{x}.
    }
  }
}

\section{Accessor methods}{
  In the code snippet below, \code{x} is an \link{RNAString} object.

  \describe{
    \item{}{
      \code{alphabet(x)}:
      If \code{x} is an RNAString object, then return the RNA
      alphabet (see above).
      See the corresponding man pages when \code{x} is a
      \link{BString}, \link{DNAString} or \link{AAString} object.
    }
  }
}

\author{H. Pages}

\seealso{
  \code{\link{IUPAC_CODE_MAP}},
  \code{\link{letter}},
  \link{XString-class},
  \link{DNAString-class},
  \code{\link{reverseComplement}},
  \code{\link{alphabetFrequency}}
}

\examples{
  d <- DNAString("TTGAAAA-CTC-N")
  r <- RNAString(d)
  r
  alphabet(r)
  RNA_ALPHABET

  ## When comparing an RNAString object with a DNAString object,
  ## U and T are considered equals:
  r == d
}

\keyword{methods}
\keyword{classes}
