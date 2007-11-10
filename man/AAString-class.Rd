\name{AAString-class}
\docType{class}

% Classes
\alias{class:AAString}
\alias{AAString-class}

% Constants:
\alias{AA_ALPHABET}

% Constructor-like functions and generics:
\alias{initialize,AAString-method}
\alias{AAString}

% Other functions and generics:
\alias{alphabet,AAString-method}


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
  based on the AAString alphabet (see below).
}

\section{The AAString alphabet}{
  This alphabet contains all letters from the
  Single-Letter Amino Acid Code (see \code{\link{AMINO_ACID_CODE}})
  + the stop letter \code{"*"}.
  It is stored in the \code{AA_ALPHABET} constant (character vector).
  The \code{alphabet} method also returns
  \code{AA_ALPHABET} when applied to a \code{AAString}
  object and is provided for convenience only.
}

\section{Constructor-like functions and generics}{
  In the code snippet below,
  \code{src} can be a character vector
  or a \code{\link{BString}} (or derived) object.

  \describe{
    \item{}{
      \code{AAString(src)}:
      Tries to convert \code{src} into an \code{AAString} object.
      \code{src} can be a single string (character vector of length 1)
      or a \link{BString} object.
    }
  }
}

\section{Other functions and generics}{
  In the code snippets below,
  \code{x} must be a \code{\link{BString}} (or derived) object.

  \describe{
    \item{}{
      \code{alphabet(x)}:
      If \code{x} is of class \code{AAString}, return the AAString
      alphabet (see above).
      See the corresponding man pages when \code{x} is of class
      \code{\link{BString}}, \code{\link{DNAString}} or \code{\link{RNAString}}.
    }
  }
}

\author{H. Pages}

\seealso{
  \code{\link{AMINO_ACID_CODE}},
  \code{\link{subBString}},
  \code{\link{letter}},
  \link{BString-class},
  \link{BStringViews-class}
}

\examples{
  AA_ALPHABET
  a <- AAString("MARKSLEMSIR*")
  length(a)
  alphabet(a)
}

\keyword{methods}
\keyword{classes}
