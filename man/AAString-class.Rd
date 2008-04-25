\name{AAString-class}
\docType{class}

% Classes
\alias{class:AAString}
\alias{AAString-class}

% Constants:
\alias{AA_ALPHABET}

% Constructor-like functions and generics:
\alias{AAString}

% Other functions and generics:
\alias{alphabet,AAString-method}


\title{AAString objects}

\description{
  An AAString object allows efficient storage
  and manipulation of a long amino acid sequence.
}

\details{
  The AAString class is a direct \link{XString} subtype
  (with no additional slot).
  Therefore all functions and methods described in the
  \link{XString} man page also work
  with an AAString object (inheritance).

  Unlike the \link{BString} container that allows storage
  of any single string (based on a single-byte character set)
  the AAString container can only store a string
  based on the Amino Acid alphabet (see below).
}

\section{The Amino Acid alphabet}{
  This alphabet contains all letters from the
  Single-Letter Amino Acid Code (see \code{\link{AMINO_ACID_CODE}})
  + the stop letter \code{"*"}.
  It is stored in the \code{AA_ALPHABET} constant (character vector).
  The \code{alphabet} method also returns
  \code{AA_ALPHABET} when applied to an AAString
  object and is provided for convenience only.
}

\section{Constructor-like functions and generics}{
  In the code snippet below,
  \code{src} can be a single string (character vector of length 1)
  or a \link{BString} object.

  \describe{
    \item{}{
      \code{AAString(src)}:
      Tries to convert \code{src} into an AAString object by reading
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
      If \code{x} is an AAString object, then return the Amino Acid
      alphabet (see above).
      See the corresponding man pages when \code{x} is of class
      \link{BString}, \link{DNAString} or \link{RNAString}.
    }
  }
}

\author{H. Pages}

\seealso{
  \code{\link{AMINO_ACID_CODE}},
  \code{\link{subXString}},
  \code{\link{letter}},
  \link{XString-class},
  \link{XStringViews-class}
}

\examples{
  AA_ALPHABET
  a <- AAString("MARKSLEMSIR*")
  length(a)
  alphabet(a)
}

\keyword{methods}
\keyword{classes}
