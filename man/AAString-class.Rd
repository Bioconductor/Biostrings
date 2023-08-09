\name{AAString-class}
\docType{class}

% Classes
\alias{class:AAString}
\alias{AAString-class}

% Constants:
\alias{AA_ALPHABET}
\alias{AA_STANDARD}
\alias{AA_PROTEINOGENIC}

% Constructor-like functions and generics:
\alias{AAString}


\title{AAString objects}

\description{
  An AAString object allows efficient storage
  and manipulation of a long amino acid sequence.
}

\usage{
AAString(x="", start=1, nchar=NA)

## Predefined constants:
AA_ALPHABET       # full Amino Acid alphabet
AA_STANDARD       # first 20 letters only
AA_PROTEINOGENIC  # first 22 letters only
}

\arguments{
  \item{x}{
    A single string.
  }
  \item{start, nchar}{
    Where to start reading from in \code{x} and how many letters to read.
  }
}

\details{
  The AAString class is a direct \link{XString} subclass
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
  Single-Letter Amino Acid Code (see \code{?\link{AMINO_ACID_CODE}})
  plus \code{"*"} (the \emph{stop} letter), \code{"-"} (the \emph{gap}
  letter), \code{"+"} (the \emph{hard masking} letter), and \code{"."}
  (the \emph{not a letter} or \emph{not available} letter).
  It is stored in the \code{AA_ALPHABET} predefined constant (character
  vector).

  The \code{alphabet()} function returns \code{AA_ALPHABET} when
  applied to an AAString object.
}

\section{Constructor-like functions and generics}{
  In the code snippet below,
  \code{x} can be a single string (character vector of length 1)
  or a \link{BString} object.

  \describe{
    \item{\code{AAString(x="", start=1, nchar=NA)}:}{
      Tries to convert \code{x} into an AAString object by reading
      \code{nchar} letters starting at position \code{start} in \code{x}.
    }
  }
}

\section{Accessor methods}{
  In the code snippet below, \code{x} is an \link{AAString} object.

  \describe{
    \item{\code{alphabet(x)}:}{
      If \code{x} is an AAString object, then return the Amino Acid
      alphabet (see above).
      See the corresponding man pages when \code{x} is a
      \link{BString}, \link{DNAString} or \link{RNAString} object.
    }
  }
}

\section{Display}{
  The letters in an AAString object are colored when displayed by the
  \code{show()} method. Set global option \code{Biostrings.coloring}
  to FALSE to turn off this coloring.
}

\author{H. Pagès}

\seealso{
  \code{\link{AMINO_ACID_CODE}},
  \code{\link{letter}},
  \link{XString-class},
  \code{\link{alphabetFrequency}}
}

\examples{
AA_ALPHABET
a <- AAString("MARKSLEMSIR*")
length(a)
alphabet(a)
}

\keyword{methods}
\keyword{classes}
