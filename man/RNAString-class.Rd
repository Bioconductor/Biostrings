\name{RNAString-class}
\docType{class}

% Classes
\alias{class:RNAString}
\alias{RNAString-class}

% Constants:
\alias{RNA_BASES}
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
  This alphabet is the same as the DNA alphabet, except that \code{"T"}
  is replaced by \code{"U"}. See \code{?\link{DNA_ALPHABET}} for more
  information about the DNA alphabet.
  The RNA alphabet is stored in the \code{RNA_ALPHABET} predefined constant
  (character vector).

  The \code{alphabet()} function returns \code{RNA_ALPHABET} when
  applied to an RNAString object.
}

\section{Constructor-like functions and generics}{
  In the code snippet below,
  \code{x} can be a single string (character vector of length 1),
  a \link{BString} object or a \link{DNAString} object.

  \describe{
    \item{\code{RNAString(x="", start=1, nchar=NA)}:}{
      Tries to convert \code{x} into an RNAString object by reading
      \code{nchar} letters starting at position \code{start} in \code{x}.
    }
  }
}

\section{Accessor methods}{
  In the code snippet below, \code{x} is an RNAString object.

  \describe{
    \item{\code{alphabet(x, baseOnly=FALSE)}:}{
      If \code{x} is an RNAString object, then return the RNA
      alphabet (see above).
      See the corresponding man pages when \code{x} is a
      \link{BString}, \link{DNAString} or \link{AAString} object.
    }
  }
}

\section{Display}{
  The letters in an RNAString object are colored when displayed by the
  \code{show()} method. Set global option \code{Biostrings.coloring}
  to FALSE to turn off this coloring.
}

\author{H. Pagès}

\seealso{
  \itemize{
    \item The \link{RNAStringSet} class to represent a collection of
          RNAString objects.

    \item The \link{XString} and \link{DNAString} classes.

    \item \code{\link{reverseComplement}}

    \item \code{\link{alphabetFrequency}}

    \item \code{\link{IUPAC_CODE_MAP}}

    \item \code{\link{letter}}
  }
}

\examples{
RNA_BASES
RNA_ALPHABET
dna <- DNAString("TTGAAAA-CTC-N")
rna <- RNAString(dna)
rna  # 'options(Biostrings.coloring=FALSE)' to turn off coloring

alphabet(rna)                 # RNA_ALPHABET
alphabet(rna, baseOnly=TRUE)  # RNA_BASES

## When comparing an RNAString object with a DNAString object,
## U and T are considered equals:
rna == dna  # TRUE
}

\keyword{methods}
\keyword{classes}
