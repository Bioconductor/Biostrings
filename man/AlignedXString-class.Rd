\name{AlignedXString-class}
\docType{class}

% Classes
\alias{class:AlignedXString}
\alias{AlignedXString-class}
\alias{AlignedXString}

% Constructor-like functions and generics:
\alias{initialize,AlignedXString-method}

% Accesor methods:
\alias{unaligned}
\alias{unaligned,AlignedXString-method}
\alias{aligned}
\alias{aligned,AlignedXString-method}
\alias{start,AlignedXString-method}
\alias{end,AlignedXString-method}
\alias{width,AlignedXString-method}
\alias{indels}
\alias{indels,AlignedXString-method}
\alias{length,AlignedXString-method}
\alias{nchar,AlignedXString-method}
\alias{alphabet,AlignedXString-method}

% Standard generic methods:
\alias{show,AlignedXString-method}
\alias{as.character,AlignedXString-method}
\alias{toString,AlignedXString-method}


\title{AlignedXString objects}

\description{
  The \code{AlignedXString} class is a container for storing
  an alignment.
}

\details{
  Before we define the notion of alignment, we introduce the notion of
  "filled-with-gaps subsequence". A "filled-with-gaps subsequence" of
  a string string1 is obtained by inserting 0 or any number of gaps in
  a subsequence of s1. For example L-A--ND and A--N-D are
  "filled-with-gaps subsequences" of LAND. An alignment between two strings
  string1 and string2 results in two strings (align1 and align2) that have
  the same length and are "filled-with-gaps subsequences" of string1 and string2.

  For example, this is an alignment between LAND and LEAVES:
  \preformatted{
    L-A
    LEA
  }

  An alignment can be seen as a compact representation of one set of basic
  operations that transforms string1 into align1. There are 3 different kinds
  of basic operations: "insertions" (gaps in align1), "deletions" (gaps in
  align2),  "replacements".
  The above alignment represents the following basic operations:
  \preformatted{
    insert E at pos 2
    insert V at pos 4
    insert E at pos 5
    replace by S at pos 6 (N is replaced by S)
    delete at pos 7 (D is deleted)
  }
  Note that "insert X at pos i" means that all letters at a position >= i
  are moved 1 place to the right before X is actually inserted.

  There are many possible alignments between two given strings string1 and
  string2 and a common problem is to find the one (or those ones) with the
  highest score, i.e. with the lower total cost in terms of basic operations.
}

\section{Accesor methods}{
  In the code snippets below,
  \code{x} is a \code{AlignedXString} object.

  \describe{
    \item{}{
      \code{unaligned(x)}:
      The original string.
    }
    \item{}{
      \code{aligned(x)}:
      The "filled-with-gaps subsequence" representing the aligned substring.
    }
    \item{}{
      \code{start(x)}:
      The start of the aligned substring.
    }
    \item{}{
      \code{end(x)}:
      The end of the aligned substring.
    }
    \item{}{
      \code{width(x)}:
      The width of the aligned substring, ignoring gaps.
    }
    \item{}{
      \code{indels(x)}:
      The positions, in the form of an \code{IRanges} object, of the insertions or
      deletions (depending on what the \code{AlignedXString} object represents).
    }
    \item{}{
      \code{length(x)} or \code{nchar(x)}:
      The length of the \code{aligned(x)}.
    }
    \item{}{
      \code{alphabet(x)}:
      Equivalent to \code{alphabet(unaligned(x))}.
    }
    \item{}{
      \code{as.character(x)} and \code{toString(x)}:
      Converts \code{aligned(x)} to a character vector.
    }
  }
}

\author{P. Aboyoun and H. Pages}

\seealso{
  \code{\link{pairwiseAlignment}},
  \code{\link{PairwiseAlignment-class}},
  \code{\link{XString-class}}
}

\examples{
  pattern <- AAString("LAND")
  subject <- AAString("LEAVES")
  nw1 <- pairwiseAlignment(pattern, subject, substitutionMatrix = "BLOSUM50", gapOpening = -3, gapExtension = -1)
  alignedPattern <- pattern(nw1)
  unaligned(alignedPattern)
  aligned(alignedPattern)
  as.character(alignedPattern)
  length(alignedPattern)
}

\keyword{methods}
\keyword{classes}
