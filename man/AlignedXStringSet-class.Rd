\name{AlignedXStringSet-class}
\docType{class}

% Classes
\alias{class:AlignedXStringSet0}
\alias{AlignedXStringSet0-class}
\alias{AlignedXStringSet0}
\alias{class:AlignedXStringSet}
\alias{AlignedXStringSet-class}
\alias{AlignedXStringSet}
\alias{class:QualityAlignedXStringSet}
\alias{QualityAlignedXStringSet-class}
\alias{QualityAlignedXStringSet}

% Accessor methods:
\alias{unaligned}
\alias{unaligned,AlignedXStringSet0-method}
\alias{aligned}
\alias{aligned,AlignedXStringSet0-method}
\alias{start,AlignedXStringSet0-method}
\alias{end,AlignedXStringSet0-method}
\alias{width,AlignedXStringSet0-method}
\alias{indel}
\alias{indel,AlignedXStringSet0-method}
\alias{nindel}
\alias{nindel,AlignedXStringSet0-method}
\alias{nchar,AlignedXStringSet0-method}
\alias{seqtype,AlignedXStringSet0-method}
\alias{ranges,AlignedXStringSet0-method}

% Standard generic methods:
\alias{show,AlignedXStringSet0-method}
\alias{as.character,AlignedXStringSet0-method}
\alias{toString,AlignedXStringSet0-method}

% Internal:
\alias{parallel_slot_names,AlignedXStringSet0-method}
\alias{parallelVectorNames,AlignedXStringSet0-method}

\title{AlignedXStringSet and QualityAlignedXStringSet objects}

\description{
  The \code{AlignedXStringSet} and \code{QualityAlignedXStringSet} classes are
  containers for storing an aligned \code{XStringSet}.
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

\section{Accessor methods}{
  In the code snippets below,
  \code{x} is a \code{AlignedXStringSet} or \code{QualityAlignedXStringSet} object.

  \describe{
    \item{\code{unaligned(x)}:}{
      The original string.
    }
    \item{\code{aligned(x, degap = FALSE)}:}{
      If \code{degap = FALSE}, the "filled-with-gaps subsequence" representing
      the aligned substring. If \code{degap = TRUE}, the "gap-less subsequence"
      representing the aligned substring.
    }
    \item{\code{ranges(x)}:}{ The bounds of the aligned substring.
    }
    \item{\code{start(x)}:}{
      The start of the aligned substring.
    }
    \item{\code{end(x)}:}{
      The end of the aligned substring.
    }
    \item{\code{width(x)}:}{
      The width of the aligned substring, ignoring gaps.
    }
    \item{\code{indel(x)}:}{
      The positions, in the form of an \code{IRanges} object, of the insertions or
      deletions (depending on what \code{x} represents).
    }
    \item{\code{nindel(x)}:}{
      A two-column matrix containing the length and sum of the widths for each of
      the elements returned by \code{indel}.
    }
    \item{\code{length(x)}:}{
      The length of the \code{aligned(x)}.
    }
    \item{\code{nchar(x)}:}{
      The nchar of the \code{aligned(x)}.
    }
    \item{\code{alphabet(x)}:}{
      Equivalent to \code{alphabet(unaligned(x))}.
    }
    \item{\code{as.character(x)}:}{
      Converts \code{aligned(x)} to a character vector.
    }
    \item{\code{toString(x)}:}{
      Equivalent to \code{toString(as.character(x))}.
    }
  }
}

\section{Subsetting methods}{
  \describe{
    \item{\code{x[i]}:}{
      Returns a new \code{AlignedXStringSet} or \code{QualityAlignedXStringSet}
      object made of the selected elements.
    }
    \item{\code{rep(x, times)}:}{
      Returns a new \code{AlignedXStringSet} or \code{QualityAlignedXStringSet}
      object made of the repeated elements.
    }
  }
}

\author{P. Aboyoun}

\seealso{
  \code{\link{pairwiseAlignment}},
  \code{\link{PairwiseAlignments-class}},
  \code{\link{XStringSet-class}}
}

\examples{
pattern <- AAString("LAND")
subject <- AAString("LEAVES")
pa1 <- pairwiseAlignment(pattern, subject, substitutionMatrix="BLOSUM50",
                         gapOpening=3, gapExtension=1)

alignedPattern <- pattern(pa1)
class(alignedPattern)  # AlignedXStringSet object

unaligned(alignedPattern)
aligned(alignedPattern)
as.character(alignedPattern)
nchar(alignedPattern)
}

\keyword{methods}
\keyword{classes}
