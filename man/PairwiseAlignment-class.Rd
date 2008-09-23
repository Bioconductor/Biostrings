\name{PairwiseAlignment-class}
\docType{class}

% Classes
\alias{class:PairwiseAlignment}
\alias{PairwiseAlignment-class}
\alias{PairwiseAlignment}
\alias{class:PairwiseAlignmentSummary}
\alias{PairwiseAlignmentSummary-class}
\alias{PairwiseAlignmentSummary}

% Constructor-like functions and generics:
\alias{initialize,PairwiseAlignment-method}

% Accesor methods:
\alias{pattern,PairwiseAlignment-method}
\alias{subject,PairwiseAlignment-method}
\alias{type}
\alias{type,PairwiseAlignment-method}
\alias{type,PairwiseAlignmentSummary-method}
\alias{score}
\alias{score,PairwiseAlignment-method}
\alias{score,PairwiseAlignmentSummary-method}
\alias{nindel,PairwiseAlignment-method}
\alias{nindel,PairwiseAlignmentSummary-method}
\alias{length,PairwiseAlignment-method}
\alias{length,PairwiseAlignmentSummary-method}
\alias{nchar,PairwiseAlignment-method}
\alias{nchar,PairwiseAlignmentSummary-method}
\alias{alphabet,PairwiseAlignment-method}
\alias{Views,PairwiseAlignment-method}

% Standard generic methods:
\alias{show,PairwiseAlignment-method}
\alias{show,PairwiseAlignmentSummary-method}
\alias{summary,PairwiseAlignment-method}
\alias{as.character,PairwiseAlignment-method}
\alias{toString,PairwiseAlignment-method}
\alias{as.matrix,PairwiseAlignment-method}
\alias{[,PairwiseAlignment-method}
\alias{[<-,PairwiseAlignment-method}
\alias{rep,PairwiseAlignment-method}


\title{PairwiseAlignment and PairwiseAlignmentSummary objects}

\description{
  The \code{PairwiseAlignment} class is a container for storing an alignment.
  The \code{PairwiseAlignmentSummary} class is a container for storing the
  summary of an alignment.
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
  \code{x} is a \code{PairwiseAlignment} object, except otherwise noted.

  \describe{
    \item{}{
      \code{pattern(x)}:
      The \code{AlignedXStringSet} object for the pattern.
    }
    \item{}{
      \code{subject(x)}:
      The \code{AlignedXStringSet} object for the subject.
    }
    \item{}{
      \code{type(x)}:
      The type of the alignment (\code{"global"}, \code{"local"},
      \code{"overlap"}, \code{"patternOverlap"}, or \code{"subjectOverlap"}).
      There is a method for \code{PairwiseAlignmentSummary} as well.
    }
    \item{}{
      \code{score(x)}:
      The score of the alignment (integer).
      There is a method for \code{PairwiseAlignmentSummary} as well.
    }
    \item{}{
      \code{nindel(x)}:
      An \code{InDel} object containing the number of insertions and deletions.
    }
    \item{}{
      \code{length(x)}:
      The length of the \code{aligned(pattern(x))} and \code{aligned(subject(x))}.
      There is a method for \code{PairwiseAlignmentSummary} as well.
    }
    \item{}{
      \code{nchar(x)}:
      The nchar of the \code{aligned(pattern(x))} and \code{aligned(subject(x))}.
      There is a method for \code{PairwiseAlignmentSummary} as well.
    }
    \item{}{
      \code{alphabet(x)}:
      Equivalent to \code{alphabet(unaligned(subject(x)))}.
    }
    \item{}{
      \code{summary(object, ...)}:
      Generates a summary for the \code{PairwiseAlignment}.
    }
    \item{}{
      \code{as.character(x)} and \code{toString(x)}:
      Converts \code{x} to a named two row matrix containing the results of
      \code{as.character(pattern(x))} and \code{as.character(subject(x))}.
    }
    \item{}{
      \code{as.matrix(x)}:
      Returns a character matrix with \code{nrow = length(x)} and
      \code{ncol = nchar(unaligned(subject(x)))} containing the "exploded"
      representation of the aligned patterns. This operation ignores insertions.
    }
    \item{}{
      \code{Views(subject, start=NA, end=NA, names=NULL)}:
      The \code{XStringViews} object that represents the pairwise alignments along \code{unaligned(subject(subject))}.
      The \code{start} and \code{end} arguments must be either \code{NA} or an integer vector of length 1 that
      denotes the offset from \code{start(subject(subject))}.
    }
  }
}

\section{Subsetting methods}{
  \describe{
    \item{}{
      \code{x[i]}:
      Returns a new \code{PairwiseAlignment} object made of the selected elements.
    }
    \item{}{
      \code{rep(x, times)}:
      Returns a new \code{PairwiseAlignment} object made of the repeated elements.
    }
  }
}

\author{P. Aboyoun and H. Pages}

\seealso{
  \code{\link{pairwiseAlignment}},
  \code{\link{AlignedXStringSet-class}},
  \code{\link{XString-class}},
  \code{\link{XStringViews-class}},
  \code{\link{match-utils}}
}

\examples{
  pattern <- AAStringSet(c("HLDNLKGTF", "HVDDMPNAL"))
  subject <- AAString("SMDDTEKMSMKL")
  nw1 <- pairwiseAlignment(pattern, subject, substitutionMatrix = "BLOSUM50", gapOpening = -3, gapExtension = -1)
  pattern(nw1)
  subject(nw1)
  as.character(nw1)
  as.matrix(nw1)
  nchar(nw1)
  score(nw1)
  nw1
}

\keyword{methods}
\keyword{classes}
