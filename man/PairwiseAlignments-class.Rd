\name{PairwiseAlignments-class}
\docType{class}

% Classes
\alias{class:PairwiseAlignments}
\alias{PairwiseAlignments-class}
\alias{PairwiseAlignments}
\alias{class:PairwiseAlignmentsSingleSubject}
\alias{PairwiseAlignmentsSingleSubject-class}
\alias{PairwiseAlignmentsSingleSubject}
\alias{class:PairwiseAlignmentsSingleSubjectSummary}
\alias{PairwiseAlignmentsSingleSubjectSummary-class}
\alias{PairwiseAlignmentsSingleSubjectSummary}

% Constructor-like functions and generics:
\alias{PairwiseAlignments,XString,XString-method}
\alias{PairwiseAlignments,XStringSet,missing-method}
\alias{PairwiseAlignments,character,missing-method}
\alias{PairwiseAlignments,character,character-method}
\alias{PairwiseAlignmentsSingleSubject,XString,XString-method}
\alias{PairwiseAlignmentsSingleSubject,XStringSet,missing-method}
\alias{PairwiseAlignmentsSingleSubject,character,missing-method}
\alias{PairwiseAlignmentsSingleSubject,character,character-method}

% Accessor methods:
\alias{pattern,PairwiseAlignments-method}
\alias{subject,PairwiseAlignments-method}
\alias{type}
\alias{type,PairwiseAlignments-method}
\alias{score,PairwiseAlignments-method}
\alias{insertion,PairwiseAlignments-method}
\alias{deletion,PairwiseAlignments-method}
\alias{indel,PairwiseAlignments-method}
\alias{nindel,PairwiseAlignments-method}
\alias{length,PairwiseAlignments-method}
\alias{nchar,PairwiseAlignments-method}
\alias{seqtype,PairwiseAlignments-method}

% Standard generic methods:
\alias{show,PairwiseAlignments-method}
\alias{[,PairwiseAlignments-method}
\alias{[<-,PairwiseAlignments-method}
\alias{rep,PairwiseAlignments-method}

% Methods for PairwiseAlignmentsSingleSubject:
\alias{summary,PairwiseAlignmentsSingleSubject-method}
\alias{Views,PairwiseAlignmentsSingleSubject-method}
\alias{aligned,PairwiseAlignmentsSingleSubject-method}
\alias{as.character,PairwiseAlignmentsSingleSubject-method}
\alias{toString,PairwiseAlignmentsSingleSubject-method}
\alias{as.matrix,PairwiseAlignmentsSingleSubject-method}

% Methods for PairwiseAlignmentsSingleSubjectSummary:
\alias{type,PairwiseAlignmentsSingleSubjectSummary-method}
\alias{score,PairwiseAlignmentsSingleSubjectSummary-method}
\alias{nindel,PairwiseAlignmentsSingleSubjectSummary-method}
\alias{length,PairwiseAlignmentsSingleSubjectSummary-method}
\alias{nchar,PairwiseAlignmentsSingleSubjectSummary-method}
\alias{show,PairwiseAlignmentsSingleSubjectSummary-method}


\title{PairwiseAlignments, PairwiseAlignmentsSingleSubject, and
  PairwiseAlignmentsSingleSubjectSummary objects}

\description{
  The \code{PairwiseAlignments} class is a container for storing
  a set of pairwise alignments.

  The \code{PairwiseAlignmentsSingleSubject} class is a container for storing
  a set of pairwise alignments with a single subject.

  The \code{PairwiseAlignmentsSingleSubjectSummary} class is a container for storing
  the summary of a set of pairwise alignments.
}

\usage{
## Constructors:
## When subject is missing, pattern must be of length 2
\S4method{PairwiseAlignments}{XString,XString}(pattern, subject,
  type = "global", substitutionMatrix = NULL, gapOpening = 0, gapExtension = 1)
\S4method{PairwiseAlignments}{XStringSet,missing}(pattern, subject,
  type = "global", substitutionMatrix = NULL, gapOpening = 0, gapExtension = 1)
\S4method{PairwiseAlignments}{character,character}(pattern, subject,
  type = "global", substitutionMatrix = NULL, gapOpening = 0, gapExtension = 1,
  baseClass = "BString")
\S4method{PairwiseAlignments}{character,missing}(pattern, subject,
  type = "global", substitutionMatrix = NULL, gapOpening = 0, gapExtension = 1,
  baseClass = "BString")
}

\arguments{
  \item{pattern}{a character vector of length 1 or 2, an \code{\link{XString}},
    or an \code{\link{XStringSet}} object of length 1 or 2.}
  \item{subject}{a character vector of length 1 or an \code{\link{XString}} object.}
  \item{type}{type of alignment. One of \code{"global"}, \code{"local"},
    \code{"overlap"}, \code{"global-local"}, and \code{"local-global"} where
      \code{"global"} = align whole strings with end gap penalties,
      \code{"local"} = align string fragments,
      \code{"overlap"} = align whole strings without end gap penalties,
      \code{"global-local"} = align whole strings in \code{pattern} with
        consecutive subsequence of \code{subject},
      \code{"local-global"} = align consecutive subsequence of \code{pattern}
        with whole strings in \code{subject}.}
  \item{substitutionMatrix}{substitution matrix for the alignment. If NULL,
     the diagonal values and off-diagonal values are set to 0 and 1 respectively.}
  \item{gapOpening}{the cost for opening a gap in the alignment.}
  \item{gapExtension}{the incremental cost incurred along the length of the gap
     in the alignment.}
  \item{baseClass}{the base \code{\link{XString}} class to use in the alignment.}
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

\section{Object extraction methods}{
  In the code snippets below,
  \code{x} is a \code{PairwiseAlignments} object, except otherwise noted.

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
      \code{summary(object, ...)}:
      Generates a summary for the \code{PairwiseAlignments}.
    }
  }
}
\section{General information methods}{
  In the code snippets below,
  \code{x} is a \code{PairwiseAlignments} object, except otherwise noted.

  \describe{
    \item{}{
      \code{alphabet(x)}:
      Equivalent to \code{alphabet(unaligned(subject(x)))}.
    }
    \item{}{
      \code{length(x)}:
      The length of the \code{aligned(pattern(x))} and
    \code{aligned(subject(x))}. 
      There is a method for \code{PairwiseAlignmentsSingleSubjectSummary} as well.
    }
    \item{}{
      \code{type(x)}:
      The type of the alignment (\code{"global"}, \code{"local"},
      \code{"overlap"}, \code{"global-local"}, or \code{"local-global"}).
      There is a method for \code{PairwiseAlignmentsSingleSubjectSummary} as well.
    }
  }
}
\section{Aligned sequence methods}{
  In the code snippets below,
  \code{x} is a \code{PairwiseAlignmentsSingleSubject} object, except
    otherwise noted. 

  \describe{
    \item{}{
      \code{aligned(x, degap = FALSE, gapCode="-", endgapCode="-")}:
      If \code{degap = FALSE}, "align" the alignments by returning an
      \code{XStringSet} object containing the aligned patterns without
      insertions. If \code{degap = TRUE}, returns
      \code{aligned(pattern(x), degap=TRUE)}.
      The \code{gapCode} and \code{endgapCode} arguments denote the code in the
      appropriate \code{\link{alphabet}} to use for the internal and end gaps.
    }
    \item{}{
      \code{as.character(x)}:
      Converts \code{aligned(x)} to a character vector.
    }
    \item{}{
      \code{as.matrix(x)}:
      Returns an "exploded" character matrix representation of \code{aligned(x)}.
    }
    \item{}{
      \code{toString(x)}:
      Equivalent to \code{toString(as.character(x))}.
    }
  }
}
\section{Subject position methods}{
  In the code snippets below,
  \code{x} is a \code{PairwiseAlignmentsSingleSubject} object, except otherwise noted.

  \describe{
    \item{}{
      \code{consensusMatrix(x, as.prob=FALSE, baseOnly=FALSE, gapCode="-",
        endgapCode="-")}
      See `\link{consensusMatrix}` for more information.
    }
    \item{}{
      \code{consensusString(x)}
      See `\link{consensusString}` for more information.
    }
    \item{}{
      \code{coverage(x, shift=0L, width=NULL, weight=1L)}
      See `\link{coverage,PairwiseAlignmentsSingleSubject-method}` for more
      information.
    }
    \item{}{
      \code{Views(subject, start=NULL, end=NULL, width=NULL, names=NULL)}:
      The \code{XStringViews} object that represents the pairwise alignments
      along \code{unaligned(subject(subject))}. The \code{start} and \code{end}
      arguments must be either \code{NULL}/\code{NA} or an integer vector
      of length 1 that denotes the offset from \code{start(subject(subject))}.
    }
  }
}
\section{Numeric summary methods}{
  In the code snippets below,
  \code{x} is a \code{PairwiseAlignments} object, except otherwise noted.

  \describe{
    \item{}{
      \code{nchar(x)}:
      The nchar of the \code{aligned(pattern(x))} and \code{aligned(subject(x))}.
      There is a method for \code{PairwiseAlignmentsSingleSubjectSummary} as well.
    }
    \item{}{
      \code{insertion(x)}:
      An \code{\link[IRanges:IRangesList-class]{CompressedIRangesList}}
    object containing the 
      locations of the insertions from the perspective of the \code{pattern}.
    }
    \item{}{
      \code{deletion(x)}:
      An \code{\link[IRanges:IRangesList-class]{CompressedIRangesList}}
    object containing the 
      locations of the deletions from the perspective of the \code{pattern}.
    }
    \item{}{
      \code{indel(x)}:
      An \code{InDel} object containing the locations of the insertions and
      deletions from the perspective of the \code{pattern}.
    }
    \item{}{
      \code{nindel(x)}:
      An \code{InDel} object containing the number of insertions and deletions.
    }
    \item{}{
      \code{score(x)}:
      The score of the alignment.
      There is a method for \code{PairwiseAlignmentsSingleSubjectSummary}
      as well.
    }
  }
}

\section{Subsetting methods}{
  \describe{
    \item{}{
      \code{x[i]}:
      Returns a new \code{PairwiseAlignments} object made of the selected
      elements.
    }
    \item{}{
      \code{rep(x, times)}:
      Returns a new \code{PairwiseAlignments} object made of the repeated
      elements.
    }
  }
}

\author{P. Aboyoun}

\seealso{
  \code{\link{pairwiseAlignment}},
  \code{\link{writePairwiseAlignments}},
  \link{AlignedXStringSet-class},
  \link{XString-class},
  \link{XStringViews-class},
  \link{align-utils},
  \code{\link{pid}}
}

\examples{
  PairwiseAlignments("-PA--W-HEAE", "HEAGAWGHE-E")
  pattern <- AAStringSet(c("HLDNLKGTF", "HVDDMPNAL"))
  subject <- AAString("SMDDTEKMSMKL")
  nw1 <- pairwiseAlignment(pattern, subject, substitutionMatrix = "BLOSUM50",
    gapOpening = 3, gapExtension = 1)
  pattern(nw1)
  subject(nw1)
  aligned(nw1)
  as.character(nw1)
  as.matrix(nw1)
  nchar(nw1)
  score(nw1)
  nw1
}

\keyword{methods}
\keyword{classes}
