\name{PairwiseAlignedXStringSet-class}
\docType{class}

% Classes
\alias{class:PairwiseAlignedXStringSet}
\alias{PairwiseAlignedXStringSet-class}
\alias{PairwiseAlignedXStringSet}
\alias{class:PairwiseAlignedFixedSubject}
\alias{PairwiseAlignedFixedSubject-class}
\alias{PairwiseAlignedFixedSubject}
\alias{class:PairwiseAlignedFixedSubjectSummary}
\alias{PairwiseAlignedFixedSubjectSummary-class}
\alias{PairwiseAlignedFixedSubjectSummary}

% Constructor-like functions and generics:
\alias{PairwiseAlignedXStringSet,XString,XString-method}
\alias{PairwiseAlignedXStringSet,XStringSet,missing-method}
\alias{PairwiseAlignedXStringSet,character,missing-method}
\alias{PairwiseAlignedXStringSet,character,character-method}
\alias{PairwiseAlignedFixedSubject,XString,XString-method}
\alias{PairwiseAlignedFixedSubject,XStringSet,missing-method}
\alias{PairwiseAlignedFixedSubject,character,missing-method}
\alias{PairwiseAlignedFixedSubject,character,character-method}

% Accessor methods:
\alias{pattern,PairwiseAlignedXStringSet-method}
\alias{subject,PairwiseAlignedXStringSet-method}
\alias{type}
\alias{type,PairwiseAlignedXStringSet-method}
\alias{score,PairwiseAlignedXStringSet-method}
\alias{insertion,PairwiseAlignedXStringSet-method}
\alias{deletion,PairwiseAlignedXStringSet-method}
\alias{indel,PairwiseAlignedXStringSet-method}
\alias{nindel,PairwiseAlignedXStringSet-method}
\alias{length,PairwiseAlignedXStringSet-method}
\alias{nchar,PairwiseAlignedXStringSet-method}
\alias{xsbasetype,PairwiseAlignedXStringSet-method}

% Standard generic methods:
\alias{show,PairwiseAlignedXStringSet-method}
\alias{[,PairwiseAlignedXStringSet-method}
\alias{[<-,PairwiseAlignedXStringSet-method}
\alias{rep,PairwiseAlignedXStringSet-method}

% Methods for PairwiseAlignedFixedSubject:
\alias{summary,PairwiseAlignedFixedSubject-method}
\alias{Views,PairwiseAlignedFixedSubject-method}
\alias{aligned,PairwiseAlignedFixedSubject-method}
\alias{as.character,PairwiseAlignedFixedSubject-method}
\alias{toString,PairwiseAlignedFixedSubject-method}
\alias{as.matrix,PairwiseAlignedFixedSubject-method}

% Methods for PairwiseAlignedFixedSubjectSummary:
\alias{type,PairwiseAlignedFixedSubjectSummary-method}
\alias{score,PairwiseAlignedFixedSubjectSummary-method}
\alias{nindel,PairwiseAlignedFixedSubjectSummary-method}
\alias{length,PairwiseAlignedFixedSubjectSummary-method}
\alias{nchar,PairwiseAlignedFixedSubjectSummary-method}
\alias{show,PairwiseAlignedFixedSubjectSummary-method}


\title{PairwiseAlignedXStringSet, PairwiseAlignedFixedSubject, and
  PairwiseAlignedFixedSubjectSummary objects}

\description{
  The \code{PairwiseAlignedXStringSet} class is a container for storing an
  elementwise pairwise alignment.
  The \code{PairwiseAlignedFixedSubject} class is a container for storing a
  pairwise alignment with a single subject.
  The \code{PairwiseAlignedFixedSubjectSummary} class is a container for storing
  the summary of an alignment.
}

\usage{
  ## Constructors:
  ## When subject is missing, pattern must be of length 2
  \S4method{PairwiseAlignedXStringSet}{XString,XString}(pattern, subject,
    type = "global", substitutionMatrix = NULL, gapOpening = 0, gapExtension = -1)
  \S4method{PairwiseAlignedXStringSet}{XStringSet,missing}(pattern, subject,
    type = "global", substitutionMatrix = NULL, gapOpening = 0, gapExtension = -1)
  \S4method{PairwiseAlignedXStringSet}{character,character}(pattern, subject,
    type = "global", substitutionMatrix = NULL, gapOpening = 0, gapExtension = -1,
    baseClass = "BString")
  \S4method{PairwiseAlignedXStringSet}{character,missing}(pattern, subject,
    type = "global", substitutionMatrix = NULL, gapOpening = 0, gapExtension = -1,
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
  \code{x} is a \code{PairwiseAlignedXStringSet} object, except otherwise noted.

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
      Generates a summary for the \code{PairwiseAlignedXStringSet}.
    }
  }
}
\section{General information methods}{
  In the code snippets below,
  \code{x} is a \code{PairwiseAlignedXStringSet} object, except otherwise noted.

  \describe{
    \item{}{
      \code{alphabet(x)}:
      Equivalent to \code{alphabet(unaligned(subject(x)))}.
    }
    \item{}{
      \code{length(x)}:
      The length of the \code{aligned(pattern(x))} and
    \code{aligned(subject(x))}. 
      There is a method for \code{PairwiseAlignedFixedSubjectSummary} as well.
    }
    \item{}{
      \code{type(x)}:
      The type of the alignment (\code{"global"}, \code{"local"},
      \code{"overlap"}, \code{"global-local"}, or \code{"local-global"}).
      There is a method for \code{PairwiseAlignedFixedSubjectSummary} as well.
    }
  }
}
\section{Aligned sequence methods}{
  In the code snippets below,
  \code{x} is a \code{PairwiseAlignedFixedSubject} object, except
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
  \code{x} is a \code{PairwiseAlignedFixedSubject} object, except otherwise noted.

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
      See `\link{coverage,PairwiseAlignedFixedSubject-method}` for more
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
  \code{x} is a \code{PairwiseAlignedXStringSet} object, except otherwise noted.

  \describe{
    \item{}{
      \code{nchar(x)}:
      The nchar of the \code{aligned(pattern(x))} and \code{aligned(subject(x))}.
      There is a method for \code{PairwiseAlignedFixedSubjectSummary} as well.
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
      There is a method for \code{PairwiseAlignedFixedSubjectSummary} as well.
    }
  }
}

\section{Subsetting methods}{
  \describe{
    \item{}{
      \code{x[i]}:
      Returns a new \code{PairwiseAlignedXStringSet} object made of the selected
      elements.
    }
    \item{}{
      \code{rep(x, times)}:
      Returns a new \code{PairwiseAlignedXStringSet} object made of the repeated
      elements.
    }
  }
}

\author{P. Aboyoun}

\seealso{
  \code{\link{pairwiseAlignment}},
  \code{\link{AlignedXStringSet-class}},
  \code{\link{XString-class}},
  \code{\link{XStringViews-class}},
  \code{\link{align-utils}},
  \code{\link{pid}}
}

\examples{
  PairwiseAlignedXStringSet("-PA--W-HEAE", "HEAGAWGHE-E")
  pattern <- AAStringSet(c("HLDNLKGTF", "HVDDMPNAL"))
  subject <- AAString("SMDDTEKMSMKL")
  nw1 <- pairwiseAlignment(pattern, subject, substitutionMatrix = "BLOSUM50",
    gapOpening = -3, gapExtension = -1)
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
