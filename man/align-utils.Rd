\name{align-utils}

\alias{align-utils}

\alias{mismatch,AlignedXStringSet0,missing-method}
\alias{nmatch,PairwiseAlignments,missing-method}
\alias{nmatch,PairwiseAlignmentsSingleSubjectSummary,missing-method}
\alias{nmismatch,AlignedXStringSet0,missing-method}
\alias{nmismatch,PairwiseAlignments,missing-method}
\alias{nmismatch,PairwiseAlignmentsSingleSubjectSummary,missing-method}
\alias{nedit}
\alias{nedit,PairwiseAlignments-method}
\alias{nedit,PairwiseAlignmentsSingleSubjectSummary-method}

\alias{mismatchTable}
\alias{mismatchTable,AlignedXStringSet0-method}
\alias{mismatchTable,QualityAlignedXStringSet-method}
\alias{mismatchTable,PairwiseAlignments-method}
\alias{mismatchSummary}
\alias{mismatchSummary,AlignedXStringSet0-method}
\alias{mismatchSummary,QualityAlignedXStringSet-method}
\alias{mismatchSummary,PairwiseAlignmentsSingleSubject-method}
\alias{mismatchSummary,PairwiseAlignmentsSingleSubjectSummary-method}

\alias{coverage,AlignedXStringSet0-method}
\alias{coverage,PairwiseAlignmentsSingleSubject-method}
\alias{coverage,PairwiseAlignmentsSingleSubjectSummary-method}

\alias{compareStrings}
\alias{compareStrings,character,character-method}
\alias{compareStrings,XString,XString-method}
\alias{compareStrings,XStringSet,XStringSet-method}
\alias{compareStrings,AlignedXStringSet0,AlignedXStringSet0-method}
\alias{compareStrings,PairwiseAlignments,missing-method}
\alias{consensusMatrix,PairwiseAlignmentsSingleSubject-method}


\title{Utility functions related to sequence alignment}

\description{
  A variety of different functions used to deal with sequence alignments.
}

\usage{
nedit(x) # also nmatch and nmismatch

mismatchTable(x, shiftLeft=0L, shiftRight=0L, \dots)
mismatchSummary(x, \dots)

\S4method{coverage}{AlignedXStringSet0}(x, shift=0L, width=NULL, weight=1L)

\S4method{coverage}{PairwiseAlignmentsSingleSubject}(x, shift=0L, width=NULL, weight=1L)

compareStrings(pattern, subject)

\S4method{consensusMatrix}{PairwiseAlignmentsSingleSubject}(x,
                as.prob=FALSE, shift=0L, width=NULL,
                baseOnly=FALSE, gapCode="-", endgapCode="-")
}

\details{
  \code{mismatchTable}:  a data.frame containing the positions and substrings
  of the mismatches for the \code{AlignedXStringSet} or
  \code{PairwiseAlignments} object.

  \code{mismatchSummary}:  a list of data.frame objects containing counts and
  frequencies of the mismatches for the \code{AlignedXStringSet} or
  \code{PairwiseAlignmentsSingleSubject} object.

  \code{compareStrings} combines two equal-length strings that are assumed to be
  aligned into a single character string containing that replaces mismatches
  with \code{"?"}, insertions with \code{"+"}, and deletions with \code{"-"}.
}

\arguments{
  \item{x}{
    A \code{character} vector or matrix, \code{XStringSet}, \code{XStringViews},
    \code{PairwiseAlignments}, or \code{list} of FASTA records containing
    the equal-length strings.
  }
  \item{shiftLeft, shiftRight}{
    Non-positive and non-negative integers respectively that specify how many
    preceding and succeeding characters to and from the mismatch position to
    include in the mismatch substrings.
  }
  \item{\dots}{
    Further arguments to be passed to or from other methods.
  }
  \item{shift, width}{
    See \code{?\link[IRanges]{coverage}}.
  }
  \item{weight}{
    An integer vector specifying how much each element in \code{x} counts.
  }
  \item{pattern, subject}{
    The strings to compare. Can be of type \code{character}, \code{XString},
    \code{XStringSet}, \code{AlignedXStringSet}, or, in the case of
    \code{pattern}, \code{PairwiseAlignments}.

    If the first argument of \code{compareStrings()} (\code{pattern})
    is a \code{PairwiseAlignments} object, then the second argument
    (\code{subject}) must be missing. In this case \code{compareStrings(x)}
    is equivalent to \code{compareStrings(pattern(x), subject(x))}.
  }
  \item{as.prob}{
    If \code{TRUE} then probabilities are reported,
    otherwise counts (the default).
  }
  \item{baseOnly}{
    \code{TRUE} or \code{FALSE}.
    If \code{TRUE}, the returned vector only contains frequencies for the
    letters in the "base" alphabet i.e. "A", "C", "G", "T" if \code{x}
    is a "DNA input", and "A", "C", "G", "U" if \code{x} is "RNA input".
    When \code{x} is a \link{BString} object (or an \link{XStringViews}
    object with a \link{BString} subject, or a \link{BStringSet} object),
    then the \code{baseOnly} argument is ignored.
  }
  \item{gapCode, endgapCode}{
    The codes in the appropriate \code{\link{alphabet}} to use for the internal
    and end gaps.
  }
}

\seealso{
  \code{\link{pairwiseAlignment}},
  \code{\link{consensusMatrix}},
  \link{XString-class}, \link{XStringSet-class}, \link{XStringViews-class},
  \link{AlignedXStringSet-class}, \link{PairwiseAlignments-class},
  \link{match-utils}
}

\examples{
  ## Compare two globally aligned strings
  string1 <- "ACTTCACCAGCTCCCTGGCGGTAAGTTGATC---AAAGG---AAACGCAAAGTTTTCAAG"
  string2 <- "GTTTCACTACTTCCTTTCGGGTAAGTAAATATATAAATATATAAAAATATAATTTTCATC"
  compareStrings(string1, string2)

  ## Create a consensus matrix
  nw1 <-
    pairwiseAlignment(AAStringSet(c("HLDNLKGTF", "HVDDMPNAL")), AAString("SMDDTEKMSMKL"),
      substitutionMatrix = "BLOSUM50", gapOpening = 3, gapExtension = 1)
  consensusMatrix(nw1)

  ## Examine the consensus between the bacteriophage phi X174 genomes
  data(phiX174Phage)
  phageConsmat <- consensusMatrix(phiX174Phage, baseOnly = TRUE)
  phageDiffs <- which(apply(phageConsmat, 2, max) < length(phiX174Phage))
  phageDiffs
  phageConsmat[,phageDiffs]
}

\keyword{methods}
