\name{align-utils}

\alias{compareStrings}
\alias{compareStrings,character,character-method}
\alias{compareStrings,XString,XString-method}
\alias{compareStrings,XStringSet,XStringSet-method}
\alias{compareStrings,AlignedXStringSet,AlignedXStringSet-method}
\alias{compareStrings,PairwiseAlignment,missing-method}
\alias{consensusMatrix}
\alias{consensusMatrix,character-method}
\alias{consensusMatrix,matrix-method}
\alias{consensusMatrix,list-method}
\alias{consensusMatrix,XStringSet-method}
\alias{consensusMatrix,XStringViews-method}
\alias{consensusMatrix,PairwiseAlignment-method}
\alias{consensusString}
\alias{consensusString,ANY-method}

% Deprecated:
\alias{consmat}
\alias{consmat,ANY-method}


\title{Utility functions related to sequence alignment}

\description{
  A variety of different functions used to deal with sequence alignments.
}

\usage{
  compareStrings(pattern, subject)
  \S4method{consensusMatrix}{character}(x, freq=FALSE)
  \S4method{consensusMatrix}{XStringSet}(x, baseOnly=FALSE, freq=FALSE)
  consensusString(x)
}

\details{
  \code{compareStrings} combines two equal-length strings that are assumed to be aligned
  into a single character string containing that replaces mismatches with \code{"?"},
  insertions with \code{"+"}, and deletions with \code{"-"}.

  \code{consensusMatrix} computes a consensus matrix for a set of equal-length strings that
  are assumed to be aligned.

  \code{consensusString} creates the string based on a 50\% + 1 vote from the consensus
  matrix with unknowns labeled with \code{"?"}.
}

\arguments{
  \item{pattern, subject}{
    The strings to compare. Can be of type \code{character}, \code{XString},
    \code{XStringSet}, \code{AlignedXStringSet}, or, in the case of
    \code{pattern}, \code{PairwiseAlignment}. If \code{pattern} is a
    \code{PairwiseAlignment} object, then \code{subject} must be missing.
  }
  \item{x}{
    A \code{character} vector or matrix, \code{XStringSet}, \code{XStringViews},
    \code{PairwiseAlignment}, or \code{list} of FASTA records containing the equal-length
    strings.
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
  \item{freq}{
    If \code{TRUE}, then letter frequencies (per position) are reported, otherwise counts.
  }
}

\seealso{
  \link{XString-class}, \link{XStringSet-class}, \link{XStringViews-class}, 
  \link{AlignedXStringSet-class}, \link{PairwiseAlignment-class}
}

\examples{
  ## Compare two globally aligned strings
  string1 <- "ACTTCACCAGCTCCCTGGCGGTAAGTTGATC---AAAGG---AAACGCAAAGTTTTCAAG"
  string2 <- "GTTTCACTACTTCCTTTCGGGTAAGTAAATATATAAATATATAAAAATATAATTTTCATC"
  compareStrings(string1, string2)

  ## Create a consensus matrix
  nw1 <-
    pairwiseAlignment(AAStringSet(c("HLDNLKGTF", "HVDDMPNAL")), AAString("SMDDTEKMSMKL"),
      substitutionMatrix = "BLOSUM50", gapOpening = -3, gapExtension = -1)
  consensusMatrix(nw1)

  ## Examine the consensus between the bacteriophage phi X174 genomes
  data(phiX174Phage)
  phageConsmat <- consensusMatrix(phiX174Phage, baseOnly = TRUE)
  phageDiffs <- which(apply(phageConsmat, 2, max) < length(phiX174Phage))
  phageDiffs
  phageConsmat[,phageDiffs]

  ## Read in ORF data
  file <- system.file("extdata", "someORF.fa", package="Biostrings")
  orf <- read.DNAStringSet(file, "fasta")

  ## To illustrate, the following example assumes the ORF data
  ## to be aligned for the first 10 positions (patently false):
  orf10 <- DNAStringSet(orf, end=10)
  consensusMatrix(orf10, baseOnly=TRUE, freq=TRUE)
  consensusString(sort(orf10)[1:5])

  ## For the character matrix containing the "exploded" representation
  ## of the strings, do:
  as.matrix(orf10, use.names=FALSE)
}

\keyword{methods}

