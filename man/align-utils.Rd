\name{align-utils}

\alias{compareStrings}
\alias{compareStrings,character,character-method}
\alias{compareStrings,AlignedXStringSet,AlignedXStringSet-method}
\alias{compareStrings,PairwiseAlignment,missing-method}
\alias{consmat}
\alias{consmat,character-method}
\alias{consmat,matrix-method}
\alias{consmat,list-method}
\alias{consmat,XStringSet-method}
\alias{consmat,XStringViews-method}
\alias{consmat,PairwiseAlignment-method}


\title{Utility functions related to sequence alignment}

\description{
  A variety of different functions used to deal with sequence alignments.
}

\usage{
  compareStrings(pattern, subject)
  consmat(x, freq=TRUE)
}

\details{
  \code{compareStrings} combines two equal-length strings that are assumed to be aligned
  into a single character string containing that replaces mismatches with \code{"?"},
  insertions with \code{"+"}, and deletions with \code{"-"}.

  \code{consmat} computes a consensus matrix for a set of equal-length strings that
  are assumed to be aligned.
}

\arguments{
  \item{pattern, subject}{
    The strings to compare. Can be of type \code{character}, \code{AlignedXStringSet}, or,
    in the case of \code{pattern}, \code{PairwiseAlignment}. If \code{pattern} is a
    \code{PairwiseAlignment} object, then \code{subject} must be missing.
  }
  \item{x}{
    A \code{character} vector, \code{XStringSet}, \code{XStringView}, \code{PairwiseAlignment},
    or \code{list} of FASTA records containing the equal-length strings.
  }
  \item{freq}{
    If \code{TRUE}, then letter frequencies (per position) are reported, otherwise counts.
  }
}

\seealso{
  \link{XStringSet-class}, \link{XStringViews-class}, \link{PairwiseAlignment-class}
}

\examples{
  ## Compare two globally aligned strings
  string1 <- "ACTTCACCAGCTCCCTGGCGGTAAGTTGATC---AAAGG---AAACGCAAAGTTTTCAAG"
  string2 <- "GTTTCACTACTTCCTTTCGGGTAAGTAAATATATAAATATATAAAAATATAATTTTCATC"
  compareStrings(string1, string2)

  ## Read in ORF data
  file <- system.file("extdata", "someORF.fa", package="Biostrings")
  orf <- read.DNAStringSet(file, "fasta")

  ## To illustrate, the following example assumes the ORF data
  ## to be aligned for the first 10 positions (patently false):
  orf10 <- DNAStringSet(orf, end=10)
  consmat(orf10)

  ## For the character matrix containing the "exploded" representation
  ## of the views, do:
  as.matrix(orf10, mode="character")
}

\keyword{methods}

