\name{align-utils}

\alias{consmat}
\alias{consmat,character-method}
\alias{consmat,matrix-method}
\alias{consmat,list-method}
\alias{consmat,XStringSet-method}
\alias{consmat,BStringViews-method}
\alias{consmat,XStringAlign-method}


\title{Utility functions related to sequence alignment}

\description{
  A variety of different functions used to deal with sequence alignments.
}

\usage{
  consmat(x, freq=TRUE)
}

\details{
  \code{consmat} computes a consensus matrix for a set of equal-length strings
  assumed to be aligned.
}

\arguments{
  \item{x}{
    A character vector, or an instance of either the BStringView, or the 
   XStringAlign classes.
  }
  \item{freq}{
    If \code{TRUE} then letter frequencies (per position) are reported, 
    otherwise counts.
  }
}

\seealso{
  \link{BStringViews-class}, \link{XStringAlign}
}

\examples{
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

