\name{align-utils}

\alias{consmat}
\alias{consmat,matrix-method}
\alias{consmat,character-method}
\alias{consmat,list-method}
\alias{consmat,BStringViews-method}

\title{Utility functions related to sequence alignment}

\description{
  [TODO]
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
    [TODO]
  }
  \item{freq}{
    [TODO]
  }
}

\seealso{
  \link{BStringViews-class}
}

\examples{
  file <- system.file("Exfiles", "someORF.fsa", package="Biostrings")
  orf <- read.BStringViews(file, "fasta", "DNAString")

  ## To illustrate, the following example assumes the ORF data
  ## to be aligned for the first 10 positions (patently false):
  orf10 <- subviews(orf, 1, 10)
  consmat(orf10)

  ## For the character matrix containing the "exploded" representation
  ## of the views, do:
  as.matrix(orf10, mode="character")
}

\keyword{methods}

