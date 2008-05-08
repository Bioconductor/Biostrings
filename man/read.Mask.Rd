\name{read.Mask}

\alias{read.Mask}
\alias{read.rmMask}
\alias{read.trfMask}


\title{Read a mask from a file}

\description{
  Two functions, \code{read.rmMask} and \code{read.trfMask}, for reading
  a mask from a RepeatMasker .out file or from a Tandem Repeats Finder .bed
  file.
}

\usage{
  read.rmMask(file, width, use.IDs=FALSE)
  read.trfMask(file, width)
}

\arguments{
  \item{file}{
    Either a character string naming a file or a connection open
    for reading.
  }
  \item{width}{
    The width of the mask to return.
    See \code{?`\link{MaskCollection-class}`} for more information about
    the width of a \link{MaskCollection} object.
  }
  \item{use.IDs}{
    Whether or not the repeat IDs provided in the RepeatMasker .out file
    should be used to name the ranges constituing the returned mask.
    See \code{?`\link{IRanges-class}`} for more information about the
    names of an \link{IRanges} object.
  }
}

\details{
  Only RepeatMasker .out files and Tandem Repeats Finder .bed files
  are supported for now.
}

\seealso{
  \code{\link{maskMotif}},
  \link{MaskCollection-class},
  \link{MaskedXString-class},
  \link{IRanges-class}
}

\examples{
  ## Before you can read a mask, you need to know the length of the
  ## sequence that you're going to put the mask on:
  library(BSgenome.Celegans.UCSC.ce2)
  chrM <- Celegans$chrM
  mask_width <- length(chrM)

  ## Read the RepeatMasker .out file for chrM in ce2:
  file1 <- system.file("extdata", "ce2chrM.fa.out", package="Biostrings")
  mask1 <- read.rmMask(file1, width=mask_width)
  mask1
  names(mask1) <- "RepeatMasker"
  mask1

  ## Read the Tandem Repeats Finder .bed file for chrM in ce2:
  file2 <- system.file("extdata", "ce2chrM.bed", package="Biostrings")
  mask2 <- read.trfMask(file2, width=mask_width)
  mask2
  names(mask2) <- "Tandem Repeats Finder (period<=12)"
  mask2

  ## Put the 2 masks on chrM:
  chrM <- NULL  # drop current masks (if any)
  masks(chrM) <- mask1
  masks(chrM) <- append(masks(chrM), mask2)
  chrM
}

\keyword{manip}
