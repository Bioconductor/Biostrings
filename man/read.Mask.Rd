\name{read.Mask}

\alias{read.Mask}
\alias{read.liftMask}
\alias{read.rmMask}
\alias{read.trfMask}


\title{Read a mask from a file}

\description{
  Three functions -- \code{read.liftMask}, \code{read.rmMask} and
  \code{read.trfMask} -- for reading a mask from an UCSC "lift" file
  (i.e. a file containing offsets of contigs within sequences),
  a RepeatMasker .out file, or a Tandem Repeats Finder .bed file.
}

\usage{
  read.liftMask(file, seqname=NA, width=NA)
  read.rmMask(file, width, use.IDs=FALSE)
  read.trfMask(file, width)
}

\arguments{
  \item{file}{
    Either a character string naming a file or a connection open
    for reading.
  }
  \item{seqname}{
    The name of the sequence for which the mask of "inter-contig gaps" 
    must be extracted.
    If no sequence is specified (i.e. \code{seqname=NA}) then an error is
    raised and the list of valid sequences is displayed.
    If the file doesn't contain any information for the specified sequence,
    then a warning is issued and an empty mask of width \code{width}
    is returned.
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
  Only UCSC "lift" files, RepeatMasker .out files, and Tandem Repeats Finder
  .bed files are supported for now.
}

\seealso{
  \code{\link{maskMotif}},
  \link{MaskCollection-class},
  \link{MaskedXString-class},
  \link{IRanges-class}
}

\examples{
  ## Read the UCSC "lift" file for hg18 and extract the mask of 
  ## "inter-contig gaps" for chr1:
  file1 <- system.file("extdata", "hg18liftAll.lft", package="Biostrings")
  mask1 <- read.liftMask(file1, seqname="chr1")
  mask1

  ## Before you can read a mask with read.rmMask() or read.trfMask(), you
  ## need to know the length of the sequence that you're going to put the
  ## mask on:
  library(BSgenome.Celegans.UCSC.ce2)
  chrM <- Celegans$chrM
  mask_width <- length(chrM)

  ## Read the RepeatMasker .out file for chrM in ce2:
  file2 <- system.file("extdata", "ce2chrM.fa.out", package="Biostrings")
  mask2 <- read.rmMask(file2, width=mask_width)
  mask2
  names(mask2) <- "RepeatMasker"
  mask2

  ## Read the Tandem Repeats Finder .bed file for chrM in ce2:
  file3 <- system.file("extdata", "ce2chrM.bed", package="Biostrings")
  mask3 <- read.trfMask(file3, width=mask_width)
  mask3
  names(mask3) <- "Tandem Repeats Finder (period<=12)"
  mask3

  ## Put the 2 masks on chrM:
  masks(chrM) <- mask2  # this would drop all current masks, if any
  masks(chrM) <- append(masks(chrM), mask3)
  chrM
}

\keyword{manip}
