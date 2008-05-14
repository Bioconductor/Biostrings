\name{read.Mask}

\alias{read.Mask}
\alias{read.agpMask}
\alias{read.liftMask}
\alias{read.rmMask}
\alias{read.trfMask}


\title{Read a mask from a file}

\description{
  Four functions -- \code{read.agpMask}, \code{read.liftMask},
  \code{read.rmMask} and \code{read.trfMask} -- for extracting a mask from
  an NCBI "agp" file, an UCSC "lift" file (i.e. a file containing offsets of
  contigs within sequences), a RepeatMasker .out file, or a Tandem Repeats
  Finder .bed file.
}

\usage{
  read.agpMask(file, width, seqname="?", gap.types=NULL, use.gap.types=TRUE)
  read.liftMask(file, seqname="?", width=NA)
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
  \item{seqname}{
    The name of the sequence for which the mask must be extracted.
    If no sequence is specified (i.e. \code{seqname="?"}) then an error is
    raised and the sequence names found in the file are displayed.
    If the file doesn't contain any information for the specified sequence,
    then a warning is issued and an empty mask of width \code{width}
    is returned.
  }
  \item{gap.types}{
    \code{NULL} or a character vector containing gap types.
    Use this argument to filter the N-gaps that are to be extracted from
    the "agp" file based on their type. Common types are \code{"contig"},
    \code{"clone"}, \code{"centromere"} and \code{"heterochromatin"}.
    With \code{gap.types=NULL}, all the N-gaps described in the file are
    extracted.
    With \code{gap.types="?"}, an error is raised and the gap types found
    in the file for the specified sequence are displayed.
  }
  \item{use.gap.types}{
    Whether or not the gap types provided in the "agp" file should be used
    to name the ranges constituing the returned mask.
    See \code{?`\link{IRanges-class}`} for more information about the
    names of an \link{IRanges} object.
  }
  \item{use.IDs}{
    Whether or not the repeat IDs provided in the RepeatMasker .out file
    should be used to name the ranges constituing the returned mask.
    See \code{?`\link{IRanges-class}`} for more information about the
    names of an \link{IRanges} object.
  }
}

\details{
  Only NCBI "agp" file, UCSC "lift" files, RepeatMasker .out files, and
  Tandem Repeats Finder .bed files are supported for now.
}

\seealso{
  \code{\link{maskMotif}},
  \link{MaskCollection-class},
  \link{MaskedXString-class},
  \link{IRanges-class}
}

\examples{
  ## ---------------------------------------------------------------------
  ## A. Extract a mask of "inter-contig gaps" with read.agpMask()
  ## ---------------------------------------------------------------------
  ## Note: The hs_b36v3_chrY.agp file was obtained by downloading,
  ## extracting and renaming the hs_ref_chrY.agp.gz file from
  ##
  ##   ftp://ftp.ncbi.nih.gov/genomes/H_sapiens/Assembled_chromosomes/
  ##     hs_ref_chrY.agp.gz      5 KB  24/03/08  04:33:00 PM
  ##
  ## on May 9, 2008.

  file1 <- system.file("extdata", "hs_b36v3_chrY.agp", package="Biostrings")
  mask1 <- read.agpMask(file1, 57772954, seqname="chrY")
  mask1
  mask1[[1]]

  mask11 <- read.agpMask(file1, 57772954, seqname="chrY",
                         gap.types=c("centromere", "heterochromatin"))
  mask11[[1]]

  ## ---------------------------------------------------------------------
  ## B. Extract a mask of "inter-contig gaps" with read.liftMask()
  ## ---------------------------------------------------------------------
  ## Note: The hg18liftAll.lft file was obtained by downloading,
  ## extracting and renaming the liftAll.zip file from
  ##
  ##   http://hgdownload.cse.ucsc.edu/goldenPath/hg18/bigZips/
  ##     liftAll.zip             03-Feb-2006 11:35  5.5K
  ##
  ## on May 8, 2008.

  file2 <- system.file("extdata", "hg18liftAll.lft", package="Biostrings")
  mask2 <- read.liftMask(file2, seqname="chr1")
  mask2
  if (interactive()) {
      ## contigs 7 and 8 for chrY are adjacent
      read.liftMask(file2, seqname="chrY")

      ## displays the sequence names found in the file
      read.liftMask(file2)

      ## specify an unknown sequence name
      read.liftMask(file2, seqname="chrZ", width=300)
  }

  ## ---------------------------------------------------------------------
  ## C. Extract an "rm" or a "trf" mask with read.rmMask() or
  ##    read.trfMask()
  ## ---------------------------------------------------------------------
  ## Note: The ce2chrM.fa.out and ce2chrM.bed files were obtained by
  ## downloading, extracting and renaming the chromOut.zip and
  ## chromTrf.zip files from
  ##
  ##   http://hgdownload.cse.ucsc.edu/goldenPath/ce2/bigZips/
  ##     chromOut.zip            21-Apr-2004 09:05  2.6M
  ##     chromTrf.zip            21-Apr-2004 09:07  182K
  ##
  ## on May 7, 2008.

  ## Before you can extract a mask with read.rmMask() or read.trfMask(), you
  ## need to know the length of the sequence that you're going to put the
  ## mask on:
  library(BSgenome.Celegans.UCSC.ce2)
  chrM <- Celegans$chrM
  mask_width <- length(chrM)

  ## Read the RepeatMasker .out file for chrM in ce2:
  file3 <- system.file("extdata", "ce2chrM.fa.out", package="Biostrings")
  mask3 <- read.rmMask(file3, mask_width)
  mask3
  names(mask3) <- "RepeatMasker"
  mask3

  ## Read the Tandem Repeats Finder .bed file for chrM in ce2:
  file4 <- system.file("extdata", "ce2chrM.bed", package="Biostrings")
  mask4 <- read.trfMask(file4, mask_width)
  mask4
  names(mask4) <- "Tandem Repeats Finder (period<=12)"
  mask4

  ## Put the 2 masks on chrM:
  masks(chrM) <- mask3  # this would drop all current masks, if any
  masks(chrM) <- append(masks(chrM), mask4)
  chrM
}

\keyword{manip}
