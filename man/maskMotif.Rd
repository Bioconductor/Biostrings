\name{maskMotif}

\alias{maskMotif}
\alias{maskMotif,MaskedXString,XString-method}
\alias{maskMotif,MaskedXString,character-method}
\alias{maskMotif,XString,ANY-method}


\title{Masking by content}

\description{
  Generic function for masking a sequence by content.
}

\usage{
  maskMotif(x, motif, min.block.width=1)
}

\arguments{
  \item{x}{
    The sequence to mask.
  }
  \item{motif}{
    The motif to mask in the sequence.
  }
  \item{min.block.width}{
    The minimum width of the blocks to mask.
  }
}

\author{H. Pages}

\seealso{
  \link{XString-class},
  \link{MaskedXString-class},
  \link{XStringViews-class},
  \link{MaskCollection-class}
}

\examples{
  ## ---------------------------------------------------------------------
  ## EXAMPLE 1
  ## ---------------------------------------------------------------------

  maskMotif(BString("AbcbcbDE"), "bcb")

  ## ---------------------------------------------------------------------
  ## EXAMPLE 2
  ## ---------------------------------------------------------------------

  x <- DNAString("ACACAACTAGATAGNACTNNGAGAGACGC")

  ## Mask the N-blocks
  x1 <- maskMotif(x, "N")
  x1
  as(x1, "XStringViews")
  gaps(x1)
  as(gaps(x1), "XStringViews")

  ## Mask the AC-blocks 
  x2 <- maskMotif(x1, "AC")
  x2
  gaps(x2)

  ## Mask the GA-blocks
  x3 <- maskMotif(x2, "GA", min.block.width=5)
  x3  # masks 2 and 3 overlap
  gaps(x3)

  ## ---------------------------------------------------------------------
  ## EXAMPLE 3
  ## ---------------------------------------------------------------------

  library(BSgenome.Dmelanogaster.FlyBase.r51)
  chrU <- Dmelanogaster[["U"]]
  chrU
  alphabetFrequency(chrU)
  chrU <- maskMotif(chrU, "N")
  chrU
  alphabetFrequency(chrU)
  as(chrU, "XStringViews")
  as(gaps(chrU), "XStringViews")

  mask2 <- Mask(mask.width=length(chrU), start=c(50000, 350000, 543900), width=25000)
  names(mask2) <- "some ugly regions"
  masks(chrU) <- append(masks(chrU), mask2)
  chrU
  as(chrU, "XStringViews")
  as(gaps(chrU), "XStringViews")
}

\keyword{methods}
\keyword{manip}
