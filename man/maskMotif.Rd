\name{maskMotif}

\alias{maskMotif}
\alias{maskMotif,MaskedXString,XString-method}
\alias{maskMotif,MaskedXString,character-method}
\alias{maskMotif,XString,ANY-method}

\alias{mask}


\title{Masking by content (or by position)}

\description{
  Functions for masking a sequence by content (or by position).
}

\usage{
maskMotif(x, motif, min.block.width=1, ...)
mask(x, start=NA, end=NA, pattern)
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
  \item{...}{
    Additional arguments for \code{matchPattern}.
  }
  \item{start}{
    An integer vector containing the starting positions
    of the regions to mask.
  }
  \item{end}{
    An integer vector containing the ending positions
    of the regions to mask.
  }
  \item{pattern}{
    The motif to mask in the sequence.
  }
}

\value{
  A \link{MaskedXString} object for \code{maskMotif}
  and an \link{XStringViews} object for \code{mask}.
}

\author{H. Pages}

\seealso{
  \code{\link[IRanges]{read.Mask}},
  \code{\link{matchPattern}},
  \link{XString-class},
  \link{MaskedXString-class},
  \link{XStringViews-class},
  \link[IRanges]{MaskCollection-class}
}

\examples{
  ## ---------------------------------------------------------------------
  ## EXAMPLE 1
  ## ---------------------------------------------------------------------

  maskMotif(BString("AbcbbcbEEE"), "bcb")
  maskMotif(BString("AbcbcbEEE"), "bcb")

  ## maskMotif() can be used in an incremental way to mask more than 1
  ## motif. Note that maskMotif() does not try to mask again what's
  ## already masked (i.e. the new mask will never overlaps with the
  ## previous masks) so the order in which the motifs are masked actually
  ## matters as it will affect the total set of masked positions.
  x0 <- BString("AbcbEEEEEbcbbEEEcbbcbc")
  x1 <- maskMotif(x0, "E")
  x1
  x2 <- maskMotif(x1, "bcb")
  x2
  x3 <- maskMotif(x2, "b")
  x3
  ## Note that inverting the order in which "b" and "bcb" are masked would
  ## lead to a different final set of masked positions.
  ## Also note that the order doesn't matter if the motifs to mask don't
  ## overlap (we assume that the motifs are unique) i.e. if the prefix of
  ## each motif is not the suffix of any other motif. This is of course
  ## the case when all the motifs have only 1 letter.

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

  library(BSgenome.Dmelanogaster.UCSC.dm3)
  chrU <- Dmelanogaster$chrU
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

  ## ---------------------------------------------------------------------
  ## EXAMPLE 4
  ## ---------------------------------------------------------------------
  ## Note that unlike maskMotif(), mask() returns an XStringViews object!

  ## masking "by position"
  mask("AxyxyxBC", 2, 6)

  ## masking "by content"
  mask("AxyxyxBC", "xyx")
  noN_chrU <- mask(chrU, "N")
  noN_chrU
  alphabetFrequency(noN_chrU, collapse=TRUE)
}

\keyword{methods}
\keyword{manip}
