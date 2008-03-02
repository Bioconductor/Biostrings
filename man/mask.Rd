\name{mask}

\alias{mask}
\alias{mask,IntIntervals-method}
\alias{mask,BStringViews-method}
\alias{mask,BString-method}
\alias{mask,character-method}


\title{Masking regions in a sequence}

\description{
  Use this function to mask regions in a sequence.
  The regions to mask can either be specified "by position" (start/end)
  or "by content" (they match a given pattern).
}

\usage{
  mask(x, start, end, ...)
}

\arguments{
  \item{x}{
    The single string or \link{BString} (or derived) object to mask.
    This can also be a \link{BStringViews} object.
  }
  \item{start}{
    [TODO]
  }
  \item{end}{
    [TODO]
  }
  \item{...}{
    Additional arguments to be passed to or from methods.
  }
}

\details{
  For \link{BString} (or derived) objects, \code{mask()} can be used in two
  different ways: by specifying the starts/ends of the regions to mask
  (masking "by position") or by specifying a pattern (masking "by content").

  For \link{BStringViews} objects, no extra argument is needed: \code{y <- mask(x)}
  is the shortest (in term of number of views) \link{BStringViews} object such that:
    (a) subject(y) == subject(x),
    (b) y views cover the regions of the subject that are not covered by x
        views,
    (c) y views are sorted from left to right,
    (d) y views are not "out of limits".
}

\value{
  A \link{BStringViews} object containing the set of "unmasked" regions.
}

\author{H. Pages}

\seealso{
  \link{BStringViews},
  \link{BString},
  \link{DNAString},
  \link{RNAString},
  \link{IntIntervals-utils}
}

\examples{
  ## masking "by position"
  mask("AxyxyxBC", 2, 6)

  ## masking "by content"
  mask("AxyxyxBC", "xyx")

  ## masking the N's in a chromosome sequence
  library(BSgenome.Dmelanogaster.FlyBase.r51)
  chrU <- Dmelanogaster[["U"]]
  alphabetFrequency(chrU)
  noN_chrU <- mask(chrU, "N")
  alphabetFrequency(noN_chrU)

  ## masking a BStringViews object
  mask(noN_chrU)
  v <- views("AxyxyxBC", c(6, 4,-1, NA), c(7, 6, 1, 1))
  mask(v)

  ## note that applying mask() again returns a BStringViews object
  ## where the original views on v have been merged, sorted from left to right
  ## and cleaned from their "out of limits" regions
  mask(mask(v))

  ## finally, note that mask(mask(mask(v))) is _always_ the same as mask(v)
  mask(mask(mask(v)))
}

\keyword{methods}

