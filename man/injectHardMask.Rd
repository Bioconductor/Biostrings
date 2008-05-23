\name{injectHardMask}

\alias{injectHardMask}
\alias{injectHardMask,XStringViews-method}
\alias{injectHardMask,MaskedXString-method}


\title{Injecting a hard mask in a sequence}

\description{
  [TODO]
}

\usage{
  injectHardMask(x, letter="+")
}

\arguments{
  \item{x}{
    A \link{MaskedXString} or \link{XStringViews} object.
  }
  \item{letter}{
    A single letter.
  }
}

\value{
  An \link{XString} object of the same length as the orignal object \code{x}
  if \code{x} is a \link{MaskedXString} object, or of the same length
  as \code{subject(x)} if it's an \link{XStringViews} object.
}

\author{H. Pages}

\seealso{
  \code{\link{replaceLetterAtLoc}},
  \code{\link{chartr}},
  \link{XString},
  \link{XStringViews-class},
  \link{MaskedXString-class}
}

\examples{
  mask0 <- Mask(mask.width=29, start=c(3, 10, 25), width=c(6, 8, 5))
  x <- DNAString("ACACAACTAGATAGNACTNNGAGAGACGC")
  masks(x) <- mask0
  x
  injectHardMask(x)
}

\keyword{utilities}
\keyword{manip}
