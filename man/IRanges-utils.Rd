\name{IRanges-utils}

\alias{IRanges-utils}

\alias{intToRanges}
\alias{intToAdjacentRanges}
\alias{restrict}
\alias{restrict,IRanges-method}
\alias{narrow}
\alias{narrow,.IRanges-method}
\alias{narrow,numeric-method}
\alias{reduce}
\alias{reduce,.IRanges-method}
\alias{toNormalIRanges}
\alias{gaps}
\alias{gaps,IRanges-method}


\title{IRanges utility functions}

\description{
  Utility functions for creating or modifying \link{IRanges} objects.
}

\usage{
  intToRanges(x, use.names=TRUE)
  intToAdjacentRanges(x, use.names=TRUE)
  restrict(x, start, end, keep.all.ranges=FALSE, use.names=TRUE)
  narrow(x, start=NA, end=NA, width=NA, use.names=TRUE)
  reduce(x, with.inframe.attrib=FALSE)
  toNormalIRanges(x)
  gaps(x, start=NA, end=NA)
}

\arguments{
  \item{x}{
    [TODO]
  }
  \item{use.names}{
    [TODO]
  }
  \item{start}{
    [TODO]
  }
  \item{end}{
    [TODO]
  }
  \item{keep.all.ranges}{
    [TODO]
  }
  \item{width}{
    [TODO]
  }
  \item{with.inframe.attrib}{
    [TODO]
  }
}

\details{
  [TODO]
}

\author{H. Pages}

\seealso{
  \link{IRanges-class},
  \link{XStringViews-class}
}

\examples{
  intToRanges(c(19, 5, 0, 8, 5))
  x <- intToAdjacentRanges(c(19, 5, 0, 8, 5))
  x
  restrict(x, 12, 34)
  y <- x[width(x) != 0]
  narrow(y, start=4, end=-2)
  narrow(y, start=-4, end=-2)
  narrow(y, end=5, width=3)

  x <- new("IRanges",
           start=c(-2L,  6L, 9L, -4L, 1L, 0L, -6L, 10L),
           width=c( 5L,  0L, 6L,  1L, 4L, 3L,  2L,  3L))
  reduce(x)
  toNormalIRanges(x)  # 3 ranges ordered from left to right and separated by
                      # gaps of width >= 1
  gaps(x)
  gaps(x, -6, 20)  # regions of the -6:20 range that are not masked by
                   # 'x'
}

\keyword{utilities}
