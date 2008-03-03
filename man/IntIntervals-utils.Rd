\name{IntIntervals-utils}

\alias{IntIntervals-utils}

\alias{intToIntervals}
\alias{intToAdjacentIntervals}
\alias{restrict}
\alias{restrict,IntIntervals-method}
\alias{narrow}
\alias{narrow,IntIntervals-method}
\alias{narrow,numeric-method}
\alias{reduce}
\alias{reduce,IntIntervals-method}


\title{IntIntervals utility functions}

\description{
  Utility functions for creating or modifying \link{IntIntervals} objects.
}

\usage{
  intToIntervals(x, use.names=TRUE)
  intToAdjacentIntervals(x, use.names=TRUE)
  restrict(x, start, end, keep.nonoverlapping=FALSE, use.names=TRUE)
  narrow(x, start=NA, end=NA, width=NA, use.names=TRUE)
  reduce(x, with.inframe.attrib=FALSE)
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
  \link{IntIntervals-class},
  \link{mask},
  \link{BStringViews}
}

\examples{
  intToIntervals(c(19, 5, 0, 8, 5))
  x <- intToAdjacentIntervals(c(19, 5, 0, 8, 5))
  x
  restrict(x, 12, 34)
  y <- x[width(x) != 0]
  narrow(y, start=4, end=-2)
  narrow(y, start=-4, end=-2)
  narrow(y, end=5, width=3)

  x <- new("IntIntervals",
           start=c(-2L,  9L, -4L, 1L, 0L, -6L, 10L),
           width=c( 5L,  6L,  1L, 4L, 3L,  2L,  3L))
  reduce(x)  # 3 intervals ordered from left to right and separated by gaps
             # of width >= 1
  mask(x, -6, 20)  # regions of the [-6, 20] interval that are not maked by
                   # 'x'
}

\keyword{utilities}
