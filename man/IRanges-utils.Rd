\name{IRanges-utils}

\alias{IRanges-utils}

\alias{intToRanges}
\alias{intToAdjacentRanges}
\alias{whichAsRanges}
\alias{slice}
\alias{shift}
\alias{shift,IRanges-method}
\alias{restrict}
\alias{restrict,IRanges-method}
\alias{narrow}
\alias{narrow,IRanges-method}
\alias{narrow,NormalIRanges-method}
\alias{narrow,numeric-method}
\alias{reduce}
\alias{reduce,IRanges-method}
\alias{toNormalIRanges}
\alias{gaps}
\alias{gaps,IRanges-method}


\title{IRanges utility functions}

\description{
  Utility functions for creating or modifying \link{IRanges} objects.
}

\usage{
  ## Create an IRanges object
  intToRanges(x, use.names=TRUE)
  intToAdjacentRanges(x, use.names=TRUE)
  whichAsRanges(x)
  slice(x, lower=-Inf, upper=Inf, includeLower=TRUE, includeUpper=TRUE)
  toNormalIRanges(x)

  ## Modify an IRanges object
  shift(x, shift, use.names=TRUE)
  restrict(x, start=NA, end=NA, keep.all.ranges=FALSE, use.names=TRUE)
  narrow(x, start=NA, end=NA, width=NA, use.names=TRUE)
  reduce(x, with.inframe.attrib=FALSE)
  gaps(x, start=NA, end=NA)
}

\arguments{
  \item{x}{
    An integer vector for \code{intToRanges} and \code{intToAdjacentRanges}.
    A logical vector for \code{whichAsRanges}.
    A numeric vector for \code{slice}.
    An \link{IRanges} object for \code{shift}, \code{restrict}, \code{reduce},
    \code{toNormalIRanges} and \code{gaps}.
    An \link{IRanges} object or an integer vector for \code{narrow}.
  }
  \item{use.names}{
    \code{TRUE} or \code{FALSE}. Should names be preserved?
  }
  \item{lower, upper}{
    The lower and upper bounds for the slice.
  }
  \item{includeLower, includeUpper}{
    Logical values that indicate whether or not to include the equality in the
    corresponding bound for the slice.
  }
  \item{shift}{
    A single integer.
  }
  \item{start}{
    A single integer or \code{NA}. The \code{narrow} function also excepts an
    integer vector of \code{length(x)}.
  }
  \item{end}{
    A single integer or \code{NA}. The \code{narrow} function also excepts an
    integer vector of \code{length(x)}.
  }
  \item{keep.all.ranges}{
    \code{TRUE} or \code{FALSE}. Should ranges that become "out of limits"
    after restriction be kept?
  }
  \item{width}{
    A single number or \code{NA}. The \code{narrow} function also excepts an
    integer vector of \code{length(x)}.
  }
  \item{with.inframe.attrib}{
    \code{TRUE} or \code{FALSE}. For internal use.
  }
}

\details{
  \code{intToRanges} returns an \link{IRanges} object of the same length
  as \code{x} where each integer in \code{x} has been converted into the
  range that starts at 1 and ends at that integer.

  \code{intToAdjacentRanges} returns an \link{IRanges} object of the same
  length as \code{x} where the ranges are adjacent and have the width of
  the corresponding integer in \code{x}. The first range in the returned
  \link{IRanges} object starts at 1.

  \code{whichAsRanges} returns an \link{IRanges} object containing all of
  the ranges where \code{x} is \code{TRUE}.

  \code{toNormalIRanges} first reduces \code{x}, then removes the
  empty ranges. The final result is returned as a \link{NormalIRanges}
  object.

  \code{shift} shifts all the ranges in \code{x}.

  \code{restrict} restricts the ranges in \code{x} to the interval
  specified by the \code{start} and \code{end} arguments.

  \code{narrow} narrows the ranges in \code{x} i.e. each range in the
  returned \link{IRanges} object is a subrange of the corresponding
  range in \code{x}. This subrange is determined by the \code{start},
  \code{end} and \code{width} arguments.
  The \code{start}, \code{end}, and \code{width} arguments must be
  numeric vectors of length 1 or \code{length(x)} with at least one of
  the corresponding \code{start}, \code{end} and \code{width} values being
  \code{NA}. If at least two of them are \code{NA}s, then \code{start=NA}
  means that the starts are unchanged and \code{end=NA} means that the ends
  are unchanged.
  A positive \code{start} or \code{end} is interpreted relatively to the
  start of the original range whereas a negative \code{start} or \code{end}
  is interpreted relatively to its end.
  Note that, if \code{width} is \code{NA}, then \code{start=NA} is
  equivalent to \code{start=1} and \code{end=NA} is equivalent to
  \code{end=-1}.
  Finally, if \code{width} is not \code{NA}, then \code{start} and
  \code{end} cannot be both \code{NA}s.

  \code{reduce} first orders the ranges in \code{x} from left to right,
  then merges the overlapping or adjacent ones.

  \code{gaps} returns the normal \link{IRanges} object describing
  the set of integers obtained by removing the set of integers described
  by \code{x} from the interval specified by the \code{start} and
  \code{end} arguments.
}

\author{H. Pages}

\seealso{
  \link{IRanges-class},
  \link{MaskCollection-class},
  \code{\link{reverse}}
  \code{\link{coverage}}
}

\examples{
  vec <- c(19, 5, 0, 8, 5)
  intToRanges(vec)
  whichAsRanges(vec >= 5)
  slice(vec, lower = 5, upper = 8)
  x <- intToAdjacentRanges(vec)
  x
  shift(x, -3)
  restrict(x, start=12, end=34)
  y <- x[width(x) != 0]
  narrow(y, start=4, end=-2)
  narrow(y, start=-4, end=-2)
  narrow(y, end=5, width=3)
  narrow(y, start = c(3, 4, 2, 3), end = c(12, 5, 7, 4))

  x <- IRanges(start=c(-2L, 6L, 9L, -4L, 1L, 0L, -6L, 10L),
               width=c( 5L, 0L, 6L,  1L, 4L, 3L,  2L,  3L))
  reduce(x)
  toNormalIRanges(x)   # 3 ranges ordered from left to right and separated by
                       # gaps of width >= 1.
  gaps(x)
  gaps(x, start=-6, end=20)  # Regions of the -6:20 range that are not masked by 'x'.
}

\keyword{utilities}
