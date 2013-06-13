\name{padAndClip}

\alias{padAndClip}
\alias{stackStrings}


\title{Pad and clip strings}

\description{
  \code{padAndClip} first conceptually pads the supplied strings with an
  infinite number of padding letters on both sides, then clip them.

  \code{stackStrings} is a convenience wrapper to \code{padAndClip}
  that turns a variable-width set of strings into a rectangular
  (i.e. constant-width) set, by padding and clipping the strings,
  after conceptually shifting them horizontally.
}

\usage{
padAndClip(x, views, padding.letter)
stackStrings(x, from, to, padding.letter, shift=0L)
}

\arguments{
  \item{x}{
    An \link{XStringSet} object.
  }
  \item{views}{
    TODO
  }
  \item{padding.letter}{
    TODO
  }
  \item{from, to}{
    TODO
  }
  \item{shift}{
    TODO
  }
}

\value{
  For \code{padAndClip}:

  For \code{stackStrings}:
}

\author{H. Pages}

\seealso{
  \itemize{
    \item The \code{\link[Rsamtools]{stackStringsFromBam}} function in
          the \pkg{Rsamtools} package for stacking the read sequences (or
          their quality strings) stored in a BAM file on a region of interest.

    \item The \link{XStringViews} class to formally represent a set of views
          on a single string.

    \item The \code{\link{extractAt}} and \code{\link{replaceAt}} functions
          for extracting/replacing arbitrary substrings from/in a string or
          set of strings.

    \item The \link{XStringSet} class.

    \item The \link[IRanges]{Ranges} class in the \pkg{IRanges} package.
  }
}

\examples{
}

\keyword{methods}
\keyword{manip}
