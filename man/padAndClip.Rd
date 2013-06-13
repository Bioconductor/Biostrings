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
padAndClip(x, views, Lpadding.letter=" ", Rpadding.letter=" ")
stackStrings(x, from, to, shift=0L, Lpadding.letter=" ", Rpadding.letter=" ")
}

\arguments{
  \item{x}{
    An \link{XStringSet} object containing the strings to pad and clip.
  }
  \item{views}{
    A \link[IRanges]{Ranges} object (recycled to the length of \code{x} if
    necessary) defining the region to keep for each string. Because the
    strings are first conceptually padded with an infinite number of
    padding letters on both sides, regions can go beyond string limits.
  }
  \item{Lpadding.letter}{
    A single letter to use for padding on the left.
  }
  \item{Rpadding.letter}{
    A single letter to use for padding on the right.
  }
  \item{from, to}{
    Another way to specify the region to keep for each string, but with
    the restriction that \code{from} and \code{to} must be single integers.
    So only 1 region can be specified, and the same region is used for all
    the strings.
  }
  \item{shift}{
    An integer vector (recycled to the length of \code{x} if necessary)
    specifying the amount of shifting (in number of letters) to apply to
    each string before doing pad and clip. Positive values shift to the
    right and negative values to the left.
  }
}

\value{
  For \code{padAndClip}: An \link{XStringSet} object of the same length
  and names as \code{x}. Its "shape" (which is described by the integer
  vector returned by \code{width}) is the same as the shape of the
  \code{views} argument after recycling.

  The class of the returned object is always a direct concrete subclass of
  the \link{XStringSet} virtual class, that is, depending on the class of
  \code{x}, it is one of the 4 following classes: \link{BStringSet},
  \link{DNAStringSet}, \link{RNAStringSet}, and \link{AAStringSet}.
  More precisely, if \code{x} is an \emph{instance} of one of those 4 classes,
  then the returned object has the same class as \code{x} (i.e. in that case,
  \code{padAndClip} acts as an endomorphism). But if \code{x} \emph{derives}
  from one of those 4 classes, then the returned object is downgraded to
  the class \code{x} derives from. In that case, \code{padAndClip} doesn't
  act as an endomorphism.

  For \code{stackStrings}: Same as \code{padAndClip}. In addition it is
  guaranteed to have a rectangular shape i.e. to be a constant-width
  \link{XStringSet} object.
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
x <- BStringSet(c(seq1="ABCD", seq2="abcdefghijk", seq4="", seq3="XYZ"))
padAndClip(x, IRanges(3, 8:5), Lpadding.letter="+", Rpadding.letter="-")
padAndClip(x, IRanges(1:-2, 7), Lpadding.letter="+", Rpadding.letter="-")
stackStrings(x, from=-2, to=8)
stackStrings(x, from=-2, to=8, shift=c(0, -11, 6, 7),
             Lpadding.letter="+", Rpadding.letter="-")
}

\keyword{methods}
\keyword{manip}
