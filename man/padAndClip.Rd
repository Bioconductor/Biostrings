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
padAndClip(x, views, Lpadding.letter=" ", Rpadding.letter=" ",
           remove.out.of.view.strings=FALSE)

stackStrings(x, from, to, shift=0L,
             Lpadding.letter=" ", Rpadding.letter=" ",
             remove.out.of.view.strings=FALSE)
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
  \item{Lpadding.letter, Rpadding.letter}{
    A single letter to use for padding on the left, and another one to
    use for padding on the right. Note that the default letter (\code{" "})
    does not work if, for example, \code{x} is a \link{DNAStringSet} object,
    because the space is not a valid DNA letter (see
    \code{?\link{DNA_ALPHABET}}).
    So the \code{Lpadding.letter} and \code{Rpadding.letter} arguments
    \emph{must} be supplied if \code{x} is not a \link{BStringSet} object.
    For example, if \code{x} is a \link{DNAStringSet} object, a typical
    choice is to use \code{"+"}.
  }
  \item{remove.out.of.view.strings}{
    \code{TRUE} or \code{FALSE}. Whether or not to remove the strings
    that are out of view in the returned object.
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
  For \code{padAndClip}: An \link{XStringSet} object.
  If \code{remove.out.of.view.strings} is \code{FALSE}, it has the same
  length and names as \code{x}, and its "shape", which is described by the
  integer vector returned by \code{width()}, is the same as the shape of the
  \code{views} argument after recycling.

  The class of the returned object is the direct concrete subclass of
  \link{XStringSet} that \code{x} belongs to or derives from.
  There are 4 direct concrete subclasses of the \link{XStringSet} virtual
  class: \link{BStringSet}, \link{DNAStringSet}, \link{RNAStringSet},
  and \link{AAStringSet}. If \code{x} is an \emph{instance} of one of
  those classes, then the returned object has the same class as \code{x}
  (i.e. in that case, \code{padAndClip} acts as an endomorphism).
  But if \code{x} \emph{derives} from one of those 4 classes, then the
  returned object is downgraded to the class \code{x} derives from.
  In that case, \code{padAndClip} does not act as an endomorphism.

  For \code{stackStrings}: Same as \code{padAndClip}. In addition it is
  guaranteed to have a rectangular shape i.e. to be a constant-width
  \link{XStringSet} object.
}

\author{H. Pages}

\seealso{
  \itemize{
    \item The \code{\link[GenomicAlignments]{stackStringsFromBam}} function
          in the \pkg{GenomicAlignments} package for stacking the read
          sequences (or their quality strings) stored in a BAM file on
          a region of interest.

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
x <- BStringSet(c(seq1="ABCD", seq2="abcdefghijk", seq3="", seq4="XYZ"))

padAndClip(x, IRanges(3, 8:5), Lpadding.letter=">", Rpadding.letter="<")
padAndClip(x, IRanges(1:-2, 7), Lpadding.letter=">", Rpadding.letter="<")

stackStrings(x, 2, 8)

stackStrings(x, -2, 8, shift=c(0, -11, 6, 7),
             Lpadding.letter="#", Rpadding.letter=".")

stackStrings(x, -2, 8, shift=c(0, -14, 6, 7),
             Lpadding.letter="#", Rpadding.letter=".")

stackStrings(x, -2, 8, shift=c(0, -14, 6, 7),
             Lpadding.letter="#", Rpadding.letter=".",
             remove.out.of.view.strings=TRUE)

library(hgu95av2probe)
probes <- DNAStringSet(hgu95av2probe)
probes

stackStrings(probes, 0, 26,
             Lpadding.letter="+", Rpadding.letter="-")

options(showHeadLines=15)
stackStrings(probes, 3, 23, shift=6*c(1:5, -(1:5)),
             Lpadding.letter="+", Rpadding.letter="N",
             remove.out.of.view.strings=TRUE)
}

\keyword{methods}
\keyword{manip}
