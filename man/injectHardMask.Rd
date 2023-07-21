\name{injectHardMask}

\alias{injectHardMask}
\alias{injectHardMask,XStringViews-method}
\alias{injectHardMask,MaskedXString-method}


\title{Injecting a hard mask in a sequence}

\description{
  \code{injectHardMask} allows the user to "fill" the masked regions
  of a sequence with an arbitrary letter (typically the \code{"+"}
  letter).
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

\details{
  The name of the \code{injectHardMask} function was chosen because of the
  primary use that it is intended for: converting a pile of active "soft
  masks" into a "hard mask".
  Here the pile of active "soft masks" refers to the active masks that have
  been put on top of a sequence. In Biostrings, the original sequence and the
  masks defined on top of it are bundled together in one of the dedicated
  containers for this: the \link{MaskedBString}, \link{MaskedDNAString},
  \link{MaskedRNAString} and \link{MaskedAAString} containers (this is the
  \link{MaskedXString} family of containers).
  The original sequence is always stored unmodified in a \link{MaskedXString}
  object so no information is lost. This allows the user to activate/deactivate
  masks without having to worry about losing the letters that are in the
  regions that are masked/unmasked. Also this allows better memory
  management since the original sequence never needs to be copied, even when
  the set of active/inactive masks changes.

  However, there are situations where the user might want to \emph{really}
  get rid of the letters that are in some particular regions by replacing
  them with a junk letter (e.g. \code{"+"}) that is guaranteed to not interfer
  with the analysis that s/he is currently doing.
  For example, it's very likely that a set of motifs or short reads will not
  contain the \code{"+"} letter (this could easily be checked) so they will
  never hit the regions filled with \code{"+"}.
  In a way, it's like the regions filled with \code{"+"} were masked but we
  call this kind of masking "hard masking".

  Some important differences between "soft" and "hard" masking:

  \itemize{
    \item \code{injectHardMask} creates a (modified) copy of the original
          sequence. Using "soft masking" does not.

    \item A function that is "mask aware" like \code{alphabetFrequency} or
          \code{matchPattern} will really skip the masked regions
          when "soft masking" is used i.e. they will not walk thru the
          regions that are under active masks. This might lead to some
          speed improvements when a high percentage of the original sequence
          is masked.
          With "hard masking", the entire sequence is walked thru.

    \item Matches cannot span over masked regions with "soft masking".
          With "hard masking" they can.
  }
}

\value{
  An \link{XString} object of the same length as the orignal object \code{x}
  if \code{x} is a \link{MaskedXString} object, or of the same length
  as \code{subject(x)} if it's an \link{XStringViews} object.
}

\author{H. Pagès}

\seealso{
  \code{\link{maskMotif}},
  \link{MaskedXString-class},
  \code{\link{replaceLetterAt}},
  \code{\link{chartr}},
  \link{XString},
  \link{XStringViews-class}
}

\examples{
## ---------------------------------------------------------------------
## A. WITH AN XStringViews OBJECT
## ---------------------------------------------------------------------
v2 <- Views("abCDefgHIJK", start=c(8, 3), end=c(14, 4))
injectHardMask(v2)
injectHardMask(v2, letter="=")

## ---------------------------------------------------------------------
## B. WITH A MaskedXString OBJECT
## ---------------------------------------------------------------------
mask0 <- Mask(mask.width=29, start=c(3, 10, 25), width=c(6, 8, 5))
x <- DNAString("ACACAACTAGATAGNACTNNGAGAGACGC")
masks(x) <- mask0
x
subject <- injectHardMask(x)

## Matches can span over masked regions with "hard masking":
matchPattern("ACggggggA", subject, max.mismatch=6)
## but not with "soft masking":
matchPattern("ACggggggA", x, max.mismatch=6)
}

\keyword{utilities}
\keyword{manip}
