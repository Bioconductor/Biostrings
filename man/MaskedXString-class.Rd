\name{MaskedXString-class}
\docType{class}

% Classes:
\alias{class:MaskedXString}
\alias{MaskedXString-class}
\alias{MaskedXString}
\alias{class:MaskedBString}
\alias{MaskedBString-class}
\alias{MaskedBString}
\alias{class:MaskedDNAString}
\alias{MaskedDNAString-class}
\alias{MaskedDNAString}
\alias{class:MaskedRNAString}
\alias{MaskedRNAString-class}
\alias{MaskedRNAString}
\alias{class:MaskedAAString}
\alias{MaskedAAString-class}
\alias{MaskedAAString}

% Accessor-like methods:
\alias{unmasked}
\alias{unmasked,MaskedXString-method}
\alias{unmasked,XString-method}
\alias{masks}
\alias{masks,XString-method}
\alias{masks,MaskedXString-method}
\alias{length,MaskedXString-method}
\alias{maskedwidth,MaskedXString-method}
\alias{maskedratio,MaskedXString-method}
\alias{nchar,MaskedXString-method}

% Coercion:
\alias{seqtype,MaskedXString-method}
\alias{seqtype<-,MaskedXString-method}

\alias{coerce,MaskedXString,MaskedBString-method}
\alias{coerce,MaskedXString,MaskedDNAString-method}
\alias{coerce,MaskedXString,MaskedRNAString-method}
\alias{coerce,MaskedXString,MaskedAAString-method}

\alias{coerce,BString,MaskedBString-method}
\alias{coerce,DNAString,MaskedDNAString-method}
\alias{coerce,RNAString,MaskedRNAString-method}
\alias{coerce,AAString,MaskedAAString-method}

\alias{coerce,MaskedBString,BString-method}
\alias{coerce,MaskedDNAString,DNAString-method}
\alias{coerce,MaskedRNAString,RNAString-method}
\alias{coerce,MaskedAAString,AAString-method}

\alias{coerce,MaskedXString,MaskCollection-method}
\alias{coerce,MaskedXString,NormalIRanges-method}
\alias{coerce,MaskedXString,XStringViews-method}
\alias{coerce,MaskedXString,Views-method}

% Other methods:
\alias{collapse,MaskedXString-method}
\alias{gaps,MaskedXString-method}
\alias{subseq,MaskedXString-method}
\alias{as.character,MaskedXString-method}
\alias{toString,MaskedXString-method}
\alias{show,MaskedXString-method}

% "masks<-" methods:
\alias{masks<-}
\alias{masks<-,MaskedXString,NULL-method}
\alias{masks<-,MaskedXString,MaskCollection-method}
\alias{masks<-,XString,NULL-method}
\alias{masks<-,XString,ANY-method}

% "Views" method:
\alias{Views,MaskedXString-method}


\title{MaskedXString objects}

\description{
  The MaskedBString, MaskedDNAString, MaskedRNAString and MaskedAAString
  classes are containers for storing masked sequences.

  All those containers derive directly (and with no additional slots)
  from the MaskedXString virtual class.
}

\details{
  In Biostrings, a pile of masks can be put on top of a sequence.
  A pile of masks is represented by a
  \link[IRanges:MaskCollection-class]{MaskCollection} 
  object and the sequence by an \link{XString} object.
  A MaskedXString object is the result of bundling them together in a
  single object.

  Note that, no matter what masks are put on top of it, the original
  sequence is always stored unmodified in a MaskedXString object.
  This allows the user to activate/deactivate masks without having to
  worry about losing the information stored in the masked/unmasked
  regions.
  Also this allows efficient memory management since the original
  sequence never needs to be copied (modifying it would require to
  make a copy of it first - sequences cannot and should never be modified
  in place in Biostrings), even when the set of active/inactive masks
  changes.
}

\section{Accessor methods}{
  In the code snippets below, \code{x} is a MaskedXString object.
  For \code{masks(x)} and \code{masks(x) <- y}, it can also be
  an \link{XString} object and \code{y} must be \code{NULL} or
  a \link[IRanges:MaskCollection-class]{MaskCollection} object.

  \describe{
    \item{\code{unmasked(x)}:}{
      Turns \code{x} into an \link{XString} object by dropping the masks.
    }
    \item{\code{masks(x)}:}{
      Turns \code{x} into a
    \link[IRanges:MaskCollection-class]{MaskCollection} object by 
      dropping the sequence.
    }
    \item{\code{masks(x) <- y}:}{
      If \code{x} is an \link{XString} object and \code{y} is \code{NULL},
      then this doesn't do anything.

      If \code{x} is an \link{XString} object and \code{y} is a
      \link[IRanges:MaskCollection-class]{MaskCollection} object, then
      this turns \code{x} into a MaskedXString object by putting the
      masks in \code{y} on top of it.

      If \code{x} is a MaskedXString object and \code{y} is \code{NULL},
      then this is equivalent to \code{x <- unmasked(x)}.

      If \code{x} is a MaskedXString object and \code{y} is a
      \link[IRanges:MaskCollection-class]{MaskCollection} object, then
      this replaces the masks currently on top of \code{x} by the masks
      in \code{y}. 
    }
    \item{\code{alphabet(x)}:}{
      Equivalent to \code{alphabet(unmasked(x))}.
      See \code{?\link{alphabet}} for more information.
    }
    \item{\code{length(x)}:}{
      Equivalent to \code{length(unmasked(x))}.
      See
      \code{?`\link[Biostrings:XString-class]{length,XString-method}`}
      for more information.
    }
  }
}

\section{"maskedwidth" and related methods}{
  In the code snippets below, \code{x} is a MaskedXString object.

  \describe{
    \item{\code{maskedwidth(x)}:}{
      Get the number of masked letters in \code{x}. A letter is considered
      masked iff it's masked by at least one active mask.
    }
    \item{\code{maskedratio(x)}:}{
      Equivalent to \code{maskedwidth(x) / length(x)}.
    }
    \item{\code{nchar(x)}:}{
      Equivalent to \code{length(x) - maskedwidth(x)}.
    }
  }
}

\section{Coercion}{
  In the code snippets below, \code{x} is a MaskedXString object.

  \describe{
    \item{\code{as(x, "Views")}:}{
      Turns \code{x} into a \link[IRanges]{Views} object where the
      views are the unmasked regions of the original sequence
      ("unmasked" means not masked by at least one active mask).
    }
  }
}

\section{Other methods}{
  In the code snippets below, \code{x} is a MaskedXString object.

  \describe{
    \item{\code{collapse(x)}:}{
      Collapses the set of masks in \code{x} into a single mask made of
      all active masks.
    }
    \item{\code{gaps(x)}:}{
      Reverses all the masks i.e. each mask is replaced by a mask
      where previously unmasked regions are now masked and
      previously masked regions are now unmasked.
    }
  }
}

\author{H. Pagès}

\seealso{
  \itemize{
    \item \code{\link{maskMotif}}
    \item \code{\link{injectHardMask}}
    \item \code{\link{alphabetFrequency}}
    \item \code{\link{reverseComplement}}
    \item \link{XString-class}
    \item \link[IRanges]{MaskCollection-class}
    \item \link[IRanges]{Views-class}
  }
}

\examples{
## ---------------------------------------------------------------------
## A. MASKING BY POSITION
## ---------------------------------------------------------------------
mask0 <- Mask(mask.width=29, start=c(3, 10, 25), width=c(6, 8, 5))
x <- DNAString("ACACAACTAGATAGNACTNNGAGAGACGC")
length(x)  # same as width(mask0)
nchar(x)   # same as length(x)
masks(x) <- mask0
x
length(x)  # has not changed
nchar(x)   # has changed
gaps(x)

## Prepare a MaskCollection object of 3 masks ('mymasks') by running the
## examples in the man page for these objects:
example(MaskCollection, package="IRanges")

## Put it on 'x':
masks(x) <- mymasks
x
alphabetFrequency(x)

## Deactivate all masks:
active(masks(x)) <- FALSE
x

## Activate mask "C":
active(masks(x))["C"] <- TRUE
x

## Turn MaskedXString object into a Views object:
as(x, "Views")

## Drop the masks:
masks(x) <- NULL
x
alphabetFrequency(x)


## ---------------------------------------------------------------------
## B. MASKING BY CONTENT
## ---------------------------------------------------------------------
## See ?maskMotif for masking by content
}

\keyword{methods}
\keyword{classes}
