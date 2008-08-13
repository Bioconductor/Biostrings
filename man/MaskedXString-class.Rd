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

% Accesor methods:
\alias{unmasked}
\alias{unmasked,MaskedXString-method}
\alias{masks}
\alias{masks,XString-method}
\alias{masks,MaskedXString-method}
\alias{alphabet,MaskedXString-method}
\alias{length,MaskedXString-method}

% "maskedwidth" and related methods:
\alias{maskedwidth,MaskedXString-method}
\alias{maskedratio,MaskedXString-method}
\alias{nchar,MaskedXString-method}

% Coercion:
\alias{coerce,BString,MaskedBString-method}
\alias{coerce,DNAString,MaskedDNAString-method}
\alias{coerce,RNAString,MaskedRNAString-method}
\alias{coerce,AAString,MaskedAAString-method}

\alias{coerce,MaskedBString,BString-method}
\alias{coerce,MaskedDNAString,DNAString-method}
\alias{coerce,MaskedRNAString,RNAString-method}
\alias{coerce,MaskedAAString,AAString-method}
\alias{XString,MaskedXString-method}

\alias{coerce,MaskedXString,MaskCollection-method}
\alias{coerce,MaskedXString,NormalIRanges-method}
\alias{coerce,MaskedXString,XStringViews-method}

% Other methods:
\alias{reduce,MaskedXString-method}
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


\title{MaskedXString objects}

\description{
  The MaskedXString class is a container for storing a masked sequence.
}

\details{
  [TODO]
}

\section{Accesor methods}{
  In the code snippets below, \code{x} is a MaskedXString object.

  \describe{
    \item{}{
      \code{unmasked(x)}:
      Turns \code{x} into an \link{XString} object by dropping the masks.
    }
    \item{}{
      \code{masks(x)}:
      Turns \code{x} into a \link[IRanges]{MaskCollection} object by dropping the sequence.
    }
    \item{}{
      \code{alphabet(x)}:
      Equivalent to \code{alphabet(unmasked(x))}.
      See \code{?\link{alphabet}} for more information.
    }
    \item{}{
      \code{length(x)}:
      Equivalent to \code{length(unmasked(x))}.
      See \code{?`\link[Biostrings]{length},XString-method`} for more information.
    }
  }
}

\section{"maskedwidth" and related methods}{
  In the code snippets below, \code{x} is a MaskedXString object.

  \describe{
    \item{}{
      \code{maskedwidth(x)}:
      Get the number of masked letters in \code{x}. A letter is considered
      masked iff it's masked by at least one active mask.
    }
    \item{}{
      \code{maskedratio(x)}:
      Equivalent to \code{maskedwidth(x) / length(x)}.
    }
    \item{}{
      \code{nchar(x)}:
      Equivalent to \code{length(x) - maskedwidth(x)}.
    }
  }
}

\section{Coercion}{
  TODO
}

\section{Other methods}{
  TODO
}

\author{H. Pages}

\seealso{
  \code{\link{maskMotif}},
  \code{\link{injectHardMask}},
  \code{\link{alphabetFrequency}},
  \code{\link{reverse}},
  \link{XString-class},
  \link[IRanges]{MaskCollection-class},
  \link{XStringViews-class},
  \link[IRanges]{IRanges-utils}
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

  ## Prepare a MaskCollection object of 3 masks by running the examples in
  ## the man page for these objects:
  example(MaskCollection, package="IRanges")

  ## Put it on 'x':
  masks(x) <- mymasks
  x
  alphabetFrequency(x)

  ## Deactivate all masks:
  active(masks(x)) <- FALSE
  x

  ## Activate mask 3:
  active(masks(x))[3] <- FALSE
  x

  ## Turn MaskedXString object into an XStringViews object:
  as(x, "XStringViews")

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
