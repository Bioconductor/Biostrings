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

\author{H. Pages}

\seealso{
  \code{\link{maskMotif}},
  \code{\link{alphabetFrequency}},
  \code{\link{reverse}},
  \link{XString-class},
  \link{MaskCollection-class},
  \link{XStringViews-class},
  \link[IRanges]{IRanges-utils}
}

\examples{
  ## Masking by position
  mask0 <- Mask(mask.width=29, start=c(3, 10, 25), width=c(6, 8, 5))
  x <- DNAString("ACACAACTAGATAGNACTNNGAGAGACGC")
  masks(x) <- mask0
  x
  gaps(x)

  ## See ?maskMotif for masking by content
}

\keyword{methods}
\keyword{classes}
