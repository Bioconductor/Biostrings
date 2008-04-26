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

% Coercion:
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

% Other methods:
\alias{nchar,MaskedXString-method}
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
\alias{masks<-,MaskedXString,XString-method}
\alias{masks<-,MaskedXString,character-method}
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
  \link{XString-class},
  \link{MaskCollection-class},
  \link{XStringViews-class},
  \link{IRanges-utils}
}

\examples{
  ## masking "by position"
  #mask("AxyxyxBC", 2, 6)

  ## masking "by content"
  #mask("AxyxyxBC", "xyx")

  ## masking the N's in a chromosome sequence
  #library(BSgenome.Dmelanogaster.FlyBase.r51)
  #chrU <- Dmelanogaster[["U"]]
  #alphabetFrequency(chrU)
  #noN_chrU <- mask(chrU, "N")
  #alphabetFrequency(noN_chrU)

  ## masking a BStringViews object
  #mask(noN_chrU)
  #v <- views("AxyxyxBC", c(6, 4,-1, NA), c(7, 6, 1, 1))
  #mask(v)

  ## note that applying mask() again returns a BStringViews object
  ## where the original views on v have been merged, sorted from left to right
  ## and cleaned from their "out of limits" regions
  #mask(mask(v))

  ## finally, note that mask(mask(mask(v))) is _always_ the same as mask(v)
  #mask(mask(mask(v)))
}

\keyword{methods}
\keyword{classes}
