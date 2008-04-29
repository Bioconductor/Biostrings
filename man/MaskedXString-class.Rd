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
  ## -------------------------------------------------------------------------
  ## A. MASKING BY POSITION
  ## -------------------------------------------------------------------------

  mask0 <- Mask(mask.width=29, start=c(3, 10, 25), width=c(6, 8, 5))
  x <- DNAString("ACTCAAACTAGGATAGNACTTNNCCAGGC")
  masks(x) <- mask0
  x
  gaps(x)

  ## -------------------------------------------------------------------------
  ## B. MASKING BY CONTENT
  ## -------------------------------------------------------------------------

  masks(x) <- "N"
  x
  as(x, "XStringViews")
  gaps(x)
  as(gaps(x), "XStringViews")

  library(BSgenome.Dmelanogaster.FlyBase.r51)
  chrU <- Dmelanogaster[["U"]]
  chrU
  alphabetFrequency(chrU)
  masks(chrU) <- "N"
  chrU
  alphabetFrequency(chrU)
  as(chrU, "XStringViews")
  as(gaps(chrU), "XStringViews")

  ## -------------------------------------------------------------------------
  ## C. MASKING WITH MORE THAN ONE MASK
  ## -------------------------------------------------------------------------

  mask2 <- Mask(mask.width=length(chrU), start=c(50000, 350000, 543900), width=25000)
  names(mask2) <- "some ugly regions"
  masks(chrU) <- append(masks(chrU), mask2)
  chrU
  as(chrU, "XStringViews")
  as(gaps(chrU), "XStringViews")
}

\keyword{methods}
\keyword{classes}
