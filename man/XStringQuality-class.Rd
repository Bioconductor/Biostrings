\name{XStringQuality-class}
\docType{class}

% Classes:
\alias{class:XStringQuality}
\alias{XStringQuality-class}
\alias{XStringQuality}

\alias{class:PhredQuality}
\alias{PhredQuality-class}
\alias{PhredQuality}
\alias{class:SolexaQuality}
\alias{SolexaQuality-class}
\alias{SolexaQuality}

% Coercion:
\alias{coerce,character,PhredQuality-method}
\alias{coerce,BString,PhredQuality-method}
\alias{coerce,BStringSet,PhredQuality-method}
\alias{coerce,integer,PhredQuality-method}
\alias{coerce,numeric,PhredQuality-method}
\alias{coerce,PhredQuality,integer-method}
\alias{as.integer,PhredQuality-method}
\alias{coerce,PhredQuality,numeric-method}
\alias{as.numeric,PhredQuality-method}

\alias{coerce,character,SolexaQuality-method}
\alias{coerce,BString,SolexaQuality-method}
\alias{coerce,BStringSet,SolexaQuality-method}
\alias{coerce,integer,SolexaQuality-method}
\alias{coerce,numeric,SolexaQuality-method}
\alias{coerce,SolexaQuality,integer-method}
\alias{as.integer,SolexaQuality-method}
\alias{coerce,SolexaQuality,numeric-method}
\alias{as.numeric,SolexaQuality-method}

\alias{as.matrix,XStringQuality-method}
\alias{coerce,XStringQuality,matrix-method}


\title{PhredQuality and SolexaQuality objects}

\description{
  Objects for storing string quality measures.
}

\usage{
## Constructors:
PhredQuality(x)
SolexaQuality(x)
}

\arguments{
  \item{x}{
    Either a character vector, \link{BString}, \link{BStringSet}, integer vector,
    or number vector of error probabilities.
  }
}

\details{
  \code{PhredQuality} objects store characters that are interpreted as [0 - 99] quality
  measures by subtracting 33 from their ASCII decimal representation
  (e.g. ! = 0, " = 1, \# = 2, ...).

  \code{SolexaQuality} objects store characters are interpreted as [-5 - 99] quality
  measures by subtracting 64 from their ASCII decimal representation
  (e.g. ; = -5, < = -4, = = -3, ...).
}

\author{P. Aboyoun}

\seealso{
  \code{\link{pairwiseAlignment}},
  \link{PairwiseAlignedXStringSet-class},
  \link{DNAString-class},
  \link{BStringSet-class}
}

\examples{
  PhredQuality(0:40)
  SolexaQuality(0:40)

  PhredQuality(seq(1e-4,0.5,length=10))
  SolexaQuality(seq(1e-4,0.5,length=10))

  x <- SolexaQuality(BStringSet(c(a="@ABC", b="abcd")))
  as.matrix(x)
}

\keyword{methods}
\keyword{classes}
