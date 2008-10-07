\name{QualityScaledXStringSet-class}
\docType{class}

% Classes:
\alias{class:QualityScaledXStringSet}
\alias{QualityScaledXStringSet-class}
\alias{QualityScaledXStringSet}

\alias{class:QualityScaledBStringSet}
\alias{QualityScaledBStringSet-class}
\alias{QualityScaledBStringSet}
\alias{class:QualityScaledDNAStringSet}
\alias{QualityScaledDNAStringSet-class}
\alias{QualityScaledDNAStringSet}
\alias{class:QualityScaledRNAStringSet}
\alias{QualityScaledRNAStringSet-class}
\alias{QualityScaledRNAStringSet}
\alias{class:QualityScaledAAStringSet}
\alias{QualityScaledAAStringSet-class}
\alias{QualityScaledAAStringSet}

% Accesor methods:
\alias{quality}
\alias{quality,QualityScaledXStringSet-method}

% Standard generic methods:
\alias{show,QualityScaledXStringSet-method}
\alias{[,QualityScaledXStringSet,ANY,ANY,ANY-method}


\title{QualityScaledBStringSet, QualityScaledDNAStringSet, QualityScaledRNAStringSet and QualityScaledAAStringSet objects}

\description{
  The QualityScaledBStringSet class is a container for storing a
  \code{\link{BStringSet}} object with an \code{\link{XStringQuality}} object.

  Similarly, the QualityScaledDNAStringSet (or QualityScaledRNAStringSet, or
  QualityScaledAAStringSet) class is a container for storing a \code{\link{DNAStringSet}}
  (or \code{\link{RNAStringSet}}, or \code{\link{AAStringSet}}) objects with
  an \code{\link{XStringQuality}} object.
}

\usage{
  ## Constructors:
  QualityScaledBStringSet(x, quality)
  QualityScaledDNAStringSet(x, quality)
  QualityScaledRNAStringSet(x, quality)
  QualityScaledAAStringSet(x, quality)
}

\arguments{
  \item{x}{
    Either a character vector, or an \link{XString}, \link{XStringSet}
    or \link{XStringViews} object.
  }
  \item{quality}{
    An \link{XStringQuality} object.
  }
}

\details{
  The \code{QualityScaledBStringSet}, \code{QualityScaledDNAStringSet},
  \code{QualityScaledRNAStringSet} and \code{QualityScaledAAStringSet}
  functions are constructors that can be used to "naturally" turn
  \code{x} into an QualityScaledXStringSet object of the desired subtype.
}

\section{Accesor methods}{
  The QualityScaledXStringSet class derives from the \link{XStringSet}
  class hence all the accessor methods defined for an \link{XStringSet}
  object can also be used on an QualityScaledXStringSet object. Common
  methods include (in the code snippets below, \code{x} is an
  QualityScaledXStringSet object):

  \describe{
    \item{}{
      \code{length(x)}:
      The number of sequences in \code{x}.
    }
    \item{}{
      \code{width(x)}:
      A vector of non-negative integers containing the number
      of letters for each element in \code{x}.
    }
    \item{}{
      \code{nchar(x)}:
      The same as \code{width(x)}.
    }
    \item{}{
      \code{names(x)}:
      \code{NULL} or a character vector of the same length as \code{x} containing
      a short user-provided description or comment for each element in \code{x}.
    }
    \item{}{
      \code{quality(x)}:
      The quality of the strings.
    }
  }
}

\section{Subsetting and appending}{
  In the code snippets below,
  \code{x} and \code{values} are XStringSet objects,
  and \code{i} should be an index specifying the elements to extract.

  \describe{
    \item{}{
      \code{x[i]}:
      Return a new QualityScaledXStringSet object made of the selected elements.
    }
  }
}

\author{P. Aboyoun}

\seealso{
  \link{BStringSet-class},
  \link{DNAStringSet-class},
  \link{RNAStringSet-class},
  \link{AAStringSet-class},
  \link{XStringQuality-class}
}

\examples{
  x1 <- DNAStringSet(c("TTGA", "CTCN"))
  q1 <- PhredQuality(c("*+,-", "6789"))
  qx1 <- QualityScaledDNAStringSet(x1, q1)
  qx1
}

\keyword{methods}
\keyword{classes}

