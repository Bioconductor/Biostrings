\name{MaskCollection-class}
\docType{class}

% Classes
\alias{class:MaskCollection}
\alias{MaskCollection-class}
\alias{MaskCollection}

% Basic accessor methods
\alias{nirlist}
\alias{nirlist,MaskCollection-method}
\alias{length,MaskCollection-method}
\alias{width,MaskCollection-method}
\alias{names,MaskCollection-method}

% Other accessor-like methods
\alias{isEmpty,MaskCollection-method}
\alias{max,MaskCollection-method}
\alias{min,MaskCollection-method}
\alias{maskedwidth,MaskCollection-method}
\alias{maskedratio,MaskCollection-method}

% "show" method and related
\alias{reduce,MaskCollection-method}
\alias{show,MaskCollection-method}

% Subsetting
\alias{[[,MaskCollection-method}
\alias{[[<-,MaskCollection-method}
\alias{[,MaskCollection-method}
\alias{[<-,MaskCollection-method}


\title{MaskCollection objects}

\description{
  The MaskCollection class is a container for storing a collection of masks
  that can be used to mask regions in a sequence.
}

\details{
  [TODO]
}

\section{Accesor methods}{
  In the code snippets below, \code{x} is a MaskCollection object.

  \describe{
    \item{}{
      \code{nirlist(x)}:
      [TODO]
    }
    \item{}{
      \code{length(x)}:
      [TODO]
    }
    \item{}{
      \code{width(x)}:
      [TODO]
    }
    \item{}{
      \code{names(x)}:
      [TODO]
    }
  }
}

\section{Other methods}{
  In the code snippets below, \code{x} is a MaskCollection object.

  \describe{
    \item{}{
      \code{isEmpty(x)}:
      [TODO]
    }
    \item{}{
      \code{max(x)}:
      [TODO]
    }
    \item{}{
      \code{min(x)}:
      [TODO]
    }
    \item{}{
      \code{maskedwidth(x)}:
      [TODO]
    }
    \item{}{
      \code{maskedratio(x)}:
      [TODO]
    }
    \item{}{
      \code{reduce(x)}:
      [TODO]
    }
  }
}

\section{Subsetting}{
  In the code snippets below, \code{x} is a MaskCollection object.

  \describe{
    \item{}{
      \code{x[[i]]}:
      [TODO]
    }
    \item{}{
      \code{x[i]}:
      [TODO]
    }
  }
}

\author{H. Pages}

\seealso{
  \link{MaskedXString-class},
  \link{NormalIRanges-class}
}

\examples{
  x1 <- new("NormalIRanges", start=c(11, 25, 28), width=c(5, 2, 2))
  x2 <- new("NormalIRanges", start=c(3, 10, 27), width=c(5, 8, 1))
  x3 <- new("NormalIRanges", start=7, width=2)
  x <- new("MaskCollection", nirlist=list(x1, x2, x3), width=29L)
  x
  length(x)
  width(x)
  reduce(x)
  x[[2]]
  length(x[[2]])
  x[[2]][-3]
}

\keyword{methods}
\keyword{classes}
