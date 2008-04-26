\name{MaskCollection-class}
\docType{class}

% Classes:
\alias{class:MaskCollection}
\alias{MaskCollection-class}
\alias{MaskCollection}

% Basic accessor methods:
\alias{nirlist}
\alias{nirlist,MaskCollection-method}
\alias{length,MaskCollection-method}
\alias{width,MaskCollection-method}
\alias{names,MaskCollection-method}

% Other methods:
\alias{isEmpty,MaskCollection-method}
\alias{max,MaskCollection-method}
\alias{min,MaskCollection-method}
\alias{maskedwidth}
\alias{maskedwidth,MaskCollection-method}
\alias{maskedratio}
\alias{maskedratio,MaskCollection-method}

% Transforming methods:
\alias{restrict,MaskCollection-method}
\alias{reduce,MaskCollection-method}
\alias{gaps,MaskCollection-method}

% "show" method :
\alias{show,MaskCollection-method}

% Subsetting:
\alias{[[,MaskCollection-method}
\alias{[[<-,MaskCollection-method}
\alias{[,MaskCollection-method}
\alias{[<-,MaskCollection-method}

% Coercion:
\alias{coerce,MaskCollection,NormalIRanges-method}

% "append" method:
\alias{append,MaskCollection-method}


\title{MaskCollection objects}

\description{
  The MaskCollection class is a container for storing a collection of masks
  that can be used to mask regions in a sequence.
}

\details{
  [TODO]
}

\section{Basic accesor methods}{
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
  }
}

\section{Transforming methods}{
  In the code snippets below, \code{x} is a MaskCollection object.

  \describe{
    \item{}{
      \code{restrict(x, start=NA, end=NA, use.names=TRUE)}:
      [TODO]
    }
    \item{}{
      \code{reduce(x)}:
      [TODO]
    }
    \item{}{
      \code{gaps(x)}:
      [TODO]
    }
  }
}

\section{Subsetting and appending}{
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
    \item{}{
      \code{append(x, values, after=length(x))}:
      [TODO]
    }
  }
}

\author{H. Pages}

\seealso{
  \link{MaskedXString-class},
  \link{NormalIRanges-class},
  \link{IRanges-utils}
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
  append(x[-2], gaps(x[2]))
  y <- restrict(x, start=8)
  y
  y[[2]]
}

\keyword{methods}
\keyword{classes}
