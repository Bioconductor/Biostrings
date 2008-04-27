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
  In the context of the Biostrings package, a mask is a set of regions
  in a sequence that need to be excluded from some computation.
  For example, when applying \code{\link{alphabetFrequency}}
  or \code{\link{matchPattern}} on a chromosome sequence, you might want
  to exclude some regions like the centromere or the repeat regions.
  This can be achieved by putting one or several masks on the sequence
  before \code{\link{alphabetFrequency}} is called on it.

  A MaskCollection object is a vector-like object that represents
  such set of masks.
  Like standard R vectors, it has a "length" which is the number of
  masks contained in it, and a "width" which determines the length
  of the sequence it can be "put on". For example, a MaskCollection
  object of width 20000 can only be put on an XString object of 20000
  letters.
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
  \code{\link{alphabetFrequency}},
  \code{\link{matchPattern}},
  \link{NormalIRanges-class},
  \link{IRanges-utils}
}

\examples{
  ### Making a MaskCollection object:
  nir1 <- new("NormalIRanges", start=c(11, 25, 28), width=c(5, 2, 2))
  nir2 <- new("NormalIRanges", start=c(3, 10, 27), width=c(5, 8, 1))
  nir3 <- new("NormalIRanges", start=7, width=2)
  mymasks <- new("MaskCollection", nirlist=list(nir1, nir2, nir3), width=29L)
  mymasks
  length(mymasks)
  width(mymasks)

  ### Putting a MaskCollection object on a sequence:
  x <- DNAString("ACTCAAACTAGGATAGNACTTAACCAGGC")
  x
  length(x)  # same as width(mymasks)
  nchar(x)   # same as length(x)
  masks(x) <- mymasks
  x
  length(x)  # has not changed
  nchar(x)   # has changed
  alphabetFrequency(x)

  ### Removing the masks:
  masks(x) <- NULL
  x
  alphabetFrequency(x)

  ### Some operations on the MaskCollection object itself:
  reduce(mymasks)
  mymasks[[2]]
  length(mymasks[[2]])
  mymasks[[2]][-3]
  append(mymasks[-2], gaps(mymasks[2]))
  mymasks2 <- restrict(mymasks, start=8)
  mymasks2
  mymasks2[[2]]
}

\keyword{methods}
\keyword{classes}
