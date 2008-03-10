\name{IRanges-class}
\docType{class}

% Classes
\alias{class:IRanges}
\alias{IRanges-class}
\alias{IRanges}

% Accesor methods
\alias{length,IRanges-method}
\alias{start,IRanges-method}
\alias{width}
\alias{width,IRanges-method}
\alias{end,IRanges-method}
\alias{names,IRanges-method}
\alias{desc}
\alias{desc,ANY-method}

% Initialization
\alias{initialize,IRanges-method}

% Replacement methods
\alias{start<-}
\alias{start<-,IRanges-method}
\alias{width<-}
\alias{width<-,IRanges-method}
\alias{end<-}
\alias{end<-,IRanges-method}
\alias{names<-,IRanges-method}
\alias{desc<-}
\alias{desc<-,ANY-method}
\alias{update,IRanges-method}

% Show method
\alias{as.data.frame,IRanges-method}
\alias{show,IRanges-method}

% Subsetting
\alias{[,IRanges-method}
\alias{[<-,IRanges-method}

% Other methods
\alias{as.matrix,IRanges-method}

% Deprecated methods
\alias{first}
\alias{first,IRanges-method}
\alias{last}
\alias{last,IRanges-method}


\title{IRanges objects}

\description{
  The IRanges class is a simple container for storing a set of integer
  ranges.
}

\details{
  An IRanges object is a data frame-like object where each row describes
  a range of integers. Each range is defined by its start and
  its length. The start can be any integer (see \code{start} below) but the
  length must be a nonnegative integer (see \code{width} below).

  An important difference with a standard R data frame object is that an
  IRanges object can only be subsetted by row whereas a standard R data
  frame object can be subsetted by row and by column. Also, the length of an
  IRanges object is its number of rows whereas the length of a standard R
  data frame object is its number of columns.

  Note that it is unlikely that the user will have to create or manipulate
  directly an IRanges instance when using the Biostrings package. However
  the IRanges class being the superclass of the \link{BStringViews}
  class, any \link{BStringViews} object is also an IRanges object and
  can be manipulated as such. Therefore all the methods described here also
  work with a \link{BStringViews} object.
}

\section{Accesor methods}{
  In the code snippets below, \code{x} is an IRanges object.

  \describe{
    \item{}{
      \code{length(x)}:
      The number of ranges in \code{x}.
    }
    \item{}{
      \code{start(x)}:
      The start values of the ranges.
      This is a vector of \code{length(x)} positive integers.
    }
    \item{}{
      \code{width(x)}:
      The number of integers of the ranges.
      This is a vector of \code{length(x)} nonnegative integers.
    }
    \item{}{
      \code{end(x)}:
      \code{start(x) + width(x) - 1L}
    }
    \item{}{
      \code{names(x)}:
      A character vector or \code{NULL} of the same length as \code{x}.
      The names are the only data in an IRanges object that can safely
      be changed by the user. All other data are immutable!
      As a general recommendation, the user should never try to modify
      an object by accessing its slots directly.
    }
    \item{}{
      \code{desc(x)}:
      \code{desc} is an alias for \code{names}.
    }
  }
}

\section{Other methods}{
  In the code snippets below, \code{x} is an IRanges objects,
  and \code{i} can be a numeric or logical vector.

  \describe{
    \item{}{
      \code{x[i]}:
      Return a new IRanges object made of the selected ranges.
      \code{i} can be a numeric vector, a logical vector, \code{NULL}
      or missing.
    }
    \item{}{
      \code{as.data.frame(x, row.names=NULL, optional=FALSE, ...)}:
      Converts \code{x} into a standard R data frame object.
    }
    \item{}{
      \code{as.matrix(x, ...)}:
      Converts \code{x} into a 2-column integer matrix containing
      \code{start(x)} and \code{width(x)}.
    }
  }
}

\section{Deprecated methods}{
    \item{}{
      \code{first(x)}:
      deprecated. Use 'start' instead.
    }
    \item{}{
      \code{last(x)}:
      deprecated. Use 'end' instead.
    }
}

\author{H. Pages}

\seealso{
  \link{BStringViews-class}
}

\examples{
  x <- new("IRanges", start=4:1, width=0:3)
  x
  start(x)
  width(x)
  end(x)
  names(x)
  names(x) <- c("range1", "range2")
  x

  ## Make a new IRanges object by selecting some ranges only:
  x[4:2]               # 3 ranges
  x[-1]                # 3 ranges
  x[c(2,4)]            # 2 ranges
  x[width(x) != 0]     # 3 ranges
  x[is.na(names(x))]   # 2 ranges
  x[!is.na(names(x))]  # 2 ranges
  x[3]                 # 1 ranges
  x[FALSE]             # 0 ranges
}

\keyword{methods}
\keyword{classes}
