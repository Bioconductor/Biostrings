\name{IntIntervals-class}
\docType{class}

% Classes:
\alias{class:IntIntervals}
\alias{IntIntervals-class}
\alias{IntIntervals}

% Constructor-like functions and generics:
\alias{initialize,IntIntervals-method}

% Accesor methods:
\alias{length,IntIntervals-method}
\alias{start,IntIntervals-method}
\alias{nchar,IntIntervals-method}
\alias{end,IntIntervals-method}
\alias{names,IntIntervals-method}
\alias{names<-,IntIntervals-method}
\alias{width}
\alias{desc}
\alias{desc,ANY-method}
\alias{desc<-}
\alias{desc<-,ANY-method}

% Other methods:
\alias{show,IntIntervals-method}
\alias{[,IntIntervals-method}
\alias{[<-,IntIntervals-method}

\alias{as.data.frame,IntIntervals-method}
\alias{as.matrix,IntIntervals-method}

% Deprecated methods:
\alias{first}
\alias{first,Views-method}
\alias{last}
\alias{last,Views-method}


\title{IntIntervals objects}

\description{
  The IntIntervals class is a simple container for storing a set of integer
  intervals.
}

\details{
  An IntIntervals object is a data frame-like object where each row describes
  an interval (or range) of integers. Each interval is defined by its start and
  its length. The start can be any integer (see \code{start} below) but the
  length must be a nonnegative integer (see \code{nchar} below).

  An important difference with a standard R data frame object is that an
  IntIntervals object can only be subsetted by row whereas a standard R data
  frame object can be subsetted by row and by column. Also, the length of an
  IntIntervals object is its number of rows whereas the length of a standard R
  data frame object is its number of columns.

  Note that it is unlikely that the user will have to create or manipulate
  directly an IntIntervals instance when using the Biostrings package. However
  the IntIntervals class being the superclass of the \code{\link{BStringSet}}
  class, any \code{\link{BStringSet}} object is also an IntIntervals object and
  can be manipulated as such. Therefore all the methods described below also
  work with a \code{\link{BStringSet}} object.
}

\section{Accesor methods}{
  In the code snippets below, \code{x} is an IntIntervals object.

  \describe{
    \item{}{
      \code{length(x)}:
      The number of intervals in \code{x}.
    }
    \item{}{
      \code{start(x)}:
      The start values of the intervals.
      This is a vector of \code{length(x)} positive integers.
    }
    \item{}{
      \code{nchar(x)}:
      The lengths of the intervals (not to be confunded with the length of
      \code{x}). This is a vector of \code{length(x)} nonnegative integers.
    }
    \item{}{
      \code{width(x)}:
      \code{width} is an alias for \code{nchar}.
    }
    \item{}{
      \code{end(x)}:
      \code{start(x) + nchar(x) - 1L}
    }
    \item{}{
      \code{names(x)}:
      A character vector or \code{NULL} of the same length as \code{x}.
      The names are the only data in an IntIntervals object that can safely
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
  In the code snippets below, \code{x} is an IntIntervals objects,
  and \code{i} can be a numeric or logical vector.

  \describe{
    \item{}{
      \code{x[i]}:
      Return a new IntIntervals object made of the selected intervals.
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
      \code{start(x)} and \code{nchar(x)}.
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
  \link{BStringSet-class}
}

\examples{
  inters <- new("IntIntervals", start=4:1, nchar=0:3)
  inters
  start(inters)
  nchar(inters)
  end(inters)
  names(inters)
  names(inters) <- c("inter1", "inter2")
  inters

  ## Make a new IntIntervals object by selecting some intervals only:
  inters[4:2]                    # 3 intervals
  inters[-1]                     # 3 intervals
  inters[c(2,4)]                 # 2 intervals
  inters[nchar(inters) != 0]     # 3 intervals
  inters[is.na(names(inters))]   # 2 intervals
  inters[!is.na(names(inters))]  # 2 intervals
  inters[3]                      # 1 intervals
  inters[FALSE]                  # 0 intervals
}

\keyword{methods}
\keyword{classes}
