\name{SeqLocs-class}
\docType{class}

% Classes:
\alias{class:SeqLocs}
\alias{SeqLocs-class}
\alias{SeqLocs}

% Constructor-like functions and generics:
\alias{initialize,SeqLocs-method}

% Accesor methods:
\alias{length,SeqLocs-method}
\alias{start,SeqLocs-method}
\alias{nchar,SeqLocs-method}
\alias{end,SeqLocs-method}
\alias{names,SeqLocs-method}
\alias{names<-,SeqLocs-method}

% Other methods:

\alias{show,SeqLocs-method}
\alias{[,SeqLocs-method}
\alias{[<-,SeqLocs-method}

\alias{as.data.frame,SeqLocs-method}
\alias{as.matrix,SeqLocs-method}


\title{SeqLocs objects}

\description{
  The SeqLocs class is the basic container for storing a set of locations
  defined by their start/nchar relatively to an unspecified sequence.
}

\details{
  A SeqLocs object is a data frame-like object where each row describes a
  location in an unspecified sequence. Each location is defined by a starting
  position (a positive integer, see \code{start} below) and a number of
  letters (a nonnegative integer, see \code{nchar} below).

  An important difference with a standard R data frame object is that a SeqLocs
  object can only be subsetted by row whereas a standard R data frame object
  can be subsetted by row and by column. Also, the length of a SeqLocs object
  is its number of rows whereas the length of a standard R data frame object
  is its number of columns.

  Note that it is unlikely that the user will have to create or manipulate
  directly a SeqLocs instance when using the Biostrings package. However the
  SeqLocs class being the superclass of the \code{\link{BStringSet}} class,
  any \code{\link{BStringSet}} object is also a SeqLocs object and can be
  manipulated as such. Therefore all the methods described below also work
  with a \code{\link{BStringSet}} object.
}

\section{Accesor methods}{
  In the code snippets below, \code{x} is a SeqLocs object.

  \describe{
    \item{}{
      \code{length(x)}:
      The number of sequence locations in \code{x}.
    }
    \item{}{
      \code{start(x)}:
      The starting position of each location in the (unspecified) sequence.
      This is a vector of \code{length(x)} positive integers.
    }
    \item{}{
      \code{nchar(x)}:
      The number of letters of each location in the (unspecified) sequence.
      This is a vector of \code{length(x)} nonnegative integers.
    }
    \item{}{
      \code{end(x)}:
      \code{start(x) + nchar(x) - 1L}
    }
    \item{}{
      \code{names(x)}:
      A character vector or \code{NULL} of the same length as \code{x}.
      The names are the only data in a SeqLocs object that can safely
      be changed by the user. All other data are immutable!
      As a general recommendation, the user should never try to modify
      an object by accessing its slots directly.
    }
  }
}

\section{Other methods}{
  In the code snippets below, \code{x} is a SeqLocs objects,
  and \code{i} can be a numeric or logical vector.

  \describe{
    \item{}{
      \code{x[i]}:
      Return a new SeqLocs object made of the selected locations.
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

\author{H. Pages}

\seealso{
  \link{BStringSet-class}
}

\examples{
  locs <- new("SeqLocs", start=4:1, nchar=0:3)
  locs
  start(locs)
  nchar(locs)
  end(locs)
  names(locs)
  names(locs) <- c("loc1", "loc2")
  locs

  ## Make a new SeqLocs object by selecting some locations only:
  locs[4:2]                  # 3 locations
  locs[-1]                   # 3 locations
  locs[c(2,4)]               # 2 locations
  locs[nchar(locs) != 0]     # 3 locations
  locs[is.na(names(locs))]   # 2 locations
  locs[!is.na(names(locs))]  # 2 locations
  locs[3]                    # 1 locations
  locs[FALSE]                # 0 locations
}

\keyword{methods}
\keyword{classes}
