\name{Views-class}
\docType{class}

% Classes:
\alias{class:Views}
\alias{Views-class}
\alias{Views}

% Constructor-like functions and generics:
\alias{initialize,Views-method}

% Accesor methods:
\alias{length,Views-method}
\alias{nview}
\alias{nview,Views-method}
\alias{start,Views-method}
\alias{end,Views-method}
\alias{first}
\alias{first,Views-method}
\alias{last}
\alias{last,Views-method}
\alias{width}
\alias{width,Views-method}
\alias{desc,Views-method}
\alias{desc<-,Views-method}

% Other methods:

\alias{show,Views-method}
\alias{[,Views-method}
\alias{[<-,Views-method}

\alias{as.data.frame,Views-method}
\alias{as.matrix,Views-method}


\title{The Views class}

\description{
  The Views class is the basic container for storing a set of start/end
  locations.
}

\details{
  A Views object is a data frame-like object containing a set of start/end
  locations. Each row in the data frame describes a start/end location
  called "a view".
  An important difference with a standard R data frame object is that a Views
  object can only be subsetted by row whereas a standard R data frame object
  can be subsetted by row and by column. Also, the length of a Views object
  is its number of rows whereas the length of a standard R data frame object
  is its number of columns.
  In a Views object, each view is defined by its start and end locations: both
  are integers such that start <= end.
}

\section{Accesor methods}{
  In the code snippets below, \code{x} is a Views object.

  \describe{
    \item{}{
      \code{length(x)}:
      The number of views in \code{x}.
    }
    \item{}{
      \code{nview(x)}:
      Same as \code{length(x)}.
    }
    \item{}{
      \code{start(x)}:
      The starting positions of the views.
      An integer vector of the same length as \code{x}.
    }
    \item{}{
      \code{end(x)}:
      The ending positions of the views.
      An integer vector of the same length as \code{x}.
    }
    \item{}{
      \code{first(x)}:
      deprecated. Use 'start' instead.
    }
    \item{}{
      \code{last(x)}:
      deprecated. Use 'end' instead.
    }
    \item{}{
      \code{width(x)}:
      The width of the views.
      Equivalent to \code{end(x) - start(x) + 1L}.
      A vector of positive integers of the same length as \code{x}.
    }
    \item{}{
      \code{desc(x)}:
      \code{NULL} or a character vector of the same length as \code{x}
      containing a short user-provided description or comment for each view.
      These are the only data in a Views object that can safely
      be changed by the user. All the other data are immutable!
      As a general recommendation, the user should never try to modify
      an object by accessing its slots directly.
    }
  }
}

\section{Other methods}{
  In the code snippets below, \code{x} is a Views objects,
  and \code{i} can be a numeric or logical vector.

  \describe{
    \item{}{
      \code{x[i]}:
      Return a new Views object made of the selected views.
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
      \code{start(x)} and \code{end(x)}.
    }
  }
}

\author{H. Pages}

\seealso{
  \link{BStringViews-class}
}

\examples{
  v <- new("Views", start=3:0, end=5:8)
  v
  start(v)
  end(v)
  width(v)
  desc(v)
  desc(v) <- c("view1", "view2")
  v

  ## Make a new Views object by selecting some views only:
  v[4:2]                # 3 views
  v[-1]                 # 3 views
  v[c(2,4)]             # 2 views
  v[width(v) >= 5]      # 3 views
  v[is.na(desc(v))]     # 2 views
  v[!is.na(desc(v))]    # 2 views
  v[3]                  # 1 view
  v[FALSE]              # 0 view
}

\keyword{methods}
\keyword{classes}
