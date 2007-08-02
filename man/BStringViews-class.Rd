\name{BStringViews-class}
\docType{class}

% Classes:
\alias{class:BStringViews}
\alias{BStringViews-class}

% Accesor methods:
\alias{subject}
\alias{subject,BStringViews-method}
\alias{length,BStringViews-method}
\alias{start}
\alias{start,BStringViews-method}
\alias{first}
\alias{first,BStringViews-method}
\alias{end}
\alias{end,BStringViews-method}
\alias{last}
\alias{last,BStringViews-method}
\alias{width}
\alias{width,BStringViews-method}
\alias{desc}
\alias{desc,BStringViews-method}

% Standard generic methods:
\alias{show,BStringViews-method}
\alias{[,BStringViews-method}
\alias{[<-,BStringViews-method}
\alias{[[,BStringViews-method}
\alias{[[<-,BStringViews-method}
\alias{==,BStringViews,BString-method}
\alias{==,BString,BStringViews-method}
\alias{==,BStringViews,ANY-method}
\alias{==,ANY,BStringViews-method}
\alias{==,BStringViews,BStringViews-method}
\alias{!=,BStringViews,BString-method}
\alias{!=,BString,BStringViews-method}
\alias{!=,BStringViews,ANY-method}
\alias{!=,ANY,BStringViews-method}
\alias{!=,BStringViews,BStringViews-method}
\alias{as.character,BStringViews-method}
\alias{toString,BStringViews-method}
\alias{nchar,BStringViews-method}
\alias{as.matrix,BStringViews-method}
\alias{as.list,BStringViews-method}


\title{The BStringViews class}

\description{
  A container for a set of views on the same \link{BString} (or
  derived) object.
}

\details{
  A BStringViews object contains a set of views on the same
  \link{BString} (or derived) object called "the subject string"
  or "the subject sequence" or simply "the subject".
  Each view is defined by its start and end locations: both are
  integers such that start <= end.
  Note that two views can overlap and that a view can be "out of limits"
  i.e. it starts before the first letter of the subject or/and ends
  after its last letter. 
}

\section{Accesor methods}{
  In the code snippets below,
  \code{x} is a BStringViews object.

  \describe{
    \item{}{
      \code{subject(x)}:
      The subject of \code{x}. This is always a \link{BString}
      (or derived) object.
    }
    \item{}{
      \code{length(x)}:
      The number of views in \code{x}.
    }
    \item{}{
      \code{start(x)}:
      The starting positions of the views.
      An integer vector of the same length as \code{x}.
    }
    \item{}{
      \code{first(x)}:
      deprecated. Use 'start' instead.
    }
    \item{}{
      \code{end(x)}:
      The ending positions of the views.
      An integer vector of the same length as \code{x}.
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
      containing a short user provided description or comment for each view.
      These are the only data in a BStringViews object that can safely
      be changed by the user. All the other data are immutable!
      As a general recommendation, the user should never try to modify
      an object by accessing its slots directly.
    }
  }
}

\section{Standard generic methods}{
  In the code snippets below,
  \code{x}, \code{object}, \code{e1} and \code{e2} are BStringViews objects,
  and \code{i} can be a numeric or logical vector.

  \describe{
    \item{}{
      \code{x[i]}:
      [TODO: Document me]
    }
    \item{}{
      \code{x[[i]]}:
      [TODO: Document me]
    }
    \item{}{
      \code{e1 == e2}:
      [TODO: Document me]
    }
    \item{}{
      \code{e1 != e2}:
      [TODO: Document me]
    }
    \item{}{
      \code{as.character(x)}:
      [TODO: Document me]
    }
    \item{}{
      \code{toString(x)}:
      [TODO: Document me]
    }
    \item{}{
      \code{nchar(x)}:
      [TODO: Document me]
    }
    \item{}{
      \code{as.matrix(x)}:
      [TODO: Document me]
    }
    \item{}{
      \code{as.list(x)}:
      [TODO: Document me]
    }
  }
}

\author{H. Pages}

\seealso{
  \link{BStringViews-constructors},
  \link{BString-class},
  \link{DNAString-class},
  \link{RNAString-class},
  \link{AAString-class},
  \code{\link{letter}}
}

\examples{
  s <- DNAString("-CTC-N")

  ## One standard way to create a BStringViews object is to use
  ## the "views" constructor:
  v4 <- views(s, 3:0, 5:8)
  v4

  ## Attach a comment to view #4:
  desc(v4)[4] <- "out of limits"
}

\keyword{methods}
\keyword{classes}
