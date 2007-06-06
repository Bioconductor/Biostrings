\name{BStringViews-class}
\docType{class}

% Classes:
\alias{class:BStringViews}
\alias{BStringViews-class}

% Accesor methods:
\alias{subject}
\alias{subject,BStringViews-method}
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
\alias{length,BStringViews-method}
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
  \link{DNAString}, or \link{RNAString}) object.
}

\details{
    [TODO: Put some details here]
}

\section{Accesor methods}{
  In the code snippets below,
  \code{x} is a BStringViews object.

  \describe{
    \item{}{
      \code{subject(x)}:
      [TODO: Document me]
    }
    \item{}{
      \code{start(x)}:
      [TODO: Document me]
    }
    \item{}{
      \code{first(x)}:
      deprecated. Use 'start' instead.
    }
    \item{}{
      \code{end(x)}:
      [TODO: Document me]
    }
    \item{}{
      \code{last(x)}:
      deprecated. Use 'end' instead.
    }
    \item{}{
      \code{width(x)}:
      [TODO: Document me]
    }
    \item{}{
      \code{desc(x)}:
      [TODO: Document me]
    }
  }
}

\section{Standard generic methods}{
  In the code snippets below,
  \code{x}, \code{object}, \code{e1} and \code{e2} are BStringViews objects,
  and \code{i} can be a numeric or logical vector.

  \describe{
    \item{}{
      \code{length(x)}:
      Used to get the length of a BStringViews, i.e., its number of views
    }
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
  \link{BStringViews-constructors}, \link{BString}, \link{DNAString}, \link{RNAString}
}

\keyword{methods}
\keyword{classes}
