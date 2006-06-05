\name{BStringViews}
\docType{class}

% Classes:
\alias{class:BStringViews}
\alias{BStringViews-class}

% Accesor methods:
\alias{subject}
\alias{subject,BStringViews-method}
\alias{first}
\alias{first,BStringViews-method}
\alias{last}
\alias{last,BStringViews-method}
\alias{width}
\alias{width,BStringViews-method}
\alias{desc}
\alias{desc,BStringViews-method}

% Constructor-like functions and generics:
\alias{initialize,BStringViews-method}
\alias{views}
\alias{adjacentViews}
\alias{BStringViews}
\alias{BStringViews,file-method}
\alias{BStringViews,BString-method}
\alias{BStringViews,BStringViews-method}

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
\alias{!=,BStringViews,BString-method}
\alias{!=,BString,BStringViews-method}
\alias{!=,BStringViews,ANY-method}
\alias{!=,ANY,BStringViews-method}
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
      \code{first(x)}:
      [TODO: Document me]
    }
    \item{}{
      \code{last(x)}:
      [TODO: Document me]
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

\section{Constructor-like functions and generics}{
  In the code snippets below,
  \code{subject} is a character vector or a \link{BString}
  (or \link{DNAString} or \link{RNAString}) object
  and \code{first}, \code{last} and \code{width} are numeric vectors.

  \describe{
    \item{}{
      \code{views(subject, first, last)}:
      [TODO: Document me]
    }
    \item{}{
      \code{adjacentViews(subject, width, gapwidth)}:
      [TODO: Document me]
    }
  }

  In the code snippets below,
  \code{src} can be a character vector, a connection object,
  a \link{BString} (or \link{DNAString} or \link{RNAString}) object
  or a BStringViews object,
  \code{subjectClass} is a character-string with value \code{"BString"},
  \code{"DNAString"} or \code{"RNAString"}
  and \code{sep} is a character-string.

  \describe{
    \item{}{
      \code{BStringViews(src, subjectClass, sep)}:
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
  \link{BString}, \link{DNAString}, \link{RNAString}
}

\examples{
  ## A BStringViews object with 12 views
  v12 <- views(DNAString("TAATAATG"), -2:9, 0:11)
}

\keyword{methods}
\keyword{classes}
