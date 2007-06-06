\name{BStringViews-constructors}

% Constructor-like functions and generics:
\alias{BStringViews-constructors}
\alias{initialize,BStringViews-method}
\alias{views}
\alias{adjacentViews}
\alias{mask}
\alias{mask,BString-method}
\alias{mask,BStringViews-method}
\alias{BStringViews}
\alias{BStringViews,ANY-method}
\alias{BStringViews,file-method}
\alias{BStringViews,BString-method}
\alias{BStringViews,BStringViews-method}


\title{Constructor-like functions and generics for BStringViews objects}

\description{
  [TODO: Put a description here]
}

\usage{
  views(subject, start=NA, end=NA)
  adjacentViews(subject, width, gapwidth=0)
  mask(x, ...)
  BStringViews(src, subjectClass, sep="")
}

\arguments{
  \item{subject}{
    [TODO]
  }
  \item{start}{
    [TODO]
  }
  \item{end}{
    [TODO]
  }
  \item{width}{
    [TODO]
  }
  \item{gapwidth}{
    [TODO]
  }
  \item{x}{
    [TODO]
  }
  \item{...}{
    additional arguments to be passed to or from methods.
  }
  \item{src}{
    [TODO]
  }
  \item{subjectClass}{
    [TODO]
  }
  \item{sep}{
    [TODO]
  }
}

\details{
  In the code snippets below,
  \code{subject} is a character vector or a \link{BString}
  (or \link{DNAString} or \link{RNAString}) object
  and \code{start}, \code{end} and \code{width} are numeric vectors.

  \describe{
    \item{}{
      \code{views(subject, start, end)}:
      [TODO: Document me]
    }
    \item{}{
      \code{adjacentViews(subject, width, gapwidth)}:
      [TODO: Document me]
    }
  }

  In the code snippets below,
  \code{x} can be a \link{BString} (or \link{DNAString}
  or \link{RNAString}) object or a BStringViews object,
  and \code{start} and \code{end} are numeric vectors.

  \describe{
    \item{}{
      \code{mask(x, start, end)}:
      [TODO: Document me]
    }
    \item{}{
      \code{mask(x)}:
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

\seealso{
  \link{BStringViews-class}, \link{BString}, \link{DNAString}, \link{RNAString}
}

\examples{
  ## A BStringViews object with 12 views
  v12 <- views(DNAString("TAATAATG"), -2:9, 0:11)
}

\keyword{methods}

