\name{BStringViews-constructors}

% Constructor-like functions and generics:
\alias{BStringViews-constructors}
\alias{views}
\alias{adjacentViews}
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
  BStringViews(src, subjectClass, collapse="")
}

\arguments{
  \item{subject}{
    A \link{BString} (or derived) object or a single string.
  }
  \item{start}{
    An integer vector containing the start positions of the views.
  }
  \item{end}{
    An integer vector containing the end positions of the views.
  }
  \item{width}{
    An integer vector containing the widths of the views.
  }
  \item{gapwidth}{
    An integer vector containing the widths of the gaps between the views.
  }
  \item{src}{
    A \link{BString} (or derived) object or a character vector.
  }
  \item{subjectClass}{
    [TODO]
  }
  \item{collapse}{
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
  \code{src} can be a character vector, a \link{BString} (or derived) object
  or a BStringViews object,
  \code{subjectClass} is a character-string with value \code{"BString"},
  \code{"DNAString"}, \code{"RNAString"} or \code{"AAString"},
  and \code{collapse} is a character string.

  \describe{
    \item{}{
      \code{BStringViews(src, subjectClass, collapse)}:
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

