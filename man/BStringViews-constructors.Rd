\name{BStringViews-constructors}

% Constructor-like functions and generics:
\alias{BStringViews-constructors}
\alias{initialize,BStringViews-method}
\alias{views}
\alias{adjacentViews}
\alias{BStringViews}
\alias{BStringViews,file-method}
\alias{BStringViews,BString-method}
\alias{BStringViews,BStringViews-method}


\title{Constructor-like functions and generics for BStringViews objects}

\description{
    [TODO: Put a description here]
}

\usage{
  views(subject, first=NA, last=NA)
  adjacentViews(subject, width, gapwidth=0)
  BStringViews(src, subjectClass, sep="")
}

\arguments{
  \item{subject}{
    [TODO]
  }
  \item{first}{
    [TODO]
  }
  \item{last}{
    [TODO]
  }
  \item{width}{
    [TODO]
  }
  \item{gapwidth}{
    [TODO]
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

\author{H. Pages}

\seealso{
  \link{BStringViews-class}, \link{BString}, \link{DNAString}, \link{RNAString}
}

\examples{
  ## A BStringViews object with 12 views
  v12 <- views(DNAString("TAATAATG"), -2:9, 0:11)
}

\keyword{methods}
