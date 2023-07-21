\name{InDel-class}
\docType{class}

% Classes
\alias{class:InDel}
\alias{InDel-class}
\alias{InDel}

% Accessor methods:
\alias{insertion}
\alias{insertion,InDel-method}
\alias{deletion}
\alias{deletion,InDel-method}


\title{InDel objects}

\description{
  The \code{InDel} class is a container for storing insertion and deletion
  information.
}

\details{
  This is a generic class that stores any insertion and deletion information.
}

\section{Accessor methods}{
  In the code snippets below,
  \code{x} is a \code{InDel} object.

  \describe{
    \item{\code{insertion(x)}:}{
      The insertion information.
    }
    \item{\code{deletion(x)}:}{
      The deletion information.
    }
  }
}

\author{P. Aboyoun}

\seealso{
  \code{\link{pairwiseAlignment}},
  \code{\link{PairwiseAlignments-class}}
}

\keyword{methods}
\keyword{classes}
