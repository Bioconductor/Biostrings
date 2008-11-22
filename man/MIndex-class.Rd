\name{MIndex-class}
\docType{class}

% Classes
\alias{MIndex}
\alias{class:MIndex}
\alias{MIndex-class}
\alias{class:ByPos_MIndex}
\alias{ByPos_MIndex-class}
\alias{class:ByName_MIndex}
\alias{ByName_MIndex-class}

% Generics and methods:
\alias{length,MIndex-method}
\alias{startIndex}
\alias{startIndex,MIndex-method}
\alias{endIndex}
\alias{countIndex}
\alias{countIndex,MIndex-method}
\alias{[[,MIndex-method}
\alias{$,MIndex-method}
\alias{unlist}
\alias{unlist,MIndex-method}

\alias{names,ByPos_MIndex-method}
\alias{show,ByPos_MIndex-method}
\alias{[[,ByPos_MIndex-method}
\alias{startIndex,ByPos_MIndex-method}
\alias{endIndex,ByPos_MIndex-method}
\alias{countIndex,ByPos_MIndex-method}

\alias{names,ByName_MIndex-method}
\alias{show,ByName_MIndex-method}
\alias{[[,ByName_MIndex-method}
\alias{startIndex,ByName_MIndex-method}
\alias{endIndex,ByName_MIndex-method}
\alias{countIndex,ByName_MIndex-method}


\title{MIndex objects}

\description{
  The MIndex class is the basic container for storing the matches of
  a set of patterns in a subject sequence.
}

\details{
  THIS IS STILL WORK IN PROGRESS!

  An MIndex object contains the matches (start/end locations) of
  a set of patterns found in an \link{XString} object
  called "the subject string" or "the subject sequence" or
  simply "the subject".

  The \code{\link{matchPDict}} function returns an MIndex object.

  MORE TO COME SOON...
}

\section{Accesor methods}{
  In the code snippets below,
  \code{x} is an MIndex object
  and \code{all.names} must be \code{TRUE} or \code{FALSE}.

  \describe{
    \item{}{
      \code{length(x)}:
      The number of patterns that matches are stored for.
    }
    \item{}{
      \code{names(x)}:
      The names of the patterns that matches are stored for.
    }
    \item{}{
      \code{startIndex(x, all.names=FALSE)}:
      A list containing the starting positions of the matches
      for each pattern.
    }
    \item{}{
      \code{endIndex(x, all.names=FALSE)}:
      A list containing the ending positions of the matches
      for each pattern.
    }
    \item{}{
      \code{countIndex(x, all.names=FALSE)}:
      An integer vector containing the number of matches
      for each pattern.
    }
  }
}

\section{Subsetting methods}{
  In the code snippets below,
  \code{x} is an MIndex object.

  \describe{
    \item{}{
      \code{x[[i]]}:
      Extract the matches for the i-th pattern as an \link[IRanges]{IRanges} object.
    }
  }
}

\section{Other utility methods and functions}{
  In the code snippets below,
  \code{x} and \code{mindex} are MIndex objects
  and \code{subject} is the \link{XString} object
  containing the sequence in which the matches were found.

  \describe{
    \item{}{
      \code{unlist(x, recursive=TRUE, use.names=TRUE)}:
      Return all the matches in a single \link[IRanges]{IRanges} object.
      \code{recursive} and \code{use.names} are ignored.
    }
    \item{}{
      \code{extractAllMatches(subject, mindex)}:
      Return all the matches in a single \link{XStringViews} object.
    }
  }
}

\author{H. Pages}

\seealso{
  \code{\link{matchPDict}},
  \link{PDict-class},
  \link[IRanges]{IRanges-class},
  \link{XStringViews-class}
}

\examples{
  ## See ?matchPDict and ?`matchPDict-inexact` for some examples.
}

\keyword{methods}
\keyword{classes}
