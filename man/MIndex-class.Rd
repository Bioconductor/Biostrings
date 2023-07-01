\name{MIndex-class}
\docType{class}

% Classes
\alias{MIndex}
\alias{class:MIndex}
\alias{MIndex-class}
\alias{class:ByPos_MIndex}
\alias{ByPos_MIndex-class}

% Generics and methods:
\alias{length,MIndex-method}
\alias{width0}
\alias{width0,MIndex-method}
\alias{names,MIndex-method}
\alias{names<-,MIndex-method}
\alias{startIndex}
\alias{endIndex}
\alias{elementNROWS,MIndex-method}
\alias{unlist,MIndex-method}
\alias{coerce,MIndex,CompressedIRangesList-method}

\alias{show,MIndex-method}
\alias{[[,ByPos_MIndex-method}
\alias{startIndex,ByPos_MIndex-method}
\alias{endIndex,ByPos_MIndex-method}


\title{MIndex objects}

\description{
  The MIndex class is the basic container for storing the matches of
  a set of patterns in a subject sequence.
}

\details{
  An MIndex object contains the matches (start/end locations) of
  a set of patterns found in an \link{XString} object
  called "the subject string" or "the subject sequence" or
  simply "the subject".

   \code{\link{matchPDict}} function returns an MIndex object.
}

\section{Accessor methods}{
  In the code snippets below, \code{x} is an MIndex object.

  \describe{
    \item{\code{length(x)}:}{
      The number of patterns that matches are stored for.
    }
    \item{\code{names(x)}:}{
      The names of the patterns that matches are stored for.
    }
    \item{\code{startIndex(x)}:}{
      A list containing the starting positions of the matches
      for each pattern.
    }
    \item{\code{endIndex(x)}:}{
      A list containing the ending positions of the matches
      for each pattern.
    }
    \item{\code{elementNROWS(x)}:}{
      An integer vector containing the number of matches
      for each pattern.
    }
  }
}

\section{Subsetting methods}{
  In the code snippets below, \code{x} is an MIndex object.

  \describe{
    \item{\code{x[[i]]}:}{
      Extract the matches for the i-th pattern as an
    \link[IRanges:IRanges-class]{IRanges} object. 
    }
  }
}

\section{Coercion}{
  In the code snippets below, \code{x} is an MIndex object.

  \describe{
    \item{\code{as(x, "CompressedIRangesList")}:}{
      Turns \code{x} into an
    \link[IRanges]{CompressedIRangesList} object. 
      This coercion changes \code{x} from one
    \link[IRanges]{IntegerRangesList} 
      subtype to another with the underlying
    \link[IRanges]{IntegerRanges} values remaining unchanged.
    }
  }
}

\section{Other utility methods and functions}{
  In the code snippets below,
  \code{x} and \code{mindex} are MIndex objects
  and \code{subject} is the \link{XString} object
  containing the sequence in which the matches were found.

  \describe{
    \item{\code{unlist(x, recursive=TRUE, use.names=TRUE)}:}{
      Return all the matches in a single
    \link[IRanges:IRanges-class]{IRanges} object. 
      \code{recursive} and \code{use.names} are ignored.
    }
    \item{\code{extractAllMatches(subject, mindex)}:}{
      Return all the matches in a single \link{XStringViews} object.
    }
  }
}

\author{H. Pagès}

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
