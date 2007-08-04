\name{BStringPartialMatches-class}
\docType{class}

% Classes:
\alias{class:BStringPartialMatches}
\alias{BStringPartialMatches-class}

% Accesor methods:
\alias{subpatterns}
\alias{subpatterns,BStringPartialMatches-method}
\alias{pattern}
\alias{pattern,BStringPartialMatches-method}

% Standard generic methods:
\alias{show,BStringPartialMatches-method}
\alias{[,BStringPartialMatches-method}

\title{The BStringPartialMatches class}

\description{
  WARNING: This class is currently under development and might
  not work properly! Full documentation will come later.

  Please send questions/comments to hpages@fhcrc.org

  Thanks for your comprehension!
}

\details{
  [TODO]
}

\section{Accesor methods}{
  In the code snippets below,
  \code{x} is a BStringPartialMatches object.

  \describe{
    \item{}{
      \code{subpatterns(x)}:
      [TODO]
    }
    \item{}{
      \code{pattern(x)}:
      [TODO]
    }
  }
}

\section{Standard generic methods}{
  In the code snippets below,
  \code{x} is a BStringPartialMatches objects,
  and \code{i} can be a numeric or logical vector.

  \describe{
    \item{}{
      \code{x[i]}:
      Return a new BStringPartialMatches object made of the selected views.
      \code{i} can be a numeric vector, a logical vector, \code{NULL} or missing.
      The returned object has the same subject as \code{x}.
    }
  }
}

\author{H. Pages}

\seealso{
  \link{BStringViews-class},
  \link{BString-class},
  \link{DNAString-class},
  \link{RNAString-class},
  \link{AAString-class},
  \code{\link{letter}}
}

\examples{
  ## [TODO]
}

\keyword{methods}
\keyword{classes}
