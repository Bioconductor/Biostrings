\name{letter}

\alias{letter}
\alias{letter,character-method}
\alias{letter,XString-method}
\alias{letter,XStringViews-method}
\alias{letter,MaskedXString-method}


\title{Subsetting a string}

\description{
  Extract a substring from a string by picking up individual
  letters by their position.
}

\usage{
  letter(x, i)
}

\arguments{
  \item{x}{
    A character vector, or an \link{XString}, \link{XStringViews}
    or \link{MaskedXString} object.
  }
  \item{i}{
    An integer vector with no NAs.
  }
}

\details{
  Unlike with the \code{substr} or \code{substring} functions,
  \code{i} must contain valid positions.
}

\value{
  A character vector of length 1 when \code{x} is an \link{XString}
  or \link{MaskedXString} object (the masks are ignored for the latter).

  A character vector of the same length as \code{x} when \code{x}
  is a character vector or an \link{XStringViews} object.

  Note that, because \code{i} must contain valid positions,
  all non-NA elements in the result are guaranteed to have exactly
  \code{length(i)} characters.
}

\seealso{
  \code{\link[IRanges]{subseq}},
  \link{XString-class},
  \link{XStringViews-class},
  \link{MaskedXString-class}
}

\examples{
  x <- c("abcd", "ABC")
  i <- c(3, 1, 1, 2, 1)

  ## With a character vector:
  letter(x[1], 3:1)
  letter(x, 3)
  letter(x, i)
  #letter(x, 4)             # Error!

  ## With a BString object:
  letter(BString(x[1]), i)  # returns a character vector
  BString(x[1])[i]          # returns a BString object

  ## With an XStringViews object:
  x2 <- XStringViews(x, "BString")
  letter(x2, i)
}

\keyword{methods}
