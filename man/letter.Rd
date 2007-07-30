\name{letter}

\alias{letter}
\alias{letter,character-method}
\alias{letter,BString-method}
\alias{letter,BStringViews-method}


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
    A character vector, a \link{BString} (or derived) object
    or a \link{BStringViews} object.
  }
  \item{i}{
    An NA-free numeric vector.
  }
}

\details{
  Unlike with the \code{\link{substr}} or \code{\link{substring}} functions,
  \code{i} must contain valid positions.
}

\value{
  A character vector of length 1 when \code{x} is
  a \link{BString} (or derived) object.

  A character vector of the same length as \code{x} when \code{x}
  is a character vector or a \link{BStringViews} object.

  Note that, because \code{i} must contain valid positions,
  all non-NA elements in the result are guaranted to have exactly
  \code{length(i)} characters.
}

\seealso{
  \code{\link{subBString}},
  \link{BString-class},
  \link{BStringViews-class}
}

\examples{
  letter("abcd", 3)

  x <- c("abcd", "ABC")
  i <- c(3, 1, 1, 2, 1)

  ## With a character vector
  letter(x, i)

  ## With a BString object
  letter(BString(x), i)

  ## With a BStringViews object
  x2 <- BStringViews(x, "BString")
  letter(x2, i)

  ## Error
  #letter(x, 4)
}

\keyword{methods}
