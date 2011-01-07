\name{gregexpr2}
\alias{gregexpr2}

\title{A replacement for R standard gregexpr function}

\description{
This is a replacement for the standard gregexpr function that
does exact matching only.
Standard gregexpr() misses matches when they are overlapping.
The gregexpr2 function finds all matches but it only works in "fixed"
mode i.e. for exact matching (regular expressions are not
supported).
}

\usage{
gregexpr2(pattern, text)
}

\arguments{
  \item{pattern}{
    character string to be matched in the given character vector
  }
  \item{text}{
    a character vector where matches are sought
  }
}

\value{
  A list of the same length as \code{text} each element of
  which is an integer vector as in \code{gregexpr},
  except that the starting positions of all (even overlapping)
  matches are given.
  Note that, unlike \code{gregexpr}, \code{gregexpr2} doesn't attach a "match.length"
  attribute to each element of the returned list because, since it only works
  in "fixed" mode, then all the matches have the length of the pattern.
  Another difference with \code{gregexpr} is that with \code{gregexpr2},
  the \code{pattern} argument must be a single (non-NA, non-empty) string.
}

\author{H. Pages}

\seealso{\code{\link{gregexpr}}, \code{\link{matchPattern}}}

\examples{
    gregexpr("aa", c("XaaaYaa", "a"), fixed=TRUE)
    gregexpr2("aa", c("XaaaYaa", "a"))
}

\keyword{manip}
