\name{longestConsecutive}
\alias{longestConsecutive}
\title{Obtain the length of the longest substring containing
  only 'letter'}
\description{
  This function accepts a character vector and computes the length of the
  longest substring containing only \code{letter} for each element of \code{x}.
}
\usage{
longestConsecutive(seq, letter) 
}
\arguments{
  \item{seq}{Character vector.}
  \item{letter}{Character vector of length 1, containing one single character.}
}
\details{
  The elements of \code{x} can be in upper case, lower case
  or mixed. NAs are handled.
}

\value{
  An integer vector of the same length as \code{x}.
}

\author{W. Huber}

\examples{
 v = c("AAACTGTGFG", "GGGAATT", "CCAAAAAAAAAATT")
 longestConsecutive(v, "A")
}
\keyword{manip}
