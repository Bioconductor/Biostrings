\name{alphabetFrequency}
\alias{alphabetFrequency}
\title{ function to calculate the frequency of letters in a biological sequence }
\description{
  Given an object of class "BioString" representing a biological
  sequence, this function calculates the frequency of each letter in the
  (base) alphabet for the "BioString" object.
}
\usage{
alphabetFrequency(x, baseOnly = TRUE)
}
\arguments{
  \item{x}{ An object of class "BioString". }
  \item{baseOnly}{ A logical value. If \code{TRUE}, the frequency table
    only contains the letters in the base alphabet for \code{x}. }
}
\value{
  An integer vector with names same as the letters in the alphabet. If
  the base alphabet is used for the frequency counts, there is an extra
  element with name 'Others' and it is the number of letters in the
  string which are not in the base alphabet. Usually, a non-zero entry
  here means that the string contains some letters from a pattern
  alphabet.
}
\author{ Saikat DebRoy }
\seealso{ \code{\link{BioString-class}},
  \code{\link{BioAlphabet-class}}, \code{\link{BioPatternAlphabet-class}} }
\examples{
data(yeastSEQCHR1)
yeast1 <- DNAString(yeastSEQCHR1)
alphabetFrequency(yeast1)
alphabetFrequency(yeast1, baseOnly=FALSE)
}
\keyword{category}
