\name{anySameLetter}
\alias{anySameLetter}
\alias{anySameLetter,character,ANY-method}
\title{ Generic to find the strings which contain a particular letter}
\description{
  This method finds all strings in a "BioString" object which contains
  at least one occurance of a particular letter. Note that this does not
  do a pattern matching. So, if the letter is "N" for a BioString object
  with DNA patterns, then it would only match strings which contain at
  least one instance of "N" and not others. See the examples for an
  illustration.
}
\usage{
anySameLetter(x, letter)
}
\arguments{
  \item{x}{ An object of class "BioString" or a character vector. In the
  later case it is converted to a BioString object.}
  \item{letter}{ A single letter. It should be part of the alphabet of
    \code{x}. }
}
\value{
  A logical vector of same length as \code{x} with elements \code{TRUE}
  or \code{FALSE} accordingly as the corresponding element of \code{x}
  contains \code{letter} or not.
}
\author{ Saikat DebRoy }
\seealso{\code{\link{BioString-class}}}
\examples{
PpiI <- "GAACNNNNNCTC"
D1 <-
  DNAString("tgctgatgcatagctagctgGAACtagctCTCtcgtagctggatgctgatNNNNNNNNNNNN")
matches <- matchDNAPattern(PpiI, D1)
matches
anyN <- anySameLetter(matches, 'N')
anyN
matches[!anyN]
}
\keyword{classes}
\keyword{methods}
