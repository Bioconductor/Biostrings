\name{allSameLetter}
\alias{allSameLetter}
\alias{allSameLetter,character,ANY-method}
\title{ Generic to find the strings which are repeats of a single letter}
\description{
  This method finds all strings in a "BioString" object which consist of
  a particular letter repeating again and again. Note that this does not
  do a pattern matching. So, if the letter is "N" for a BioString object
  with DNA patterns, then it would only match strings which are repeats
  of "N" and not others. See the examples for an illustration.
}
\usage{
allSameLetter(x, letter)
}
\arguments{
  \item{x}{ An object of class "BioString" or a character vector. In the
  later case it is converted to a BioString object.}
  \item{letter}{ A single letter. It should be part of the alphabet of
    \code{x}. }
}
\value{
  A logical vector of same length as \code{x} with elements \code{TRUE}
  or \coode{FALSE} accordingly as the corresponding element of \code{x}
  is a repeat of \code{letter} or not.
}
\author{ Saikat DebRoy }
\seealso{\code{\link{BioString-class}}}
\examples{
PpiI <- "GAACNNNNNCTC"
D1 <-
  DNAString("tgctgatgcatagctagctgGAACtagctCTCtcgtagctggatgctgatNNNNNNNNNNNN")
matches <- matchDNAPattern(PpiI, D1)
matches
allN <- allSameLetter(matches, 'N')
allN
matches[!allN]
}
\keyword{classes}
\keyword{methods}
