\name{BioStringNewValues}
\alias{BioStringNewValues}
\title{generic to do create an object suitable for the values slot of
  "BioStrings" class}
\description{
  This generic creates an external vector object that contains (in its
  tag field) an R vector (of type either \code{INTSXP} or
  \code{CHARSXP}) of the given length. The actual type is determined
  according to the size of the alphabet.
}
\usage{
BioStringNewValues(alphabet, length.string)
}
\arguments{
  \item{alphabet}{ An object of class "BioAlphabet" }
  \item{length.string}{ An integer, the length of the storage in the result }
}
\value{
  An external pointer.
}
\author{Saikat DebRoy}
\seealso{\code{\link{BioString-class}} for the class which uses this
  generic during its initialization and \code{\link{BioAlphabet-class}}
  for the class of the \code{alphabet} parameter.}
\examples{
BioStringNewValues(DNAAlphabet(), 0)
}
\keyword{classes}
\keyword{internal}
