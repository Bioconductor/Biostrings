\name{NucleotideString}
\alias{NucleotideString}
\title{Function to convert a character string to a BioString object }
\description{
  Given a character string representing either an RNA or DNA sequence,
  convert it to a BioString object with appropriate alphabet.
}
\usage{
NucleotideString(src, type = c("DNA", "RNA"),
                 srctype = c("character", "connection"))
}
\arguments{
  \item{src}{ A character string }
  \item{type}{ Either "DNA" (the default) or "RNA". }
  \item{srctype}{ Currently the only valid value (and the default) is
    "character". In future, we may allow "connection" as another valid
    value denoting that the \code{src} argument is a file name.}
}
\value{
  An object of class "BioString".
}
\author{Saikat DebRoy}
\seealso{\code{\link{BioString-class}}}
\examples{
NucleotideString("ACTGAACT")
}
\keyword{classes}
