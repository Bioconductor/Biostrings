\name{gapletter<-}
\alias{gapletter<-}
\title{ function to modify the gap letter in an alphabet }
\description{
  This function modifies the gap letter in a "BioAlphabet" object.
}
\usage{
gapletter(x) <- value
}
\arguments{
  \item{x}{ An object of class "BioAlphabet" }
  \item{value}{ A single letter. }
}
\value{
  The modified object of class "BioAlphabet"
}
\author{ Saikat DebRoy }
\seealso{\code{\link{BioAlphabet-class}} }
\examples{
d <- DNAAlphabet()
gapletter(d) <- 'N'
d
r <- RNAPatternAlphabet()
gapletter(r) <- '*'
r
}
\keyword{classes}
