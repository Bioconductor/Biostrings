\name{reverseComplement}
\alias{reverseComplement}
\title{ Function to reverse a DNA or RNA sequence and complement each base}
\description{
  Given an object of class "BioString", this function complements each
  base in the underlying string (replacing A with T or U, T or U with A,
  C with G and G with C) and reverses it. It also changes the substring
  indices so that they reflect the reversing.
}
\usage{
reverseComplement(x)
}
\arguments{
  \item{x}{ an object of class "BioString" }
}
\value{
  An object of class "BioString" containg the reversed and complemented
  underlying string and substring indices.
}
\author{ Saikat DebRoy }
\examples{
reverseComplement("ATCG-AA")
}
\keyword{manip}
