\name{matchDNAPattern}
\alias{matchDNAPattern}
\title{ Generic to find all matches of a pattern in a DNA string}
\description{
  Generic that finds all matches of a pattern in a DNA string using an
  extended Boyer-Morre algorithm. The extended algorithm allows
  some wildcards in addition to the symbols for the bases and gap.
}
\usage{
matchDNAPattern(pattern, x, algorithm)
}
\arguments{
  \item{pattern}{ An object representing the pattern string. The string in
    \code{pattern} can use the letters A,C,G,T,- (the last being the gap
    character) and also the wildcards N (matching A,C,G,T), V (matching
    A,G,C), R (matching A,G) and Y (matching C,T).
  }
  \item{x}{ An object representing a DNA string. }
  \item{algorithm}{ Currently the only valid value (which is also the
    default) is "boyer-moore" }
}
\value{
  An object of class "BioString" with the same length as the number of
  matches. Each element in the "BioString" object is a match. To obtain
  the start and end points of the matches, use \code{as.matrix} on the
  return value. See documentation for the "BioString" class for more
  details.
}
\author{ Saikat DebRoy }
\seealso{ \code{BioString-class} for the type of the return value. }
\examples{
x <- DNAString("AAGCGCGATATG")
m <- matchDNAPattern("GCNNNAT", x)
m
as.matrix(m)
}
\keyword{classes}
