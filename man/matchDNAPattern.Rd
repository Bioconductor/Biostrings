\name{matchDNAPattern}
\alias{matchDNAPattern}
\alias{matchDNAPattern,character,ANY-method}
\alias{matchDNAPattern,ANY,character-method}
\title{ Generic to find all matches of a pattern in a DNA string}
\description{
  Generic that finds all matches of a pattern in a DNA string. Currently
  two algorithms are implemented. The default algorithm is an extension
  of the Boyer-Moore algorithm. The extended algorithm allows
  some wildcards in addition to the symbols for the bases and gap. The
  other algorithm is a simple forward search that examines all
  substrings of the full string of the same length as the pattern from
  the begining to end.
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
  \item{algorithm}{ Currently the only valid values are
    \code{"boyer-moore"} (which is the default) and
    \code{"forward-search"}. The forward search algorithm is often as
    fast as the more sphisticated Boyer-Moore algorithm when the
    patterns being matched are very simple.}
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
m1 <- matchDNAPattern("GCNNNAT", x)
m1
as.matrix(m1)
m2 <- matchDNAPattern("GCNNNAT", x, algorithm="forward-search")
m2
as.matrix(m2)
}
\keyword{classes}
