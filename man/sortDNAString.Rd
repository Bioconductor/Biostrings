\name{sortDNAString}
\alias{sortDNAString}
\title{ function to sort DNA strings }
\description{
  This function sorts a given object of class "BioString" representing
  a vector of DNA strings in increasing order. The sorting can be done
  based only on a given number of letters at the begining of each string.
}
\usage{
sortDNAString(x, prefixLength = max(nchar(x)))
}
\arguments{
  \item{x}{ an objecct of class "BioString" representing one or more DNA
    strings. }
  \item{prefixLength}{ the number of letters at the begining of each
    string that are used in the sorting. }
}
\value{
  An object of class "BioString", with the strings from the argument
  sorted in increasing order. For the purpose of sorting, each string is
  assumed to be padded at the end with a letter lower in order than
  any letter in the alphabet of the string.
}
\author{ Saikat DebRoy }
\seealso{ \code{\link{DNASuffixArray}} }
\examples{
data('yeastSEQCHR1')
yeast1 <- DNAString(yeastSEQCHR1)
x <- substring(yeast1, seq(1, by=10, length=30),
                    seq(10, by=10, length=30))
x
sortDNAString(x)
x <- substr(yeast1, 1, 30)
x
sortDNAString(substring(x, 1:30, 30)) # suffix array for x
}
\keyword{ utilities }
