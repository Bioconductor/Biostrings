\name{DNASuffixArray}
\alias{DNASuffixArray}
\title{ function to create a suffix array for a DNA string }
\description{
  This function creates a suffix array for the given DNA strings. The
  sorting of the suffix array can be done based on a given number of
  letters at the begining of each suffix.
}
\usage{
DNASuffixArray(x, prefixLength = max(nchar(x)))
}
\arguments{
  \item{x}{ an objecct of class "BioString" representing one or more DNA
    strings. If its length is more than one, suffix arrays for each of
    the string are created and merged. }
  \item{prefixLength}{ the number of letters at the begining of each
    suffix that are used in sorting the suffix array. }
}
\value{
  An object of class "BioString", with the suffixes from the argument
  sorted in increasing order. For the purpose of sorting, each string is
  assumed to be padded at the end with a letter lower in order than
  any letter in the alphabet of the string.
}
\references{  }
\author{ Saikat DebRoy }
\seealso{ \code{\link{sortDNAString}} for sorting the strings in an
  object of class "BioString". }
\examples{
data('yeastSEQCHR1')
yeast1 <- DNAString(yeastSEQCHR1)
DNASuffixArray(substr(yeast1, 1, 30))
x <- substring(yeast1, c(1, 16), c(15, 30))
x
DNASuffixArray(x)
}
\keyword{ utilities }
