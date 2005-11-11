\name{LongestCommonPrefix}
\alias{LongestCommonPrefix}
\title{Find the longest common prefixes for a suffix array}
\description{
  Given an object of class "BioString" representing a suffix array, this
  function returns an object (of class "LongestCommonPrefix")
  representing a matrix of the longest common prefixes between any two
  elements in the suffix array.
}
\usage{
LongestCommonPrefix(x)
}
\arguments{
  \item{x}{ An object of class "BioString" representing a suffix array. }
}
\details{
  The suffix array argument to the \code{LongestCommonPrefix} function
  is assumed to be created by the \code{DNASuffixArray} function. In
  other words, the suffixes are assumed to be sorted in increasing
  order.
}
\value{
  An object of class "LongestCommonPrefix".
}
\author{Saikat DebRoy}
\seealso{\code{\link{DNASuffixArray}}, \code{\link{LongestCommonPrefix-class}}}
\examples{
data('yeastSEQCHR1')
yeast1 <- DNAString(yeastSEQCHR1)
x <- substring(yeast1, 1, 30)
x
suf <- DNASuffixArray(x)
lcp <- LongestCommonPrefix(suf)
}
\keyword{internal}
