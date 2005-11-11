\name{LongestCommonPrefix-class}
\docType{class}
\alias{LongestCommonPrefix-class}
\alias{[,LongestCommonPrefix,matrix,missing,missing-method}
\title{Class of longest common prefix lengths for suffix array elemets}
\description{This class represents a matrix of the length of the longest
common prefixes between elements of a suffix array. The diagonal
elements are all zero.}
\section{Objects from the Class}{
  Objects can be created by calling the function
  \code{\link{LongestCommonPrefix}} with the (only) argument a suffix
  array.
}
\section{Slots}{
  \describe{
    \item{\code{abovediag}:}{Object of class \code{"integer"},
      represents the elements just above diagonal (ie. the elements of
      the matrix whose column numbers are one more than there row numbers). }
  }
}
\section{Methods}{
  \describe{
    \item{lcp[i]}{The only method for extracting elements from this
      class of objects. Here \code{lcp} is an object of class
      "LongestCommonPrefix" and \code{i} is a matrix of two columns, the
      columns representing the rwos and columns of the elements to be
      extracted. The result is an integer vector of the corresponding
      elements.}
  }
}
\references{Dan Gusfield (1997) \emph{Algorithms on strings, trees, and
    sequences}, Cambridge University Press, pp. 152}
\author{Saikat DebRoy}
\seealso{
  \code{\link{LongestCommonPrefix}}
}
\examples{
data('yeastSEQCHR1')
yeast1 <- DNAString(yeastSEQCHR1)
x <- substring(yeast1, 1, 30) # very different from c(1,30)
x
suf <- DNASuffixArray(x)
lcp <- LongestCommonPrefix(suf)
lcp[cbind(1, 1:30)]
}
\keyword{classes}
