\name{DNAAlphabet}
\alias{DNAAlphabet}
\alias{RNAAlphabet}
\title{ functions to create DNA and RNA alphabets }
\description{
  The function \code{DNAAlphabet} returns an alphabet consisting of the
  letters \sQuote{A}, \sQuote{C}, \sQuote{G}, \sQuote{T} and \sQuote{-} (the last representing a gap).

  The function \code{RNAAlphabet} returns an alphabet consisting of the
  letters \sQuote{A}, \sQuote{C}, \sQuote{G}, \sQuote{U} and \sQuote{-} (the last representing a gap).
}
\usage{
DNAAlphabet()
RNAAlphabet()
}
\value{
  An object of class "NucleotideAlphabet". See
  \code{\link{NucleotideAlphabet-class}} for details.
}
\author{ Saikat DebRoy }
\seealso{ \code{\link{DNAPatternAlphabet}} and
  \code{\link{RNAPatternAlphabet}} return alphabets for standrad
  patterns used in DNA and RNA sequences. }
\examples{
DNAAlphabet()
RNAAlphabet()
}
\keyword{classes}
