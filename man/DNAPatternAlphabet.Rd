\name{DNAPatternAlphabet}
\alias{DNAPatternAlphabet}
\alias{RNAPatternAlphabet}
\title{ functions to create DNA and RNA pattern alphabets }
\description{
  The function \code{DNAAlphabet} returns an alphabet consisting of
  standard patterns used in DNA sequences.

  The function \code{RNAAlphabet} returns an alphabet consisting of
  standard patterns used in RNA sequences.
}
\usage{
DNAPatternAlphabet()
RNAPatternAlphabet()
}
\section{Memebers of the pattern alphabets}{
  In the table below, we provide the letters from the base DNA
  alphabet that are matched by the added letters in the DNA pattern
  alphabet. The corresponding table for RNA is identical except for the
  substitution of \sQuote{U} for \sQuote{T}.

  \tabular{ll}{
    \sQuote{N}\tabmatches \sQuote{A}, \sQuote{C}, \sQuote{G}, \sQuote{T}\cr
    \sQuote{B}\tabmatches \sQuote{C}, \sQuote{G}, \sQuote{T}\cr
    \sQuote{D}\tabmatches \sQuote{A}, \sQuote{G}, \sQuote{T}\cr
    \sQuote{H}\tabmatches \sQuote{A}, \sQuote{C}, \sQuote{T}\cr
    \sQuote{K}\tabmatches \sQuote{G}, \sQuote{T}\cr
    \sQuote{M}\tabmatches \sQuote{A}, \sQuote{C}\cr
    \sQuote{R}\tabmatches \sQuote{A}, \sQuote{G}\cr
    \sQuote{S}\tabmatches \sQuote{C}, \sQuote{G}\cr
    \sQuote{V}\tabmatches \sQuote{A}, \sQuote{C}, \sQuote{G}\cr
    \sQuote{W}\tabmatches \sQuote{A}, \sQuote{T}\cr
    \sQuote{Y}\tabmatches \sQuote{C}, \sQuote{T}\cr
    \tab \cr
  }
}
\value{
  An object of class "BioPatternAlphabet". See
  \code{\link{BioPatternAlphabet-class}} for details.
}
\author{ Saikat DebRoy }
\seealso{ \code{\link{DNAAlphabet}}, \code{\link{RNAAlphabet}}}
\examples{
DNAPatternAlphabet()
RNAPatternAlphabet()
}
\keyword{classes}
