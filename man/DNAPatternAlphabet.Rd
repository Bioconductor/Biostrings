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
  substitution of 'U' for 'T'.

  \tabular{ll}{
    'N'\tabmatches 'A', 'C', 'G', 'T'\cr
    'B'\tabmatches 'C', 'G', 'T'\cr
    'D'\tabmatches 'A', 'G', 'T'\cr
    'H'\tabmatches 'A', 'C', 'T'\cr
    'K'\tabmatches 'G', 'T'\cr
    'M'\tabmatches 'A', 'C'\cr
    'R'\tabmatches 'A', 'G'\cr
    'S'\tabmatches 'C', 'G'\cr
    'V'\tabmatches 'A', 'C', 'G'\cr
    'W'\tabmatches 'A', 'T'\cr
    'Y'\tabmatches 'C', 'T'\cr
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
