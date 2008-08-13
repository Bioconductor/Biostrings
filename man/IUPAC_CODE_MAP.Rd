\name{IUPAC_CODE_MAP}

\alias{IUPAC_CODE_MAP}
\alias{mergeIUPACLetters}


\title{The IUPAC Extended Genetic Alphabet}

\description{
  IUPAC nucleotide ambiguity codes are used for representing
  sequences of nucleotides where 

  [TODO]
}

\usage{
    IUPAC_CODE_MAP
    mergeIUPACLetters(x)
}

\arguments{
  \item{x}{
    A vector of non-empty character strings made of IUPAC letters.
  }
}

\details{
  [TODO]
}

\author{H. Pages}

\references{
  IUPAC-IUB SYMBOLS FOR NUCLEOTIDE NOMENCLATURE:
  Cornish-Bowden (1985) \emph{Nucl. Acids Res.} 13: 3021-3030.
}

\seealso{
  \code{\link{DNAString}}, \code{\link{RNAString}}
}

\examples{
  IUPAC_CODE_MAP
  mergeIUPACLetters(c("A", "AC", "MAAMC", "GM", "AB", "MK"))
}

\keyword{utilities}
\keyword{data}
