\name{IUPAC_CODE_MAP}

\alias{IUPAC_CODE_MAP}
\alias{lettersToAmbiguity}


\title{The IUPAC Extended Genetic Alphabet}

\description{
  [TODO]
}

\usage{
    IUPAC_CODE_MAP
    lettersToAmbiguity(x)
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

\seealso{
  \code{\link{DNAString}}, \code{\link{RNAString}}
}

\examples{
  IUPAC_CODE_MAP
  lettersToAmbiguity(c("A", "AC", "MAAMC", "GM", "AB", "MK"))
}

\keyword{utilities}
\keyword{data}
