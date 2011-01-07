\name{IUPAC_CODE_MAP}

\alias{IUPAC_CODE_MAP}
\alias{mergeIUPACLetters}


\title{The IUPAC Extended Genetic Alphabet}

\description{
  The \code{IUPAC_CODE_MAP} named character vector contains the mapping from
  the IUPAC nucleotide ambiguity codes to their meaning.

  The \code{mergeIUPACLetters} function provides the reverse mapping.
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
  IUPAC nucleotide ambiguity codes are used for representing
  sequences of nucleotides where the exact nucleotides that occur
  at some given positions are not known with certainty.
}

\value{
  \code{IUPAC_CODE_MAP} is a named character vector where the names are
  the IUPAC nucleotide ambiguity codes and the values are their corresponding
  meanings. The meaning of each code is described by a string that enumarates
  the base letters (\code{"A"}, \code{"C"}, \code{"G"} or \code{"T"})
  associated with the code.

  The value returned by \code{mergeIUPACLetters} is an unnamed character
  vector of the same length as its argument \code{x} where each element
  is an IUPAC nucleotide ambiguity code.
}

\author{H. Pages}

\references{
  \url{http://www.chick.manchester.ac.uk/SiteSeer/IUPAC\_codes.html}

  IUPAC-IUB SYMBOLS FOR NUCLEOTIDE NOMENCLATURE:
  Cornish-Bowden (1985) \emph{Nucl. Acids Res.} 13: 3021-3030.
}

\seealso{
  \code{\link{DNAString}}, \code{\link{RNAString}}
}

\examples{
  IUPAC_CODE_MAP
  some_iupac_codes <- c("R", "M", "G", "N", "V")
  IUPAC_CODE_MAP[some_iupac_codes]
  mergeIUPACLetters(IUPAC_CODE_MAP[some_iupac_codes])

  mergeIUPACLetters(c("Ca", "Acc", "aA", "MAAmC", "gM", "AB", "bS", "mk"))
}

\keyword{utilities}
\keyword{data}
