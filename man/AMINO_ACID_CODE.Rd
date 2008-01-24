\name{AMINO_ACID_CODE}
\alias{AMINO_ACID_CODE}


\title{The Single-Letter Amino Acid Code}

\description{
  Named character vector mapping single-letter amino acid representations
  to 3-letter amino acid representations.
}

\seealso{
  \code{\link{AAString}},
  \code{\link{GENETIC_CODE}}
}

\examples{
  ## See all the 3-letter codes
  AMINO_ACID_CODE

  ## Convert an AAString object to a vector of 3-letter amino acid codes
  aa <- AAString("LANDEECQW")
  AMINO_ACID_CODE[strsplit(as.character(aa), NULL)[[1]]]
}

\keyword{utilities}
\keyword{data}
