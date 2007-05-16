\name{blosum50}
\alias{blosum50}
\docType{data}
\title{ scoring matrix for codon substitution in protein alignments}
\description{
scoring matrix for codon substitution in protein alignments
}
\usage{data(blosum50)}
\format{
  A symmetric data frame with 24 rows and columns.
  }
}
\details{
  taken from the web; many sources available
}
\examples{
data(blosum50)
blosum50["A","A"]
blosum50["A","B"]
}
\keyword{datasets}
