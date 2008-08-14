\name{GENETIC_CODE}

\alias{GENETIC_CODE}
\alias{RNA_GENETIC_CODE}

\title{The Standard Genetic Code}

\description{
  Two predefined objects (\code{GENETIC_CODE} and \code{RNA_GENETIC_CODE})
  that represent The Standard Genetic Code.
}

\usage{
    GENETIC_CODE
    RNA_GENETIC_CODE
}

\details{
  Formally, a genetic code is a mapping between tri-nucleotide sequences
  called codons, and amino acids.

  The Standard Genetic Code (aka The Canonical Genetic Code, or simply The
  Genetic Code) is the particular mapping that encodes the vast majority of
  genes in nature.

  \code{GENETIC_CODE} and \code{RNA_GENETIC_CODE} are predefined named
  character vectors that represent this mapping.
}

\value{
  \code{GENETIC_CODE} and \code{RNA_GENETIC_CODE} are both named character
  vectors of length 64 (the number of all possible tri-nucleotide sequences)
  where each element is a single letter representing either an amino acid
  or the stop codon \code{"*"} (aka termination codon).

  The names of the \code{GENETIC_CODE} vector are the DNA codons i.e. the
  tri-nucleotide sequences (directed 5' to 3') that are assumed to belong
  to the "coding DNA strand" (aka "sense DNA strand" or "non-template DNA
  strand") of the gene.
 
  The names of the \code{RNA_GENETIC_CODE} are the RNA codons i.e. the
  tri-nucleotide sequences (directed 5' to 3') that are assumed to belong
  to the mRNA of the gene.

  Note that the values in the \code{GENETIC_CODE} and \code{RNA_GENETIC_CODE}
  vectors are the same, only their names are different. The names of the
  latter are those of the former where all occurences of T (thymine) have
  been replaced by U (uracil).
}

\author{H. Pages}

\references{
  \url{http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi}
}

\seealso{
  \code{\link{AA_ALPHABET}},
  \code{\link{AMINO_ACID_CODE}},
  \code{\link{trinucleotideFrequency}},
  \code{\link{DNAString}}, \code{\link{RNAString}}, \code{\link{AAString}}
}

\examples{
  GENETIC_CODE
  RNA_GENETIC_CODE
  all(GENETIC_CODE == RNA_GENETIC_CODE)  # TRUE
}

\keyword{utilities}
\keyword{data}
