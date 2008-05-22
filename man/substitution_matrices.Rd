\name{substitution.matrices}
\docType{data}

\alias{substitution.matrices}
\alias{BLOSUM45}
\alias{BLOSUM50}
\alias{BLOSUM62}
\alias{BLOSUM80}
\alias{BLOSUM100}
\alias{PAM30}
\alias{PAM40}
\alias{PAM70}
\alias{PAM120}
\alias{PAM250}
\alias{qualitySubstitutionMatrices}


\title{Scoring matrices}

\description{
  Predefined substitution matrices for nucleotide and amino acid alignments.
}

\usage{
  data(BLOSUM45)
  data(BLOSUM50)
  data(BLOSUM62)
  data(BLOSUM80)
  data(BLOSUM100)
  data(PAM30)
  data(PAM40)
  data(PAM70)
  data(PAM120)
  data(PAM250)
  qualitySubstitutionMatrices(alphabetLength = 4L, bitScale = 1)
}

\arguments{
  \item{alphabetLength}{an integer representing the number of letters
    in the underlying string alphabet. For DNA and RNA, this would be 4L.
    For Amino Acids, this could be 20L.}
  \item{bitScale}{a numeric value to scale the quality-based substitution
    matrices. By default, this is 1, representing bit-scale scoring.}
}

\format{
  The BLOSUM and PAM matrices are square symmetric matrices with integer
  coefficients, whose row and column names are identical and unique: each
  name is a single letter representing a nucleotide or an amino acid.

  \code{qualitySubstitutionMatrices} produces a two element list of 100 x 100
  numeric square symmetric matrices, one for matches and one for mismatches.
  The rows and columns should be thought of as indexed from [0 - 99] to
  correspond to the [0 - 99] range of quality measures.
}

\details{
  The BLOSUM and PAM matrices are not unique. For example, the definition of
  the widely used BLOSUM62 matrix varies depending on the source, and even a
  given source can provide different versions of "BLOSUM62" without keeping
  track of the changes over time. NCBI provides many matrices here
  ftp://ftp.ncbi.nih.gov/blast/matrices/ but their definitions don't match
  those of the matrices bundled with their stand-alone BLAST software
  available here ftp://ftp.ncbi.nih.gov/blast/

  The BLOSUM45, BLOSUM62, BLOSUM80, PAM30 and PAM70 matrices were taken
  from NCBI stand-alone BLAST software.

  The BLOSUM50, BLOSUM100, PAM40, PAM120 and PAM250 matrices were taken
  from ftp://ftp.ncbi.nih.gov/blast/matrices/


  The quality matrices computed in \code{qualitySubstitutionMatrices} are
  based on the paper by Ketil Malde. For quality measure \eqn{Q} in \eqn{[0, 99]},
  the corresponding error probability is given by \eqn{\epsilon_i = 10^{-Q/10}}.
  Assuming independence within and between base reads, the combined error
  probability of a mismatch when the underlying bases do match is
  \eqn{\epsilon_c = \epsilon_1 + \epsilon_2 - (n/(n-1)) * \epsilon_1 * \epsilon_2},
  where \eqn{n} is the number of letters in the underlying alphabet. Using
  \eqn{\epsilon_c}, the substitution score when two bases match is given by
  \eqn{b * \log_2((1 - \epsilon_c) * n)} and the substitution score when two bases
  don't match is given by \eqn{b * \log_2(\epsilon_c * (n/(n-1)))}, where \eqn{b}
  is the bit-scaling for the scoring.
}

\references{
K. Malde, The effect of sequence quality on sequence alignment, Bioinformatics, Feb 23, 2008.
}

\author{H. Pages and P. Aboyoun}

\seealso{
  \code{\link{pairwiseAlignment}},
  \link{PairwiseAlignment-class},
  \link{DNAString-class},
  \link{AAString-class}
}

\examples{
  ## Examine quality-based match and mismatch scores for DNA/RNA strings
  ## in pairwiseAlignment. By default patternQuality = 22L and
  ## subjectQuality = 22L.
  qualityMatrices <- qualitySubstitutionMatrices()
  qualityMatrices[["matchMatrix"]][23, 23]
  qualityMatrices[["mismatchMatrix"]][23, 23]

  s1 <- 
    DNAString("ACTTCACCAGCTCCCTGGCGGTAAGTTGATCAAAGGAAACGCAAAGTTTTCAAG")
  s2 <-
    DNAString("GTTTCACTACTTCCTTTCGGGTAAGTAAATATATAAATATATAAAAATATAATTTTCATC")
  pairwiseAlignment(s1, s2)


  ## Align two amino acid sequences with the BLOSUM62 matrix
  aa1 <- AAString("HXBLVYMGCHFDCXVBEHIKQZ")
  aa2 <- AAString("QRNYMYCFQCISGNEYKQN")
  pairwiseAlignment(aa1, aa2, substitutionMatrix = "BLOSUM62", gapOpening = -3, gapExtension = -1)

  ## See how the gap penalty influences the alignment
  pairwiseAlignment(aa1, aa2, substitutionMatrix = "BLOSUM62", gapOpening = -6, gapExtension = -2)

  ## See how the substitution matrix influences the alignment
  pairwiseAlignment(aa1, aa2, substitutionMatrix = "BLOSUM50", gapOpening = -3, gapExtension = -1)

  ## Compare our BLOSUM62 with BLOSUM62 from ftp://ftp.ncbi.nih.gov/blast/matrices/
  data(BLOSUM62)
  BLOSUM62["Q", "Z"]
  file <- "ftp://ftp.ncbi.nih.gov/blast/matrices/BLOSUM62"
  b62 <- as.matrix(read.table(file, check.names=FALSE))
  b62["Q", "Z"]
}

\keyword{utilities}
\keyword{data}
\keyword{datasets}
