\name{scoring.matrices}
\docType{data}

\alias{scoring.matrices}
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


\title{Scoring matrices}

\description{
  Predefined substitution scoring matrices for nucleotide
  and amino acid alignments.
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
}

\format{
  A square symetric matrix with integer coefficients.
  The row and column names are identical and unique: each name
  is a single letter representing a nucleotide or an amino acid.
}

\details{
  Note that there can exist different versions of a given scoring matrix.
  For example, definition of widely used BLOSUM62 matrix varies depending on
  the source. Even a given source can provide different versions of it but the
  name is always BLOSUM62 and they provide no history or versioning mechanism!
  NCBI for example provides many matrices here ftp://ftp.ncbi.nih.gov/blast/matrices/
  but their definitions don't match those of the matrices bundled with their
  standalone BLAST software available here ftp://ftp.ncbi.nih.gov/blast/

  The BLOSUM45, BLOSUM62, BLOSUM80, PAM30 and PAM70 matrices were taken
  from NCBI standalone BLAST software.

  The BLOSUM50, BLOSUM100, PAM40, PAM120 and PAM250 matrices were taken
  from ftp://ftp.ncbi.nih.gov/blast/matrices/
}

\seealso{
  \code{\link{needwunsQS}},
  \link{BStringAlign-class},
  \link{DNAString-class},
  \link{AAString-class}
}

\examples{
  ## Align 2 amino acid sequences with the BLOSUM62 matrix
  aa1 <- AAString("HXBLVYMGCHFDCXVBEHIKQZ")
  aa2 <- AAString("QRNYMYCFQCISGNEYKQN")
  needwunsQS(aa1, aa2, "BLOSUM62", gappen=3)

  ## See how the gap penalty influences the alignment
  needwunsQS(aa1, aa2, "BLOSUM62", gappen=8)

  ## See how the scoring matrix influences the alignment
  needwunsQS(aa1, aa2, "BLOSUM50", gappen=3)

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

