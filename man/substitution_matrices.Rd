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
\alias{nucleotideSubstitutionMatrix}
\alias{qualitySubstitutionMatrices}
\alias{errorSubstitutionMatrices}


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

nucleotideSubstitutionMatrix(match=1, mismatch=0, baseOnly=FALSE,
                             type="DNA", symmetric=TRUE)

qualitySubstitutionMatrices(fuzzyMatch=c(0, 1), alphabetLength=4L,
                            qualityClass="PhredQuality", bitScale=1)

errorSubstitutionMatrices(errorProbability, fuzzyMatch=c(0, 1),
                          alphabetLength=4L, bitScale=1)
}

\arguments{
  \item{match}{the scoring for a nucleotide match.}
  \item{mismatch}{the scoring for a nucleotide mismatch.}
  \item{baseOnly}{\code{TRUE} or \code{FALSE}. If \code{TRUE}, only
    uses the letters in the "base" alphabet i.e. "A", "C", "G", "T".}
  \item{type}{either "DNA" or "RNA".}
  \item{symmetric}{\code{TRUE} or \code{FALSE}. Default is \code{TRUE}. 
    If \code{FALSE}, the resulting matrix will be asymmetric.}
  \item{fuzzyMatch}{a named or unnamed numeric vector representing the base
    match probability.}
  \item{errorProbability}{a named or unnamed numeric vector representing
    the error probability.}
  \item{alphabetLength}{an integer representing the number of letters
    in the underlying string alphabet. For DNA and RNA, this would be 4L.
    For Amino Acids, this could be 20L.}
  \item{qualityClass}{a character string of \code{"PhredQuality"},
    \code{"SolexaQuality"}, or \code{"IlluminaQuality"}.}
  \item{bitScale}{a numeric value to scale the quality-based substitution
    matrices. By default, this is 1, representing bit-scale scoring.}
}

\format{
  The BLOSUM and PAM matrices are square symmetric matrices with integer
  coefficients, whose row and column names are identical and unique: each
  name is a single letter representing a nucleotide or an amino acid.

  \code{nucleotideSubstitutionMatrix} produces a substitution matrix for
  all IUPAC nucleic acid codes based upon match and mismatch parameters.

  \code{errorSubstitutionMatrices} produces a two element list of numeric
  square symmetric matrices, one for matches and one for mismatches.

  \code{qualitySubstitutionMatrices} produces the substitution matrices
  for Phred or Solexa quality-based reads.
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
  based on the paper by Ketil Malde. Let \eqn{\epsilon_i} be the probability
  of an error in the base read. For \code{"Phred"} quality measures \eqn{Q}
  in \eqn{[0, 99]}, these error probabilities are given by
  \eqn{\epsilon_i = 10^{-Q/10}}. For \code{"Solexa"} quality measures
  \eqn{Q} in \eqn{[-5, 99]}, they are given by
  \eqn{\epsilon_i = 1 - 1/(1 + 10^{-Q/10})}.
  Assuming independence within and between base reads, the combined error
  probability of a mismatch when the underlying bases do match is
  \eqn{\epsilon_c = \epsilon_1 + \epsilon_2 - (n/(n-1)) * \epsilon_1 * \epsilon_2},
  where \eqn{n} is the number of letters in the underlying alphabet. Using
  \eqn{\epsilon_c}, the substitution score is given by when two bases match is given by
  \eqn{b * \log_2(\gamma_{x,y} * (1 - \epsilon_c) * n + (1 - \gamma_{x,y}) * \epsilon_c * (n/(n-1)))},
  where \eqn{b} is the bit-scaling for the scoring and \eqn{\gamma_{x,y}} is the
  probability that characters \eqn{x} and \eqn{y} represents the same underlying
  information (e.g. using IUPAC, \eqn{\gamma_{A,A} = 1} and \eqn{\gamma_{A,N} = 1/4}.
  In the arguments listed above \code{fuzzyMatch} represents \eqn{\gamma_{x,y}}
  and \code{errorProbability} represents \eqn{\epsilon_i}.
}

\references{
K. Malde, The effect of sequence quality on sequence alignment, Bioinformatics, Feb 23, 2008.
}

\author{H. Pagès and P. Aboyoun, with contribution from Albert Vill (support
        for asymmetric matrices in \code{nucleotideSubstitutionMatrix()})}

\seealso{
  \code{\link{pairwiseAlignment}},
  \link{PairwiseAlignments-class},
  \link{DNAString-class},
  \link{AAString-class},
  \link{PhredQuality-class},
  \link{SolexaQuality-class},
  \link{IlluminaQuality-class}
}

\examples{
s1 <- DNAString("ACTTCACCAGCTCCCTGGCGGTAAGTTGATCAAAGGAAACGCAAAGTTTTCAAG")
s2 <- DNAString("GTTTCACTACTTCCTTTCGGGTAAGTAAATATATAAATATATAAAAATATAATTTTCATC")
s3 <- DNAString("NCTTCRCCAGCTCCCTGGMGGTAAGTTGATCAAAGVAAACGCAAAGTTNTCAAG")

## Fit a global pairwise alignment using edit distance scoring:
nsm <- nucleotideSubstitutionMatrix(0, -1, TRUE)
pairwiseAlignment(s1, s2, substitutionMatrix=nsm,
                          gapOpening=0, gapExtension=1)

## Align sequences using an asymmetric substitution matrix.
## The asymmetry of the matrix means that the query sequence is not
## penalized for ambiguous bases in the subject / consensus sequence:
ansm <- nucleotideSubstitutionMatrix(symmetric=FALSE)
ansm["M", c("A","C","G","T")]
#   A   C   G   T 
# 0.5 0.5 0.0 0.0 
ansm[c("A","C","G","T"), "M"]
# A C G T 
# 1 1 0 0 
ansm["M", "H"]
# 1
ansm["H", "M"]
# 0.6666667

## Due to this asymmetry, the order of the sequences is important:
pairwiseAlignment(s1, s3, substitutionMatrix=ansm)
pairwiseAlignment(s3, s1, substitutionMatrix=ansm)

## Examine quality-based match and mismatch bit scores for DNA/RNA
## strings in pairwiseAlignment. By default patternQuality and
## subjectQuality are PhredQuality(22L):
qualityMatrices <- qualitySubstitutionMatrices()
qualityMatrices["22", "22", "1"]
qualityMatrices["22", "22", "0"]

pairwiseAlignment(s1, s2)

## Get the substitution scores when the error probability is 0.1:
subscores <- errorSubstitutionMatrices(errorProbability=0.1)
submat <- matrix(subscores[ , , "0"], 4, 4)
diag(submat) <- subscores[ , , "1"]
dimnames(submat) <- list(DNA_ALPHABET[1:4], DNA_ALPHABET[1:4])
submat
pairwiseAlignment(s1, s2, substitutionMatrix=submat)

## Align two amino acid sequences with the BLOSUM62 matrix:
aa1 <- AAString("HXBLVYMGCHFDCXVBEHIKQZ")
aa2 <- AAString("QRNYMYCFQCISGNEYKQN")
pairwiseAlignment(aa1, aa2, substitutionMatrix="BLOSUM62",
                            gapOpening=3, gapExtension=1)

## See how the gap penalty influences the alignment:
pairwiseAlignment(aa1, aa2, substitutionMatrix="BLOSUM62",
                            gapOpening=6, gapExtension=2)

## See how the substitution matrix influences the alignment:
pairwiseAlignment(aa1, aa2, substitutionMatrix="BLOSUM50",
                            gapOpening=3, gapExtension=1)

if (interactive()) {
  ## Compare our BLOSUM62 with BLOSUM62 from
  ## ftp://ftp.ncbi.nih.gov/blast/matrices/:
  data(BLOSUM62)
  BLOSUM62["Q", "Z"]
  file <- "ftp://ftp.ncbi.nih.gov/blast/matrices/BLOSUM62"
  b62 <- as.matrix(read.table(file, check.names=FALSE))
  b62["Q", "Z"]
}
}

\keyword{utilities}
\keyword{data}
\keyword{datasets}
