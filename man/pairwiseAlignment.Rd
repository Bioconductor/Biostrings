\name{pairwiseAlignment}
\alias{pairwiseAlignment}

\alias{pairwiseAlignment,character,character-method}
\alias{pairwiseAlignment,character,XString-method}
\alias{pairwiseAlignment,XString,character-method}
\alias{pairwiseAlignment,XString,XString-method}

\title{Optimal Pairwise Alignment}
\description{
Simple gap implementation of Needleman-Wunsch global alignment and
Smith-Waterman local alignment algorithms.
}
\usage{
pairwiseAlignment(string1, string2, matchScoring, gapOpening = -5L,
                  gapExtension = -2L, type = "global")
}
\arguments{
  \item{string1, string2}{an R character vector of length 1 or an \link{XString} object.}
  \item{matchScoring}{matrix of alignment score values.}
  \item{gapOpening}{penalty for opening a gap in the alignment.}
  \item{gapExtension}{penalty for extending a gap in the alignment.}
  \item{type}{type of alignment (\code{"global"}, \code{"local"}, \code{"overlap"}).}
}
\details{
Follows specification of Haubold and Wiehe (2006).
}
\value{
An instance of class \code{"XStringAlign"}.
}
\references{B. Haubold, T. Wiehe, Introduction to Computational Biology, Birkhauser Verlag 2006, Chapter 2.}
\author{Patrick Aboyoun and Herve Pages.}
\seealso{
  \link{XStringAlign-class},
  \link{substitution.matrices}
}
\examples{
  ## Nucleotide global, local, and overlap alignments
  mat <- matrix(-3L, nrow = 4, ncol = 4)
  diag(mat) <- 1L
  rownames(mat) <- colnames(mat) <- DNA_ALPHABET[1:4]
  s1 <- 
    DNAString("ACTTCACCAGCTCCCTGGCGGTAAGTTGATCAAAGGAAACGCAAAGTTTTCAAG")
  s2 <-
    DNAString("GTTTCACTACTTCCTTTCGGGTAAGTAAATATATAAATATATAAAAATATAATTTTCATC")
  globalAlign <- pairwiseAlignment(s1, s2, mat)
  localAlign <- pairwiseAlignment(s1, s2, mat, type = "local")
  overlapAlign <- pairwiseAlignment(s1, s2, mat, type = "overlap")
  pairwiseAlignment(s1, s2, mat, gapOpening = -3L, gapExtension = -1L)

  ## Amino acid global alignment
  pairwiseAlignment("PAWHEAE", "HEAGAWGHEE", matchScoring = "BLOSUM50")
}
\keyword{models}
\keyword{methods}
