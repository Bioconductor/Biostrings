\name{pairwiseAlignment}
\alias{pairwiseAlignment}

\alias{pairwiseAlignment,character,character-method}
\alias{pairwiseAlignment,character,XString-method}
\alias{pairwiseAlignment,XString,character-method}
\alias{pairwiseAlignment,XString,XString-method}

\title{Optimal Pairwise Alignment}
\description{
Solves the (Needleman-Wunsch) global alignment and (Smith-Waterman) local alignment problems.
}
\usage{
pairwiseAlignment(string1, string2, quality1 = 1, quality2 = 1,
                  matchScoring, gapOpening = -5, gapExtension = -2, type = "global")
}
\arguments{
  \item{string1, string2}{an R character vector of length 1 or an \link{XString} object.}
  \item{quality1, quality2}{quality measures for \code{string1} and \code{string2}
    represented by [0, 1] numeric vectors with either one element or one element
    for each character in \code{string1} and \code{string2} respectively.}
  \item{matchScoring}{matrix of alignment score values.}
  \item{gapOpening}{penalty for opening a gap in the alignment.}
  \item{gapExtension}{penalty for extending a gap in the alignment.}
  \item{type}{type of alignment (\code{"global"}, \code{"local"}, \code{"overlap"}).}
}
\details{
Implementation based on Chapter 2 of Haubold and Wiehe (2006).
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
  globalAlign <- pairwiseAlignment(s1, s2, matchScoring = mat)
  localAlign <- pairwiseAlignment(s1, s2, matchScoring = mat, type = "local")
  overlapAlign <- pairwiseAlignment(s1, s2, matchScoring = mat, type = "overlap")
  pairwiseAlignment(s1, s2,
                    quality1 = rep(c(0.85, 0.65), times = c(36, 18)),
                    quality2 = rep(c(0.75, 0.50), times = c(40, 20)),
                    matchScoring = mat, gapOpening = -3, gapExtension = -1)

  ## Amino acid global alignment
  pairwiseAlignment("PAWHEAE", "HEAGAWGHEE", matchScoring = "BLOSUM50")
}
\keyword{models}
\keyword{methods}
