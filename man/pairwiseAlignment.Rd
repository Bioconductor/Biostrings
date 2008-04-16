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
pairwiseAlignment(string1, string2, matchScores, gapExtension = -8L,
                  gapOpening = 0L, type = "global")
}
\arguments{
  \item{string1, string2}{an R character vector of length 1 or an \link{XString} object.}
  \item{matchScores}{matrix of alignment score values.}
  \item{gapExtension}{penalty for extending a gap in the alignment.}
  \item{gapOpening}{penalty for opening a gap in the alignment. (not used)}
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
  ## nucleotide global alignment
  mat <- matrix(-5L, nrow = 4, ncol = 4)
  for (i in seq_len(4)) mat[i, i] <- 0L
  rownames(mat) <- colnames(mat) <- DNA_ALPHABET[1:4]
  s1 <- DNAString(paste(sample(DNA_ALPHABET[1:4], 1000, replace=TRUE), collapse=""))
  s2 <- DNAString(paste(sample(DNA_ALPHABET[1:4], 1000, replace=TRUE), collapse=""))
  nw0 <- pairwiseAlignment(s1, s2, mat, gapExtension = 0L)
  nw1 <- pairwiseAlignment(s1, s2, mat, gapExtension = -1L)
  nw5 <- pairwiseAlignment(s1, s2, mat, gapExtension = -5L)

  ## amino acid local alignment
  pairwiseAlignment("PAWHEAE", "HEAGAWGHEE", matchScores = "BLOSUM50", type = "local")
}
\keyword{models}
\keyword{methods}
