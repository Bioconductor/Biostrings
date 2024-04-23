\name{needwunsQS}

\alias{needwunsQS}

\title{(Defunct) Needleman-Wunsch Global Alignment}
\description{Simple gap implementation of Needleman-Wunsch global alignment algorithm.}
\usage{
needwunsQS(s1, s2, substmat, gappen = 8)
}
\arguments{
  \item{s1, s2}{an R character vector of length 1 or an \link{XString} object.}
  \item{substmat}{matrix of alignment score values.}
  \item{gappen}{penalty for introducing a gap in the alignment.}
}
\details{
Follows specification of Durbin, Eddy, Krogh, Mitchison (1998).
This function is now defunct. Please use
\code{\link[pwalign]{pairwiseAlignment}} from the \pkg{pwalign} instead.
}
\value{
An instance of class \code{"PairwiseAlignments"}.
}
\references{R. Durbin, S. Eddy, A. Krogh, G. Mitchison,
Biological Sequence Analysis, Cambridge UP 1998, sec 2.3.}
\author{Vince Carey (\email{stvjc@channing.harvard.edu}) (original author) and H. Pagès (current maintainer).}
\seealso{
  \code{\link[pwalign]{pairwiseAlignment}} and
  \link[pwalign]{PairwiseAlignments-class} in the \pkg{pwalign} package,
  \link[pwalign]{substitution_matrices}
}
\examples{
\dontrun{
  ## This function is now defunct.
  ## Please use pairwiseAlignment() from the pwalign package instead.

  ## nucleotide alignment
  mat <- matrix(-5L, nrow = 4, ncol = 4)
  for (i in seq_len(4)) mat[i, i] <- 0L
  rownames(mat) <- colnames(mat) <- DNA_ALPHABET[1:4]
  s1 <- DNAString(paste(sample(DNA_ALPHABET[1:4], 1000, replace=TRUE), collapse=""))
  s2 <- DNAString(paste(sample(DNA_ALPHABET[1:4], 1000, replace=TRUE), collapse=""))
  nw0 <- needwunsQS(s1, s2, mat, gappen = 0)
  nw1 <- needwunsQS(s1, s2, mat, gappen = 1)
  nw5 <- needwunsQS(s1, s2, mat, gappen = 5)

  ## amino acid alignment
  needwunsQS("PAWHEAE", "HEAGAWGHEE", substmat = "BLOSUM50")
}
}
\keyword{models}
\keyword{methods}
