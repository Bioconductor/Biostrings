\name{needwunsQS}
\alias{needwunsQS}

\alias{needwunsQS,character,character-method}
\alias{needwunsQS,character,XString-method}
\alias{needwunsQS,XString,character-method}
\alias{needwunsQS,XString,XString-method}
\alias{print.needwunsQS}

\title{(Deprecated) Needleman-Wunsch Global Alignment}
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
This function has been deprecated and is being replaced by
\code{pairwiseAlignment}.
}
\value{
An instance of class \code{"PairwiseAlignments"}.
}
\references{R. Durbin, S. Eddy, A. Krogh, G. Mitchison,
Biological Sequence Analysis, Cambridge UP 1998, sec 2.3.}
\author{Vince Carey (\email{stvjc@channing.harvard.edu}) (original author) and H. Pages (current maintainer).}
\seealso{
  \link{pairwiseAlignment},
  \link{PairwiseAlignments-class},
  \link{substitution.matrices}
}
\examples{
\dontrun{
  ## This function has been deprecated
  ## Use 'pairwiseAlignment' instead.

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
