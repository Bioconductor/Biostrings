\name{needwunsQS}
\alias{needwunsQS}

\alias{needwunsQS,character,character-method}
\alias{needwunsQS,character,XString-method}
\alias{needwunsQS,XString,character-method}
\alias{needwunsQS,XString,XString-method}
\alias{print.needwunsQS}
\alias{alignScore.needwunsQS}
\alias{alignScore}

\title{quadratic-space, simple gap implementation of Needleman-Wunsch
global alignment algorithm  }
\description{quadratic-space, simple gap implementation of Needleman-Wunsch
global alignment algorithm
}
\usage{
needwunsQS(s1, s2, substmat, gappen=8)
}
\arguments{
  \item{s1}{ s1: first string, an R character vector of length 1
or an \link{XString} object}
  \item{s2}{ s2: second string, same data type as \code{s1}}
  \item{substmat}{ substmat: matrix of alignment score values}
  \item{gappen}{ gappen: penalty for introducing a gap in the alignment}
}
\details{
follows specification of Durbin, Eddy, Krogh, Mitchison (1998).
}
\value{
an instance of the S3 class ``needwunsQS'', for which there are
print and alignScore methods.
}
\references{R. Durbin, S. Eddy, A. Krogh, G. Mitchison, 
Biological Sequence Analysis, Cambridge UP 1998, sec 2.3.  }
\author{Vince Carey <stvjc@channing.harvard.edu> (original author) and H. Pages (current maintainer, please send questions to me)}
%\note{  }


\seealso{
  \link{XStringAlign-class},
  \link{substitution.matrices}
}

\examples{
  ## nucleotide alignment
  mat <- matrix(-5L, nrow=4, ncol=4)
  for (i in seq_len(4)) mat[i, i] <- 0L
  rownames(mat) <- colnames(mat) <- DNA_ALPHABET[1:4]
  s1 <- DNAString(paste(sample(DNA_ALPHABET[1:4], 1000, replace=TRUE), collapse=""))
  s2 <- DNAString(paste(sample(DNA_ALPHABET[1:4], 1000, replace=TRUE), collapse=""))
  nw0 <- needwunsQS(s1, s2, mat, gappen=0)
  nw1 <- needwunsQS(s1, s2, mat, gappen=1)
  nw5 <- needwunsQS(s1, s2, mat, gappen=5)

  ## amino acid alignment
  needwunsQS("PAWHEAE", "HEAGAWGHEE", substmat="BLOSUM50")
}

\keyword{models}
\keyword{methods}
