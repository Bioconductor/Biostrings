\name{needwunsQS}
\alias{needwunsQS}

\alias{needwunsQS,character,character-method}
\alias{needwunsQS,character,BString-method}
\alias{needwunsQS,BString,character-method}
\alias{needwunsQS,BString,BString-method}
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
or a \link{BString} instance}
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
\author{Vince Carey <stvjc@channing.harvard.edu>}
%\note{  }


%\seealso{  }

\examples{
data(blosum50)
needwunsQS( "PAWHEAE", "HEAGAWGHEE", blosum50 )
}
\keyword{ models }
