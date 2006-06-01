\name{alphabetFrequency}

\alias{alphabetFrequency}
\alias{alphabetFrequency,BString-method}
\alias{alphabetFrequency,DNAString-method}
\alias{alphabetFrequency,RNAString-method}
\alias{alphabetFrequency,BStringViews-method}

\alias{alphabetFrequency}


\title{Function to calculate the frequency of letters in a biological sequence}

\description{
  This function calculates the frequency of each letter in the (base) alphabet
  for a biological sequence.
}

\usage{
  alphabetFrequency(x, baseOnly=TRUE)
}

\arguments{
  \item{x}{
    A \link{BString}, \link{DNAString}, \link{RNAString}
    or \link{BStringViews} object.
  }
  \item{baseOnly}{
    \code{TRUE} or \code{FALSE}.
    If \code{TRUE}, the returned vector only contains frequencies for the
    letters in the "base" alphabet + the gap symbol
    i.e. "A", "C", "G", "T" + "-" when \code{x} is a \link{DNAString} object
    (or a \link{BStringViews} object with a \link{DNAString} subject),
    and "U", "G", "C", "A" + "-" when \code{x} is a \link{RNAString} object
    (or a \link{BStringViews} object with a \link{RNAString} subject).
    When \code{x} is a \link{BString} object (or a \link{BStringViews}
    object with a \link{BString} subject), then the \code{baseOnly}
    argument is ignored.
  }
}

\value{
  An integer vector.
  If \code{x} is a \link{DNAString} or \link{RNAString} object
  (or a \link{BStringViews} object with a \link{DNAString}
  or \link{RNAString} subject), then the returned vector is named
  with the letters in the alphabet.
  If the \code{baseOnly} argument is TRUE, then the returned vector
  has an extra element named 'other'.
  If \code{x} is a \link{BStringViews} object, then the returned
  vector contains the cumulated frequencies found in each view.
}

\author{H. Pages}

\seealso{
  \link{BString}, \link{DNAString}, \link{RNAString}, \link{BStringViews}
}

\examples{
  data(yeastSEQCHR1)
  yeast1 <- DNAString(yeastSEQCHR1)
  alphabetFrequency(yeast1)
  alphabetFrequency(yeast1, baseOnly=FALSE)
}

\keyword{category}

