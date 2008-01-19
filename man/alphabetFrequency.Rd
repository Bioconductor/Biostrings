\name{alphabetFrequency}

\alias{alphabetFrequency}
\alias{alphabetFrequency,BString-method}
\alias{alphabetFrequency,DNAString-method}
\alias{alphabetFrequency,RNAString-method}
\alias{alphabetFrequency,BStringViews-method}

\alias{dinucleotideFrequency}
\alias{dinucleotideFrequency,DNAString-method}
\alias{dinucleotideFrequency,RNAString-method}
\alias{dinucleotideFrequency,BStringViews-method}

\alias{trinucleotideFrequency}
\alias{trinucleotideFrequency,DNAString-method}
\alias{trinucleotideFrequency,RNAString-method}
\alias{trinucleotideFrequency,BStringViews-method}

\alias{mkAllStrings}


\title{Functions to calculate the frequency of letters in a biological sequence}

\description{
  Given a biological sequence, the \code{alphabetFrequency} function will
  calculate the frequency of each letter in the (base) alphabet,
  the \code{dinucleotideFrequency} function the frequency of all possible
  dinucleotides and the \code{trinucleotideFrequency} function the frequency
  of all possible trinucleotides.
}

\usage{
  alphabetFrequency(x, baseOnly=FALSE)
  dinucleotideFrequency(x)
  trinucleotideFrequency(x)

  ## Some related utility function
  mkAllStrings(alphabet, width)
}

\arguments{
  \item{x}{
    A \link{BString}, \link{DNAString}, \link{RNAString}
    or \link{BStringViews} object.
  }
  \item{baseOnly}{
    \code{TRUE} or \code{FALSE}.
    If \code{TRUE}, the returned vector only contains frequencies for the
    letters in the "base" alphabet i.e. "A", "C", "G", "T" when \code{x}
    is a \link{DNAString} object (or a \link{BStringViews} object with a
    \link{DNAString} subject), and "A", "C", "G", "U" when \code{x} is a
    \link{RNAString} object (or a \link{BStringViews} object with a
    \link{RNAString} subject).
    When \code{x} is a \link{BString} object (or a \link{BStringViews}
    object with a \link{BString} subject), then the \code{baseOnly}
    argument is ignored.
  }
}

\value{
  An integer vector.

  For \code{alphabetFrequency}: if \code{x} is a \link{DNAString} or
  \link{RNAString} object (or a \link{BStringViews} object with a
  \link{DNAString} or \link{RNAString} subject), then the returned vector
  is named with the letters in the alphabet.
  If the \code{baseOnly} argument is TRUE, then the returned vector
  has an extra element named 'other'.

  \code{dinucleotideFrequency} and \code{trinucleotideFrequency} only
  work on a \link{DNAString} or \link{RNAString} object (or a
  \link{BStringViews} object with a \link{DNAString} or \link{RNAString}
  subject) and return a vector named with all the possible dinucleotides
  or trinucleotides.
  
  If \code{x} is a \link{BStringViews} object, then the returned
  vector contains the cumulated frequencies found in each view.
}

\seealso{
  \link{BString}, \link{DNAString}, \link{RNAString}, \link{BStringViews}
}

\examples{
  data(yeastSEQCHR1)
  yeast1 <- DNAString(yeastSEQCHR1)

  alphabetFrequency(yeast1)
  alphabetFrequency(yeast1, baseOnly=TRUE)

  dinucleotideFrequency(yeast1)
  trinucleotideFrequency(yeast1)

  mkAllStrings(LETTERS[1:3], 4)
}

\keyword{category}

