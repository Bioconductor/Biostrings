\name{alphabetFrequency}

\alias{alphabetFrequency}
\alias{alphabetFrequency,BString-method}
\alias{alphabetFrequency,DNAString-method}
\alias{alphabetFrequency,RNAString-method}
\alias{alphabetFrequency,BStringViews-method}

\alias{allOligonucleotideFrequency}
\alias{allOligonucleotideFrequency,DNAString-method}
\alias{allOligonucleotideFrequency,RNAString-method}
\alias{allOligonucleotideFrequency,BStringViews-method}

\alias{dinucleotideFrequency}
\alias{trinucleotideFrequency}

\alias{strrev}
\alias{mkAllStrings}


\title{Functions to calculate the frequency of letters in a biological sequence}

\description{
  Given a biological sequence, the \code{alphabetFrequency} function will
  calculate the frequency of each letter in the (base) alphabet,
  the \code{dinucleotideFrequency} function the frequency of all possible
  dinucleotides and the \code{trinucleotideFrequency} function the frequency
  of all possible trinucleotides.

  More generally, the \code{allOligonucleotideFrequency} function will
  calculate the frequency of all possible oligonucleotides of a given
  length (called the "width" in this particular context).
}

\usage{
  alphabetFrequency(x, baseOnly=FALSE, freq=FALSE)
  dinucleotideFrequency(x, freq=FALSE, fast.moving.side="right", as.matrix=FALSE)
  trinucleotideFrequency(x, freq=FALSE, fast.moving.side="right", as.array=FALSE)
  allOligonucleotideFrequency(x, width, freq=FALSE, fast.moving.side="right", as.array=FALSE)

  ## Some related utility functions
  strrev(x)
  mkAllStrings(alphabet, width, fast.moving.side="right")
}

\arguments{
  \item{x}{
    A \link{BString}, \link{DNAString}, \link{RNAString}
    or \link{BStringViews} object for all the \code{*Frequency}
    functions and a character vector for \code{strrev}.
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
  \item{freq}{
    If \code{TRUE} then frequencies are reported, otherwise counts.
  }
  \item{fast.moving.side}{
    Which side of the strings should move fastest?
  }
  \item{as.matrix}{
    If \code{TRUE} then return a numeric matrix, otherwise a numeric
    vector with no dim attribute.
  }
  \item{as.array}{
    If \code{TRUE} then return a numeric array, otherwise a numeric
    vector with no dim attribute.
  }
  \item{width}{
    The number of nucleotides per oligonucleotide for
    \code{allOligonucleotideFrequency}.
    The number of letters per string for \code{mkAllStrings}.
  }
  \item{alphabet}{
    The alphabet to use to make the strings.
  }
}

\value{
  All the \code{*Frequency} functions return an integer vector if \code{freq}
  is \code{FALSE} (default), otherwise double vector.
  If \code{as.matrix} or \code{as.array} is \code{TRUE}, this vector is
  formatted as a matrix or an array.

  For \code{alphabetFrequency}: if \code{x} is a \link{DNAString} or
  \link{RNAString} object (or a \link{BStringViews} object with a
  \link{DNAString} or \link{RNAString} subject), then the returned vector
  is named with the letters in the alphabet.
  If the \code{baseOnly} argument is TRUE, then the returned vector
  has an extra element named 'other'.

  \code{dinucleotideFrequency} (resp. \code{trinucleotideFrequency} and
  \code{allOligonucleotideFrequency}) only works on a \link{DNAString} or
  \link{RNAString} object (or a \link{BStringViews} object with a
  \link{DNAString} or \link{RNAString} subject)
  and returns a vector named with all the possible dinucleotides
  (resp. trinucleotides or oligonucleotides).
  
  If \code{x} is a \link{BStringViews} object, then the returned
  vector contains the cumulated frequencies found in each view.
}

\author{H. Pages}

\seealso{
  \link{BString}, \link{DNAString}, \link{RNAString}, \link{BStringViews},
  \code{\link{rev}}, \code{\link{strsplit}}
}

\examples{
  data(yeastSEQCHR1)
  yeast1 <- DNAString(yeastSEQCHR1)

  alphabetFrequency(yeast1)
  alphabetFrequency(yeast1, baseOnly=TRUE)

  dinucleotideFrequency(yeast1)
  trinucleotideFrequency(yeast1)
  allOligonucleotideFrequency(yeast1, 4)

  ## Get the less and most represented 6-mers
  f6 <- allOligonucleotideFrequency(yeast1, 6)
  f6[f6 == min(f6)]
  f6[f6 == max(f6)]

  ## Get the result as an array
  triA <- trinucleotideFrequency(yeast1, as.array=TRUE)
  triA["A", "A", "C"] # == trinucleotideFrequency(yeast1)["AAC"]

  ## Note that when dropping the dimensions of the triA array, elements
  ## in the resulting vector are ordered as if they had been obtained
  ## with fast.moving.side="left":
  triL <- trinucleotideFrequency(yeast1, fast.moving.side="left")
  all(as.vector(triA) == triL)

  ## Some related utility functions
  dict1 <- mkAllStrings(LETTERS[1:3], 4)
  dict2 <- mkAllStrings(LETTERS[1:3], 4, fast.moving.side="left")
  identical(strrev(dict1), dict2) # TRUE 
}

\keyword{category}

