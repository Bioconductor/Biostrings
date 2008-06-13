\name{alphabetFrequency}

\alias{alphabetFrequency}
\alias{alphabetFrequency,XString-method}
\alias{alphabetFrequency,DNAString-method}
\alias{alphabetFrequency,RNAString-method}
\alias{alphabetFrequency,XStringSet-method}
\alias{alphabetFrequency,DNAStringSet-method}
\alias{alphabetFrequency,RNAStringSet-method}
\alias{alphabetFrequency,XStringViews-method}
\alias{alphabetFrequency,MaskedXString-method}

\alias{oligonucleotideFrequency}
\alias{oligonucleotideFrequency,DNAString-method}
\alias{oligonucleotideFrequency,RNAString-method}
\alias{oligonucleotideFrequency,XStringSet-method}
\alias{oligonucleotideFrequency,XStringViews-method}
\alias{oligonucleotideFrequency,MaskedXString-method}

\alias{dinucleotideFrequency}
\alias{trinucleotideFrequency}
\alias{oligonucleotideTransitions}

\alias{strrev}
\alias{mkAllStrings}


\title{Functions to calculate the frequency of letters in a biological sequence}

\description{
  Given a biological sequence, the \code{alphabetFrequency} function will
  calculate the frequency of each letter in the (base) alphabet,
  the \code{dinucleotideFrequency} function the frequency of all possible
  dinucleotides and the \code{trinucleotideFrequency} function the frequency
  of all possible trinucleotides.

  More generally, the \code{oligonucleotideFrequency} function will
  calculate the frequency of all possible oligonucleotides of a given
  length (called the "width" in this particular context).

  In this man page we call "DNA input" a \link{DNAString} object,
  or a \link{DNAStringSet} object, or an \link{XStringViews} object
  with a \link{DNAString} subject, or a \link{MaskedDNAString} object.
  Similarly we call "RNA input" an \link{RNAString} object,
  or an \link{RNAStringSet} object, or an \link{XStringViews} object
  with an \link{RNAString} subject, or a \link{MaskedRNAString} object.
}

\usage{
  alphabetFrequency(x, baseOnly=FALSE, freq=FALSE, ...)
  dinucleotideFrequency(x, freq=FALSE, fast.moving.side="right", as.matrix=FALSE, with.labels=TRUE, ...)
  trinucleotideFrequency(x, freq=FALSE, fast.moving.side="right", as.array=FALSE, with.labels=TRUE, ...)
  oligonucleotideFrequency(x, width, freq=FALSE, fast.moving.side="right", as.array=FALSE, with.labels=TRUE, ...)
  oligonucleotideTransitions(x, left=1, right=1, freq=FALSE)

  ## Some related utility functions
  strrev(x)
  mkAllStrings(alphabet, width, fast.moving.side="right")
}

\arguments{
  \item{x}{
    An \link{XString}, \link{XStringSet}, \link{XStringViews}
    or \link{MaskedXString} object for all the \code{*Frequency} functions.
    A character vector for \code{strrev}.
  }
  \item{baseOnly}{
    \code{TRUE} or \code{FALSE}.
    If \code{TRUE}, the returned vector only contains frequencies for the
    letters in the "base" alphabet i.e. "A", "C", "G", "T" if \code{x}
    is a "DNA input", and "A", "C", "G", "U" if \code{x} is "RNA input".
    When \code{x} is a \link{BString} object (or an \link{XStringViews}
    object with a \link{BString} subject, or a \link{BStringSet} object),
    then the \code{baseOnly} argument is ignored.
  }
  \item{freq}{
    If \code{TRUE} then frequencies are reported, otherwise counts.
  }
  \item{...}{
    Further arguments to be passed to or from other methods.
    For the \link{XStringViews} and \link{XStringSet} methods,
    the \code{collapse} argument is accepted.
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
  \item{with.labels}{
    If \code{TRUE} then return a named vector (or array).
  }
  \item{width}{
    The number of nucleotides per oligonucleotide for
    \code{oligonucleotideFrequency}.
    The number of letters per string for \code{mkAllStrings}.
  }
  \item{left, right}{
    The number of nucleotides per oligonucleotide for the rows and columns
    respectively in the transition matrix created by \code{oligonucleotideTransitions}.
  }
  \item{alphabet}{
    The alphabet to use to make the strings.
  }
}

\details{
  \code{alphabetFrequency} and \code{oligonucleotideFrequency} are generic
  functions with methods in the Biostrings package defined for
  \link{BString}, \link{DNAString}, \link{RNAString}, \link{XStringViews}
  and \link{XStringSet} objects.
}

\value{
  All the \code{*Frequency} functions return an integer vector if \code{freq}
  is \code{FALSE} (default), otherwise a double vector.
  If \code{as.matrix} or \code{as.array} is \code{TRUE}, this vector is
  formatted as a matrix or an array.

  For \code{alphabetFrequency}: if \code{x} is a "DNA or RNA input", then the
  returned vector is named with the letters in the alphabet (unless
  \code{with.labels} is \code{FALSE}).
  If the \code{baseOnly} argument is \code{TRUE}, then the returned vector
  has only 5 elements: 4 elements corresponding to the 4 nucleotides
  + the 'other' element.

  \code{dinucleotideFrequency} (resp. \code{trinucleotideFrequency} and
  \code{oligonucleotideFrequency}) only works on "DNA or RNA input"
  and returns a vector named with all the possible dinucleotides
  (resp. trinucleotides or oligonucleotides).
  
  If \code{x} is a multiple sequence input (i.e. an \link{XStringViews} or
  \link{XStringSet} object), then the returned object is a matrix (or a list)
  with the same number of rows (or elements) as \code{x} unless \code{collapse=TRUE}
  is specified. In that case the returned vector (or array) contains the frequencies
  cumulated across all sequences in \code{x}.
}

\author{H. Pages}

\seealso{
  \code{\link{countPDict}},
  \link{XString-class},
  \link{XStringSet-class},
  \link{XStringViews-class},
  \link{MaskedXString-class},
  \code{\link{reverse}},
  \code{\link{rev}},
  \code{\link{strsplit}},
  \code{\link{GENETIC_CODE}},
  \code{\link{AMINO_ACID_CODE}}
}

\examples{
  data(yeastSEQCHR1)
  yeast1 <- DNAString(yeastSEQCHR1)

  alphabetFrequency(yeast1)
  alphabetFrequency(yeast1, baseOnly=TRUE)

  dinucleotideFrequency(yeast1)
  trinucleotideFrequency(yeast1)
  oligonucleotideFrequency(yeast1, 4)

  ## With a multiple sequence input
  library(drosophila2probe)
  x <- DNAStringSet(drosophila2probe$sequence)
  alphabetFrequency(x[1:50], baseOnly=TRUE)
  alphabetFrequency(x, baseOnly=TRUE, collapse=TRUE)

  ## Get the less and most represented 6-mers
  f6 <- oligonucleotideFrequency(yeast1, 6)
  f6[f6 == min(f6)]
  f6[f6 == max(f6)]

  ## Get the result as an array
  tri <- trinucleotideFrequency(yeast1, as.array=TRUE)
  tri["A", "A", "C"] # == trinucleotideFrequency(yeast1)["AAC"]
  tri["T", , ] # frequencies of trinucleotides starting with a "T"

  ## Get nucleotide transition matrices for yeast1
  oligonucleotideTransitions(yeast1)
  oligonucleotideTransitions(yeast1, 2, freq=TRUE)

  ## Note that when dropping the dimensions of the 'tri' array, elements
  ## in the resulting vector are ordered as if they were obtained with
  ## 'fast.moving.side="left"':
  triL <- trinucleotideFrequency(yeast1, fast.moving.side="left")
  all(as.vector(tri) == triL) # TRUE

  ## Convert the trinucleotide frequency into the amino acid frequency based on
  ## translation
  tri1 <- trinucleotideFrequency(yeast1)
  names(tri1) <- GENETIC_CODE[names(tri1)]
  sapply(split(tri1, names(tri1)), sum) # 12512 occurences of the stop codon

  ## When the returned vector is very long (e.g. width >= 10), using
  ## 'with.labels=FALSE' will improve the performance considerably (100x, 1000x
  ## or more):
  f12 <- oligonucleotideFrequency(yeast1, 12, with.labels=FALSE) # very fast!

  ## Some related utility functions
  dict1 <- mkAllStrings(LETTERS[1:3], 4)
  dict2 <- mkAllStrings(LETTERS[1:3], 4, fast.moving.side="left")
  identical(strrev(dict1), dict2) # TRUE 
}

\keyword{category}

