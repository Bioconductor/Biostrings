\name{oligonucleotideFrequency}

\alias{oligonucleotideFrequency}

\alias{oligonucleotideFrequency}
\alias{oligonucleotideFrequency,XString-method}
\alias{oligonucleotideFrequency,XStringSet-method}
\alias{oligonucleotideFrequency,XStringViews-method}
\alias{oligonucleotideFrequency,MaskedXString-method}

\alias{dinucleotideFrequency}
\alias{trinucleotideFrequency}

\alias{oligonucleotideTransitions}
\alias{mkAllStrings}


\title{Functions to calculate the frequency of letters in a biological
sequence and other related functions}

\description{
  Given a biological sequence (or a set of biological sequences),
  the \code{oligonucleotideFrequency} function computes the frequency
  of all possible oligonucleotides of a given length (called the "width"
  in this particular context).

  The \code{dinucleotideFrequency} and \code{trinucleotideFrequency}
  functions are convenient wrappers for calling \code{oligonucleotideFrequency}
  with \code{width=2} and \code{width=3}, respectively.
  
  In this man page we call "DNA input" (or "RNA input") an
  \link{XString}, \link{XStringSet}, \link{XStringViews} or
  \link{MaskedXString} object of base type DNA (or RNA).
}

\usage{
  oligonucleotideFrequency(x, width, freq=FALSE, fast.moving.side="right",
                           as.array=FALSE, with.labels=TRUE, ...)

  \S4method{oligonucleotideFrequency}{XStringSet}(x,
      width, freq=FALSE, fast.moving.side="right",
      as.array=FALSE, with.labels=TRUE, simplify.as="matrix")

  dinucleotideFrequency(x, freq=FALSE, fast.moving.side="right",
                        as.matrix=FALSE, with.labels=TRUE, ...)
  trinucleotideFrequency(x, freq=FALSE, fast.moving.side="right",
                         as.array=FALSE, with.labels=TRUE, ...)

  ## Related utility functions:
  oligonucleotideTransitions(x, left=1, right=1, freq=FALSE)
  mkAllStrings(alphabet, width, fast.moving.side="right")
}

\arguments{
  \item{x}{
    DNA or RNA input.
  }
  \item{width}{
    The number of nucleotides per oligonucleotide for
    \code{oligonucleotideFrequency}.

    The number of letters per string for \code{mkAllStrings}.
  }
  \item{freq}{
    If \code{TRUE} then relative frequencies are reported,
    otherwise counts (the default).
  }
  \item{fast.moving.side}{
    Which side of the strings should move fastest?
  }
  \item{as.array,as.matrix}{
    Controls the "shape" of the returned object.
    If \code{TRUE} then it's a numeric matrix (or array),
    otherwise it's just a "flat" numeric vector i.e. a
    vector with no dim attribute (the default).
  }
  \item{with.labels}{
    If \code{TRUE} then the returned object is named.
  }
  \item{...}{
    Further arguments to be passed to or from other methods.
  }
  \item{simplify.as}{
    Together with the \code{as.array} and \code{as.matrix}
    arguments, controls the "shape" of the returned object
    when the input \code{x} is an \link{XStringSet} or
    \link{XStringViews} object.
    Supported \code{simplify.as} values are \code{"matrix"}
    (the default), \code{"list"} and \code{"collapsed"}.
    If \code{simplify.as} is \code{"matrix"}, the returned
    object is a matrix with \code{length(x)} rows where the
    \code{i}-th row contains the frequencies for \code{x[[i]]}.
    If \code{simplify.as} is \code{"list"}, the returned
    object is a list of the same length as \code{length(x)}
    where the \code{i}-th element contains the frequencies
    for \code{x[[i]]}.
    If \code{simplify.as} is \code{"collapsed"}, then the
    the frequencies are computed for the entire object \code{x}
    as a whole (i.e. frequencies cumulated across all sequences
    in \code{x}).
  }
  \item{left, right}{
    The number of nucleotides per oligonucleotide for the rows
    and columns respectively in the transition matrix created
    by \code{oligonucleotideTransitions}.
  }
  \item{alphabet}{
    The alphabet to use to make the strings.
  }
}

\details{
  \code{oligonucleotideFrequency} is a generic function defined in the
  Biostrings package with methods for \link{XString}, \link{XStringSet},
  \link{XStringViews} and \link{MaskedXString} objects.
}

\value{
  If \code{x} is an \link{XString} or \link{MaskedXString} object,
  the \code{*Frequency} functions return a numeric vector of length
  \code{4^width}. If \code{as.array} (or \code{as.matrix}) is \code{TRUE},
  then this vector is formatted as an array (or matrix).
  If \code{x} is an \link{XStringSet} or \link{XStringViews} object,
  the returned object has the shape specified by the \code{simplify.as}
  argument.
}

\author{H. Pages and P. Aboyoun}

\seealso{
  \code{\link{alphabetFrequency}},
  \code{\link{alphabet}},
  \link{XString-class},
  \link{XStringSet-class},
  \link{XStringViews-class},
  \link{MaskedXString-class},
  \code{\link{GENETIC_CODE}},
  \code{\link{AMINO_ACID_CODE}}
}

\examples{
  ## ---------------------------------------------------------------------
  ## A. BASIC *Frequency() EXAMPLES
  ## ---------------------------------------------------------------------
  data(yeastSEQCHR1)
  yeast1 <- DNAString(yeastSEQCHR1)

  dinucleotideFrequency(yeast1)
  trinucleotideFrequency(yeast1)
  oligonucleotideFrequency(yeast1, 4)

  ## With input made of multiple sequences:
  library(drosophila2probe)
  x <- DNAStringSet(drosophila2probe$sequence)

  ## Get the less and most represented 6-mers:
  f6 <- oligonucleotideFrequency(yeast1, 6)
  f6[f6 == min(f6)]
  f6[f6 == max(f6)]

  ## Get the result as an array
  tri <- trinucleotideFrequency(yeast1, as.array=TRUE)
  tri["A", "A", "C"] # == trinucleotideFrequency(yeast1)["AAC"]
  tri["T", , ] # frequencies of trinucleotides starting with a "T"

  ## ---------------------------------------------------------------------
  ## B. oligonucleotideTransitions()
  ## ---------------------------------------------------------------------
  ## Get nucleotide transition matrices for yeast1
  oligonucleotideTransitions(yeast1)
  oligonucleotideTransitions(yeast1, 2, freq=TRUE)

  ## ---------------------------------------------------------------------
  ## C. ADVANCED *Frequency() EXAMPLES
  ## ---------------------------------------------------------------------
  ## Note that when dropping the dimensions of the 'tri' array, elements
  ## in the resulting vector are ordered as if they were obtained with
  ## 'fast.moving.side="left"':
  triL <- trinucleotideFrequency(yeast1, fast.moving.side="left")
  all(as.vector(tri) == triL) # TRUE

  ## Convert the trinucleotide frequency into the amino acid frequency based on
  ## translation
  tri1 <- trinucleotideFrequency(yeast1)
  names(tri1) <- GENETIC_CODE[names(tri1)]
  sapply(split(tri1, names(tri1)), sum) # 12512 occurrences of the stop codon

  ## When the returned vector is very long (e.g. width >= 10), using
  ## 'with.labels=FALSE' will improve the performance considerably (100x, 1000x
  ## or more):
  f12 <- oligonucleotideFrequency(yeast1, 12, with.labels=FALSE) # very fast!

  ## Some related utility functions
  dict1 <- mkAllStrings(LETTERS[1:3], 4)
  dict2 <- mkAllStrings(LETTERS[1:3], 4, fast.moving.side="left")
  identical(strrev(dict1), dict2) # TRUE 
}

\keyword{methods}
\keyword{manip}
