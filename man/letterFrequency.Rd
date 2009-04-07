\name{letterFrequency}

\alias{letterFrequency}

\alias{alphabetFrequency}
\alias{alphabetFrequency,XString-method}
\alias{alphabetFrequency,DNAString-method}
\alias{alphabetFrequency,RNAString-method}
\alias{alphabetFrequency,XStringSet-method}
\alias{alphabetFrequency,DNAStringSet-method}
\alias{alphabetFrequency,RNAStringSet-method}
\alias{alphabetFrequency,XStringViews-method}
\alias{alphabetFrequency,MaskedXString-method}

\alias{hasOnlyBaseLetters}
\alias{hasOnlyBaseLetters,DNAString-method}
\alias{hasOnlyBaseLetters,DNAStringSet-method}
\alias{hasOnlyBaseLetters,RNAString-method}
\alias{hasOnlyBaseLetters,RNAStringSet-method}
\alias{hasOnlyBaseLetters,XStringViews-method}
\alias{hasOnlyBaseLetters,MaskedDNAString-method}
\alias{hasOnlyBaseLetters,MaskedRNAString-method}

\alias{uniqueLetters}
\alias{uniqueLetters,XString-method}
\alias{uniqueLetters,XStringSet-method}
\alias{uniqueLetters,XStringViews-method}
\alias{uniqueLetters,MaskedXString-method}

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

\alias{consensusMatrix}
\alias{consensusMatrix,character-method}
\alias{consensusMatrix,matrix-method}
\alias{consensusMatrix,list-method}
\alias{consensusMatrix,XStringSet-method}
\alias{consensusMatrix,XStringViews-method}
\alias{consensusString}
\alias{consensusString,matrix-method}
\alias{consensusString,XStringSet-method}
\alias{consensusString,XStringViews-method}
\alias{consensusString,ANY-method}

% Deprecated:
\alias{consmat}
\alias{consmat,ANY-method}


\title{Functions to calculate the frequency of letters in a biological
sequence and other related functions}

\description{
  Given a biological sequence (or a set of biological sequences),
  the \code{alphabetFrequency} function computes the frequency of each
  letter in the (base) alphabet,
  the \code{dinucleotideFrequency} function the frequency of all
  possible dinucleotides,
  and the \code{trinucleotideFrequency} function the frequency
  of all possible trinucleotides.

  More generally, the \code{oligonucleotideFrequency} function
  computes the frequency of all possible oligonucleotides of a given
  length (called the "width" in this particular context).

  The \code{consensusMatrix} function computes the consensus matrix
  of a set of sequences, and the \code{consensusString} function creates
  the consensus sequence based on a 50\% + 1 vote from the consensus
  matrix (using the \code{"?"} letter to represent the lack of consensus).

  In this man page we call "DNA input" a \link{DNAString} object,
  or a \link{DNAStringSet} object, or an \link{XStringViews} object
  with a \link{DNAString} subject, or a \link{MaskedDNAString} object.
  Similarly we call "RNA input" an \link{RNAString} object,
  or an \link{RNAStringSet} object, or an \link{XStringViews} object
  with an \link{RNAString} subject, or a \link{MaskedRNAString} object.
}

\usage{
  alphabetFrequency(x, baseOnly=FALSE, freq=FALSE, ...)
  hasOnlyBaseLetters(x)
  uniqueLetters(x)

  dinucleotideFrequency(x, freq=FALSE, fast.moving.side="right",
                        as.matrix=FALSE, with.labels=TRUE, ...)
  trinucleotideFrequency(x, freq=FALSE, fast.moving.side="right",
                         as.array=FALSE, with.labels=TRUE, ...)
  oligonucleotideFrequency(x, width, freq=FALSE, fast.moving.side="right",
                           as.array=FALSE, with.labels=TRUE, ...)
  oligonucleotideTransitions(x, left=1, right=1, freq=FALSE)

  ## Some related utility functions
  strrev(x)
  mkAllStrings(alphabet, width, fast.moving.side="right")

  \S4method{consensusMatrix}{character}(x, freq=FALSE)
  \S4method{consensusMatrix}{XStringSet}(x, baseOnly=FALSE, freq=FALSE, shift=0L, width=NULL)

  \S4method{consensusString}{matrix}(x)
  \S4method{consensusString}{XStringSet}(x, shift=0L, width=NULL)
  \S4method{consensusString}{ANY}(x)
}

\arguments{
  \item{x}{
    An \link{XString}, \link{XStringSet}, \link{XStringViews}
    or \link{MaskedXString} object for the \code{*Frequency} and
    \code{uniqueLetters} functions.

    "DNA or RNA input" for \code{hasOnlyBaseLetters}.

    A character vector for \code{strrev}.

    A character vector, \link{XStringSet} or \link{XStringViews}
    object for \code{consensusMatrix}.

    A consensus matrix (as returned by \code{consensusMatrix}),
    or an \link{XStringSet} or \link{XStringViews} object
    for \code{consensusString}.
  }
  \item{baseOnly}{
    \code{TRUE} or \code{FALSE}.
    If \code{TRUE}, the returned vector (or matrix) only contains the
    frequencies of the letters that belong to the "base" alphabet
    of \code{x} i.e. to the alphabet returned by
    \code{alphabet(x, baseOnly=TRUE)}.
    Note that when \code{x} is not a "DNA or RNA input", then
    specifying \code{baseOnly} has no effect.
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
    The number of columns of the returned matrix for the \code{consensusMatrix}
    method for \link{XStringSet} objects.
    The length of the returned sequence for the \code{consensusString}
    method for \link{XStringSet} objects.
  }
  \item{left, right}{
    The number of nucleotides per oligonucleotide for the rows and columns
    respectively in the transition matrix created by \code{oligonucleotideTransitions}.
  }
  \item{alphabet}{
    The alphabet to use to make the strings.
  }
  \item{shift}{
    An integer vector (recycled to the length of \code{x}) specifying how
    each sequence in \code{x} should be (horizontally) shifted with respect
    to the first column of the consensus matrix to be returned.
    By default (\code{shift=0}), each sequence in \code{x} has its
    first letter aligned with the first column of the matrix.
    A positive \code{shift} value means that the corresponding sequence
    must be shifted to the right, and a negative \code{shift} value
    that it must be shifted to the left.
    For example, a shift of 5 means that it must be shifted 5 positions
    to the right (i.e. the first letter in the sequence must be aligned
    with the 6th column of the matrix), and a shift of -3 means that
    it must be shifted 3 positions to the left (i.e. the 4th letter in
    the sequence must be aligned with the first column of the matrix).
  }
}

\details{
  \code{alphabetFrequency}  and \code{oligonucleotideFrequency} are generic
  functions defined in the Biostrings package with methods defined for
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

  \code{hasOnlyBaseLetters} returns \code{TRUE} or \code{FALSE} indicating
  whether or not \code{x} contains only base letters (i.e. As, Cs, Gs and Ts
  for "DNA input" and As, Cs, Gs and Us for "RNA input").

  \code{uniqueLetters} returns a vector of 1-letter or empty strings. The empty
  string is used to represent the nul character if \code{x} happens to contain
  any. Note that this can only happen if the base class of \code{x}
  is \link{BString}.

  An integer matrix with letters as row names for \code{consensusMatrix}.

  A standard character string for \code{consensusString}.
}

\author{H. Pages and P. Aboyoun}

\seealso{
  \code{\link{alphabet}},
  \code{\link[IRanges]{coverage}},
  \link{XStringSet-class},
  \link{XStringViews-class},
  \link{MaskedXString-class},
  \code{\link{reverse,XString-method}},
  \code{\link{rev}},
  \code{\link{strsplit}},
  \code{\link{GENETIC_CODE}},
  \code{\link{AMINO_ACID_CODE}}
  \code{\link{countPDict}}
}

\examples{
  ## ---------------------------------------------------------------------
  ## A. BASIC *Frequency() EXAMPLES
  ## ---------------------------------------------------------------------
  data(yeastSEQCHR1)
  yeast1 <- DNAString(yeastSEQCHR1)

  alphabetFrequency(yeast1)
  alphabetFrequency(yeast1, baseOnly=TRUE)
  hasOnlyBaseLetters(yeast1)
  uniqueLetters(yeast1)

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

  ## ---------------------------------------------------------------------
  ## D. consensus*() EXAMPLES
  ## ---------------------------------------------------------------------
  ## Read in ORF data
  file <- system.file("extdata", "someORF.fa", package="Biostrings")
  orf <- read.DNAStringSet(file, "fasta")

  ## To illustrate, the following example assumes the ORF data
  ## to be aligned for the first 10 positions (patently false):
  orf10 <- DNAStringSet(orf, end=10)
  consensusMatrix(orf10, baseOnly=TRUE)

  ## The following example assumes the first 10 positions to be aligned
  ## after some incremental shifting to the right (patently false):
  consensusMatrix(orf10, baseOnly=TRUE, shift=0:6)
  consensusMatrix(orf10, baseOnly=TRUE, shift=0:6, width=10)

  ## For the character matrix containing the "exploded" representation
  ## of the strings, do:
  as.matrix(orf10, use.names=FALSE)

  ## consensusMatrix() can be used to just compute the alphabet frequency
  ## for each position in the input sequences:
  library(drosophila2probe)
  probes <- DNAStringSet(drosophila2probe$sequence)
  consensusMatrix(probes, baseOnly=TRUE)

  ## After sorting, the first 5 probes might look similar (at least on
  ## their first bases):
  consensusString(sort(probes)[1:5])

  ## ---------------------------------------------------------------------
  ## E. RELATIONSHIP BETWEEN consensusMatrix() AND coverage()
  ## ---------------------------------------------------------------------
  ## Applying colSums() on a consensus matrix gives the coverage that
  ## would be obtained by piling up (after shifting) the input sequences
  ## on top of an (imaginary) reference sequence:
  cm <- consensusMatrix(orf10, shift=0:6, width=10)
  colSums(cm)

  ## Note that this coverage can also be obtained with:
  as.integer(coverage(IRanges(rep(1, length(orf)), width(orf)), shift=0:6, width=10))
}

\keyword{methods}
\keyword{manip}
