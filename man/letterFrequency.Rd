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
\alias{hasOnlyBaseLetters,MaskedXString-method}

\alias{uniqueLetters}
\alias{uniqueLetters,XString-method}
\alias{uniqueLetters,XStringSet-method}
\alias{uniqueLetters,XStringViews-method}
\alias{uniqueLetters,MaskedXString-method}

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

% Old stuff:
\alias{consmat}
\alias{consmat,ANY-method}


\title{Calculate the frequency of letters in a biological
sequence, or the consensus matrix of a set of sequences}

\description{
  Given a biological sequence (or a set of biological sequences),
  the \code{alphabetFrequency} function computes the frequency of each
  letter in the (base) alphabet.

  The \code{consensusMatrix} function computes the consensus matrix
  of a set of sequences, and the \code{consensusString} function creates
  the consensus sequence based on a 50\% + 1 vote from the consensus
  matrix (using the \code{"?"} letter to represent the lack of consensus).

  In this man page we call "DNA input" (or "RNA input") an
  \link{XString}, \link{XStringSet}, \link{XStringViews} or
  \link{MaskedXString} object of base type DNA (or RNA).
}

\usage{
  alphabetFrequency(x, as.prob=FALSE, freq=FALSE, ...)
  hasOnlyBaseLetters(x)
  uniqueLetters(x)

  consensusMatrix(x, as.prob=FALSE, freq=FALSE, shift=0L, width=NULL, ...)

  \S4method{consensusString}{matrix}(x)
  \S4method{consensusString}{XStringSet}(x, shift=0L, width=NULL)
  \S4method{consensusString}{ANY}(x)
}

\arguments{
  \item{x}{
    An \link{XString}, \link{XStringSet}, \link{XStringViews}
    or \link{MaskedXString} object for \code{alphabetFrequency}
    and \code{uniqueLetters}.

    DNA or RNA input for \code{hasOnlyBaseLetters}.

    A character vector, or an \link{XStringSet} or \link{XStringViews}
    object for \code{consensusMatrix}.

    A consensus matrix (as returned by \code{consensusMatrix}),
    or an \link{XStringSet} or \link{XStringViews} object
    for \code{consensusString}.
  }
  \item{as.prob}{
    If \code{TRUE} then probabilities are reported,
    otherwise counts (the default).
  }
  \item{freq}{
    This argument is deprecated.
    Please use the \code{as.prob} argument instead.
  }
  \item{...}{
    Further arguments to be passed to or from other methods.

    For the \link{XStringViews} and \link{XStringSet} methods,
    the \code{collapse} argument is accepted.

    For DNA or RNA input, the \code{baseOnly} argument is accepted.
    If \code{baseOnly} is \code{TRUE}, the returned vector (or matrix)
    only contains the frequencies of the letters that belong to the
    "base" alphabet of \code{x} i.e. to the alphabet returned by
    \code{alphabet(x, baseOnly=TRUE)}.
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
  \item{width}{
    The number of columns of the returned matrix for the \code{consensusMatrix}
    method for \link{XStringSet} objects.
    When \code{width=NULL} (the default), then this method returns a matrix
    that has just enough columns to have its last column aligned
    with the rightmost letter of all the sequences in \code{x} after
    those sequences have been shifted (see the \code{shift} argument above).
    This ensures that any wider consensus matrix would be a "padded with zeros"
    version of the matrix returned when \code{width=NULL}.
    
    The length of the returned sequence for the \code{consensusString}
    method for \link{XStringSet} objects.
  }
}

\details{
  \code{alphabetFrequency} is a generic function defined in the
  Biostrings package.
}

\value{
  \code{alphabetFrequency} returns a numeric vector when \code{x} is an
  \link{XString} or \link{MaskedXString} object. When \code{x} is an
  \link{XStringSet} or \link{XStringViews} object, then it returns
  a numeric matrix with \code{length(x)} rows where the
  \code{i}-th row contains the frequencies for \code{x[[i]]}.
  If \code{x} is a DNA or RNA input, then the returned vector is named
  with the letters in the alphabet. If the \code{baseOnly} argument is
  \code{TRUE}, then the returned vector has only 5 elements: 4 elements
  corresponding to the 4 nucleotides + the 'other' element.

  \code{hasOnlyBaseLetters} returns \code{TRUE} or \code{FALSE} indicating
  whether or not \code{x} contains only base letters (i.e. As, Cs, Gs and Ts
  for DNA input and As, Cs, Gs and Us for RNA input).

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
  \code{\link{oligonucleotideFrequency}},
  \code{\link{countPDict}},
  \link{XString-class},
  \link{XStringSet-class},
  \link{XStringViews-class},
  \link{MaskedXString-class},
  \code{\link{strsplit}}
}

\examples{
  ## ---------------------------------------------------------------------
  ## A. BASIC alphabetFrequency() EXAMPLES
  ## ---------------------------------------------------------------------
  data(yeastSEQCHR1)
  yeast1 <- DNAString(yeastSEQCHR1)

  alphabetFrequency(yeast1)
  alphabetFrequency(yeast1, baseOnly=TRUE)
  hasOnlyBaseLetters(yeast1)
  uniqueLetters(yeast1)

  ## With input made of multiple sequences:
  library(drosophila2probe)
  probes <- DNAStringSet(drosophila2probe$sequence)
  alphabetFrequency(probes[1:50], baseOnly=TRUE)
  alphabetFrequency(probes, baseOnly=TRUE, collapse=TRUE)

  ## ---------------------------------------------------------------------
  ## B. consensus*() EXAMPLES
  ## ---------------------------------------------------------------------
  ## Read in ORF data:
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
  consensusMatrix(probes, baseOnly=TRUE)

  ## After sorting, the first 5 probes might look similar (at least on
  ## their first bases):
  consensusString(sort(probes)[1:5])

  ## ---------------------------------------------------------------------
  ## C. RELATIONSHIP BETWEEN consensusMatrix() AND coverage()
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
