\name{letterFrequency}

\alias{letterFrequency}
\alias{letterFrequency,MaskedXString-method}
\alias{letterFrequency,XString-method}
\alias{letterFrequency,XStringSet-method}
\alias{letterFrequency,XStringViews-method}

\alias{alphabetFrequency}
\alias{alphabetFrequency,XString-method}
\alias{alphabetFrequency,DNAString-method}
\alias{alphabetFrequency,RNAString-method}
\alias{alphabetFrequency,AAString-method}
\alias{alphabetFrequency,XStringSet-method}
\alias{alphabetFrequency,DNAStringSet-method}
\alias{alphabetFrequency,RNAStringSet-method}
\alias{alphabetFrequency,AAStringSet-method}
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

\alias{letterFrequencyInSlidingView}
\alias{letterFrequencyInSlidingView,XString-method}
%\alias{letterFrequencyInSlidingView,XStringSet-method}
%\alias{letterFrequencyInSlidingView,XStringViews-method}
%\alias{letterFrequencyInSlidingView,MaskedXString-method}

\alias{consensusMatrix}
\alias{consensusMatrix,character-method}
\alias{consensusMatrix,matrix-method}
\alias{consensusMatrix,XStringSet-method}
\alias{consensusMatrix,XStringViews-method}
\alias{consensusString}
\alias{consensusString,matrix-method}
\alias{consensusString,BStringSet-method}
\alias{consensusString,DNAStringSet-method}
\alias{consensusString,RNAStringSet-method}
\alias{consensusString,XStringViews-method}
\alias{consensusString,ANY-method}


\title{Calculate the frequency of letters in a biological
sequence, or the consensus matrix of a set of sequences}

\description{
  Given a biological sequence (or a set of biological sequences),
  the \code{alphabetFrequency} function computes the frequency of
  each letter of the relevant \link{alphabet}.

  \code{letterFrequency} is similar, but more compact if one is only
  interested in certain letters.
  It can also tabulate letters "in common".

  \code{letterFrequencyInSlidingView} is a more specialized version
  of \code{letterFrequency} for (non-masked) \link{XString} objects.
  It tallys the requested letter frequencies for a fixed-width view,
  or window, that is conceptually slid along the entire input sequence.

  The \code{consensusMatrix} function computes the consensus matrix
  of a set of sequences, and the \code{consensusString} function creates
  the consensus sequence from the consensus matrix based upon specified
  criteria.

  In this man page we call "DNA input" (or "RNA input") an
  \link{XString}, \link{XStringSet}, \link{XStringViews} or
  \link{MaskedXString} object of base type DNA (or RNA).
}

\usage{
alphabetFrequency(x, as.prob=FALSE, ...)
hasOnlyBaseLetters(x)
uniqueLetters(x)

letterFrequency(x, letters, OR="|", as.prob=FALSE, ...)
letterFrequencyInSlidingView(x, view.width, letters, OR="|", as.prob=FALSE)

consensusMatrix(x, as.prob=FALSE, shift=0L, width=NULL, ...)

\S4method{consensusString}{matrix}(x, ambiguityMap="?", threshold=0.5)
\S4method{consensusString}{DNAStringSet}(x, ambiguityMap=IUPAC_CODE_MAP,
             threshold=0.25, shift=0L, width=NULL)
\S4method{consensusString}{RNAStringSet}(x,
             ambiguityMap=
             structure(as.character(RNAStringSet(DNAStringSet(IUPAC_CODE_MAP))),
                       names=
                       as.character(RNAStringSet(DNAStringSet(names(IUPAC_CODE_MAP))))),
             threshold=0.25, shift=0L, width=NULL)
}

\arguments{
  \item{x}{
    An \link{XString}, \link{XStringSet}, \link{XStringViews}
    or \link{MaskedXString} object for \code{alphabetFrequency},
    \code{letterFrequency}, or \code{uniqueLetters}.

    DNA or RNA input for \code{hasOnlyBaseLetters}.

    An \link{XString} object for \code{letterFrequencyInSlidingView}.

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
  \item{view.width}{
    For \code{letterFrequencyInSlidingView},
    the constant (e.g. 35, 48, 1000) size of the "window" to slide along
    \code{x}.
    The specified \code{letters} are tabulated in each window of length
    \code{view.width}.
    The rows of the result (see value) correspond to the various windows.
  }
  \item{letters}{
    For \code{letterFrequency} or \code{letterFrequencyInSlidingView},
    a character vector (e.g. "C", "CG", \link{c}("C", "G")) giving the
    letters to tabulate.
    When \code{x} is DNA or RNA input, \code{letters} must come from
    \code{\link{alphabet}(x)}.
    Except with \code{OR=0}, multi-character elements of letters
    ('nchar' > 1) are taken as groupings of letters into subsets, to
    be tabulated in common ("or"'d), as if their alphabetFrequency's
    were added (\link{Arithmetic}).
    The columns of the result (see value) correspond to the individual
    and sets of letters which are counted separately.
    Unrelated (and, with some post-processing, related) counts may of
    course be obtained in separate calls.
  }
  \item{OR}{
    For \code{letterFrequency} or \code{letterFrequencyInSlidingView},
    the string (default \code{|}) to use as a separator in forming names
    for the "grouped" columns, e.g. "C|G".
    The otherwise exceptional value \code{0} (zero) disables or'ing and
    is provided for convenience, allowing a single multi-character string
    (or several strings) of letters that should be counted separately.
    If some but not all letters are to be counted separately, they must
    reside in separate elements of letters (with 'nchar' 1 unless they
    are to be grouped with other letters), and \code{OR} cannot be 0.
  }
  \item{ambiguityMap}{
    Either a single character to use when agreement is not reached or
    a named character vector where the names are the ambiguity characters
    and the values are the combinations of letters that comprise the
    ambiguity (e.g. \code{link{IUPAC_CODE_MAP}}).
    When \code{ambiguityMap} is a named character vector, occurrences of
    ambiguous letters in \code{x} are replaced with their base alphabet
    letters that have been equally weighted to sum to 1. (See Details for
    some examples.)
  }
  \item{threshold}{
    The minimum probability threshold for an agreement to be declared.
    When \code{ambiguityMap} is a single character, \code{threshold}
    is a single number in (0, 1].
    When \code{ambiguityMap} is a named character vector
    (e.g. \code{link{IUPAC_CODE_MAP}}),  \code{threshold}
    is a single number in (0, 1/sum(nchar(ambiguityMap) == 1)].
  }
  \item{...}{
    Further arguments to be passed to or from other methods.

    For the \link{XStringViews} and \link{XStringSet} methods,
    the \code{collapse} argument is accepted.

    Except for \code{letterFrequency} or \code{letterFrequencyInSlidingView},
    and with DNA or RNA input, the \code{baseOnly} argument is accepted.
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
  \code{alphabetFrequency}, \code{letterFrequency}, and
  \code{letterFrequencyInSlidingView} are
  generic functions defined in the Biostrings package.

  \code{letterFrequency} is similar to \code{alphabetFrequency} but
  specific to the letters of interest, hence more compact, especially
  with \code{OR} non-zero.

  \code{letterFrequencyInSlidingView} yields the same result, on the
  sequence \code{x}, that \code{letterFrequency} would, if applied to the
  hypothetical (and possibly huge) \code{\link{XStringViews}} object
  consisting of all the intervals of length \code{view.width} on \code{x}.
  Taking advantage of the knowledge that successive "views" are nearly
  identical, for letter counting purposes, it is both lighter and faster.

  For \code{letterFrequencyInSlidingView}, a masked (\link{MaskedXString})
  object \code{x} is only supported through a cast to an (ordinary)
  \link{XString} such as \code{\link{unmasked}} (which includes its masked
  regions).

  When \code{consensusString} is executed with a named character
  \code{ambiguityMap} argument, it weights each input string equally and
  assigns an equal probability to each of the base letters represented by
  an ambiguity letter. So for DNA and a \code{threshold} of 0.25,
  a "G" and an "R" would result in an "R" since
  1/2 "G" + 1/2 "R" = 3/4 "G" + 1/4 "A" => "R";
  two "G"'s and one "R" would result in a "G" since
  2/3 "G" + 1/3 "R" = 5/6 "G" + 1/6 "A" => "G"; and
  one "A" and one "N" would result in an "N" since
  1/2 "A" + 1/2 "N" = 5/8 "A" + 1/8 "C" + 1/8 "G" + 1/8 "T" => "N".
}

\value{
  \code{alphabetFrequency} returns an integer vector when \code{x} is an
  \link{XString} or \link{MaskedXString} object. When \code{x} is an
  \link{XStringSet} or \link{XStringViews} object, then it returns
  an integer matrix with \code{length(x)} rows where the
  \code{i}-th row contains the frequencies for \code{x[[i]]}.
  If \code{x} is a DNA, RNA, or AA input, then the returned vector is named
  with the letters in the alphabet. If the \code{baseOnly} argument is
  \code{TRUE}, then the returned vector has only 5 elements for DNA/RNA input
  (4 elements corresponding to the 4 nucleotides + the 'other' element) and
  21 elements for AA input (20 elements corresponding to the 20 base amino acids
  + the 'other' element).

  \code{letterFrequency} returns, similarly, an integer vector or matrix,
  but restricted and/or collated according to \code{letters} and \code{OR}.

  \code{letterFrequencyInSlidingView} returns, for an \link{XString}
  object \code{x} of length (\code{\link{nchar}}) L, an integer matrix
  with \code{L-view.width+1} rows, the \code{i}-th of which holding the
  letter frequencies of \code{\link{substring}(x, i, i+view.width-1)}.

  \code{hasOnlyBaseLetters} returns \code{TRUE} or \code{FALSE} indicating
  whether or not \code{x} contains only base letters (i.e. As, Cs, Gs and Ts
  for DNA input, As, Cs, Gs and Us for RNA input, or any of the 20 standard
  amino acids for AA input).

  \code{uniqueLetters} returns a vector of 1-letter or empty strings. The empty
  string is used to represent the nul character if \code{x} happens to contain
  any. Note that this can only happen if the base class of \code{x}
  is \link{BString}.

  An integer matrix with letters as row names for \code{consensusMatrix}.

  A standard character string for \code{consensusString}.
}

\author{H. Pagès and P. Aboyoun; H. Jaffee for letterFrequency and
 letterFrequencyInSlidingView}

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
## alphabetFrequency()
## ---------------------------------------------------------------------
data(yeastSEQCHR1)
yeast1 <- DNAString(yeastSEQCHR1)

alphabetFrequency(yeast1)
alphabetFrequency(yeast1, baseOnly=TRUE)

hasOnlyBaseLetters(yeast1)
uniqueLetters(yeast1)

## With input made of multiple sequences:
library(drosophila2probe)
probes <- DNAStringSet(drosophila2probe)
alphabetFrequency(probes[1:50], baseOnly=TRUE)
alphabetFrequency(probes, baseOnly=TRUE, collapse=TRUE)

## ---------------------------------------------------------------------
## letterFrequency()
## ---------------------------------------------------------------------
letterFrequency(probes[[1]], letters="ACGT", OR=0)
base_letters <- alphabet(probes, baseOnly=TRUE)
base_letters
letterFrequency(probes[[1]], letters=base_letters, OR=0)
base_letter_freqs <- letterFrequency(probes, letters=base_letters, OR=0)
head(base_letter_freqs)
GC_content <- letterFrequency(probes, letters="CG")
head(GC_content)
letterFrequency(probes, letters="CG", collapse=TRUE)

## ---------------------------------------------------------------------
## letterFrequencyInSlidingView()
## ---------------------------------------------------------------------
data(yeastSEQCHR1)
x <- DNAString(yeastSEQCHR1)
view.width <- 48
letters <- c("A", "CG")
two_columns <- letterFrequencyInSlidingView(x, view.width, letters)
head(two_columns)
tail(two_columns)
three_columns <- letterFrequencyInSlidingView(x, view.width, letters, OR=0)
head(three_columns)
tail(three_columns)
stopifnot(identical(two_columns[ , "C|G"],
                    three_columns[ , "C"] + three_columns[ , "G"]))

## Note that, alternatively, 'three_columns' can also be obtained by
## creating the views on 'x' (as a Views object) and by calling
## alphabetFrequency() on it. But, of course, that is be *much* less
## efficient (both, in terms of memory and speed) than using
## letterFrequencyInSlidingView():
v <- Views(x, start=seq_len(length(x) - view.width + 1), width=view.width)
v
three_columns2 <- alphabetFrequency(v, baseOnly=TRUE)[ , c("A", "C", "G")]
stopifnot(identical(three_columns2, three_columns))

## Set the width of the view to length(x) to get the global frequencies:
letterFrequencyInSlidingView(x, letters="ACGTN", view.width=length(x), OR=0)

## ---------------------------------------------------------------------
## consensus*()
## ---------------------------------------------------------------------
## Read in ORF data:
file <- system.file("extdata", "someORF.fa", package="Biostrings")
orf <- readDNAStringSet(file)

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
consensusString(sort(probes)[1:5], ambiguityMap = "N", threshold = 0.5)

## Consensus involving ambiguity letters in the input strings
consensusString(DNAStringSet(c("NNNN","ACTG")))
consensusString(DNAStringSet(c("AANN","ACTG")))
consensusString(DNAStringSet(c("ACAG","ACAR")))
consensusString(DNAStringSet(c("ACAG","ACAR", "ACAG")))

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
