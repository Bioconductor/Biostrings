\name{nucleotideFrequency}

\alias{oligonucleotideFrequency}

\alias{oligonucleotideFrequency}
\alias{oligonucleotideFrequency,XString-method}
\alias{oligonucleotideFrequency,XStringSet-method}
\alias{oligonucleotideFrequency,XStringViews-method}
\alias{oligonucleotideFrequency,MaskedXString-method}

\alias{dinucleotideFrequency}
\alias{trinucleotideFrequency}

\alias{nucleotideFrequencyAt}
\alias{nucleotideFrequencyAt,XStringSet-method}
\alias{nucleotideFrequencyAt,XStringViews-method}

\alias{oligonucleotideTransitions}
\alias{mkAllStrings}


\title{Calculate the frequency of oligonucleotides in a DNA
  or RNA sequence (and other related functions)}

\description{
  Given a DNA or RNA sequence (or a set of DNA or RNA sequences),
  the \code{oligonucleotideFrequency} function computes the frequency
  of all possible oligonucleotides of a given length (called the "width"
  in this particular context).

  The \code{dinucleotideFrequency} and \code{trinucleotideFrequency}
  functions are convenient wrappers for calling \code{oligonucleotideFrequency}
  with \code{width=2} and \code{width=3}, respectively.

  The \code{nucleotideFrequencyAt} function computes the frequency
  of the short sequences formed by extracting the nucleotides found
  at some fixed positions from each sequence of a set of DNA or RNA
  sequences.
 
  In this man page we call "DNA input" (or "RNA input") an
  \link{XString}, \link{XStringSet}, \link{XStringViews} or
  \link{MaskedXString} object of base type DNA (or RNA).
}

\usage{
oligonucleotideFrequency(x, width, as.prob=FALSE, as.array=FALSE,
                         fast.moving.side="right", with.labels=TRUE, ...)

\S4method{oligonucleotideFrequency}{XStringSet}(x,
    width, as.prob=FALSE, as.array=FALSE,
    fast.moving.side="right", with.labels=TRUE, simplify.as="matrix")

dinucleotideFrequency(x, as.prob=FALSE, as.matrix=FALSE,
                      fast.moving.side="right", with.labels=TRUE, ...)
trinucleotideFrequency(x, as.prob=FALSE, as.array=FALSE,
                       fast.moving.side="right", with.labels=TRUE, ...)

nucleotideFrequencyAt(x, at, as.prob=FALSE, as.array=TRUE,
                      fast.moving.side="right", with.labels=TRUE, ...)

## Some related functions:
oligonucleotideTransitions(x, left=1, right=1, as.prob=FALSE)
mkAllStrings(alphabet, width, fast.moving.side="right")
}

\arguments{
  \item{x}{
    Any DNA or RNA input for the \code{*Frequency} and
    \code{oligonucleotideTransitions} functions.

    An \link{XStringSet} or \link{XStringViews} object of base type DNA or RNA
    for \code{nucleotideFrequencyAt}.
  }
  \item{width}{
    The number of nucleotides per oligonucleotide for
    \code{oligonucleotideFrequency}.

    The number of letters per string for \code{mkAllStrings}.
  }
  \item{at}{
    An integer vector containing the positions to look at in each element
    of \code{x}.
  }
  \item{as.prob}{
    If \code{TRUE} then probabilities are reported,
    otherwise counts (the default).
  }
  \item{as.array,as.matrix}{
    Controls the "shape" of the returned object.
    If \code{TRUE} (the default for \code{nucleotideFrequencyAt})
    then it's a numeric matrix (or array),
    otherwise it's just a "flat" numeric vector i.e. a
    vector with no dim attribute (the default for the
    \code{*Frequency} functions).
  }
  \item{fast.moving.side}{
    Which side of the strings should move fastest?
    Note that, when \code{as.array} is TRUE, then the supplied value
    is ignored and the effective value is \code{"left"}.
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
  \code{\link{hasLetterAt}},
  \link{XString-class},
  \link{XStringSet-class},
  \link{XStringViews-class},
  \link{MaskedXString-class},
  \code{\link{GENETIC_CODE}},
  \code{\link{AMINO_ACID_CODE}},
  \code{\link{reverse,XString-method}},
  \code{\link{rev}}
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

  ## Get the less and most represented 6-mers:
  f6 <- oligonucleotideFrequency(yeast1, 6)
  f6[f6 == min(f6)]
  f6[f6 == max(f6)]

  ## Get the result as an array:
  tri <- trinucleotideFrequency(yeast1, as.array=TRUE)
  tri["A", "A", "C"] # == trinucleotideFrequency(yeast1)["AAC"]
  tri["T", , ] # frequencies of trinucleotides starting with a "T"

  ## With input made of multiple sequences:
  library(drosophila2probe)
  probes <- DNAStringSet(drosophila2probe)
  dfmat <- dinucleotideFrequency(probes)  # a big matrix
  dinucleotideFrequency(probes, simplify.as="collapsed")
  dinucleotideFrequency(probes, simplify.as="collapsed", as.matrix=TRUE)

  ## ---------------------------------------------------------------------
  ## B. OBSERVED DINUCLEOTIDE FREQUENCY VERSUS EXPECTED DINUCLEOTIDE
  ##    FREQUENCY
  ## ---------------------------------------------------------------------
  ## The expected frequency of dinucleotide "ab" based on the frequencies
  ## of its individual letters "a" and "b" is:
  ##    exp_Fab = Fa * Fb / N if the 2 letters are different (e.g. CG)
  ##    exp_Faa = Fa * (Fa-1) / N if the 2 letters are the same (e.g. TT)
  ## where Fa and Fb are the frequencies of "a" and "b" (respectively) and
  ## N the length of the sequence.
  
  ## Here is a simple function that implements the above formula for a
  ## DNAString object 'x'. The expected frequencies are returned in a 4x4
  ## matrix where the rownames and colnames correspond to the 1st and 2nd
  ## base in the dinucleotide:
  expectedDinucleotideFrequency <- function(x)
  {
      # Individual base frequencies.
      bf <- alphabetFrequency(x, baseOnly=TRUE)[DNA_BASES]
      (as.matrix(bf) \%*\% t(bf) - diag(bf)) / length(x)
  }

  ## On Celegans chrI:
  library(BSgenome.Celegans.UCSC.ce2)
  chrI <- Celegans$chrI
  obs_df <- dinucleotideFrequency(chrI, as.matrix=TRUE)
  obs_df  # CG has the lowest frequency
  exp_df <- expectedDinucleotideFrequency(chrI)
  ## A sanity check:
  stopifnot(as.integer(sum(exp_df)) == sum(obs_df))

  ## Ratio of observed frequency to expected frequency:
  obs_df / exp_df  # TA has the lowest ratio, not CG!

  ## ---------------------------------------------------------------------
  ## C. nucleotideFrequencyAt()
  ## ---------------------------------------------------------------------
  nucleotideFrequencyAt(probes, 13)
  nucleotideFrequencyAt(probes, c(13, 20))
  nucleotideFrequencyAt(probes, c(13, 20), as.array=FALSE)

  ## nucleotideFrequencyAt() can be used to answer questions like: "how
  ## many probes in the drosophila2 chip have T, G, T, A at position
  ## 2, 4, 13 and 20, respectively?"
  nucleotideFrequencyAt(probes, c(2, 4, 13, 20))["T", "G", "T", "A"]
  ## or "what's the probability to have an A at position 25 if there is
  ## one at position 13?"
  nf <- nucleotideFrequencyAt(probes, c(13, 25))
  sum(nf["A", "A"]) / sum(nf["A", ])
  ## Probabilities to have other bases at position 25 if there is an A
  ## at position 13:
  sum(nf["A", "C"]) / sum(nf["A", ])  # C
  sum(nf["A", "G"]) / sum(nf["A", ])  # G
  sum(nf["A", "T"]) / sum(nf["A", ])  # T

  ## See ?hasLetterAt for another way to get those results.

  ## ---------------------------------------------------------------------
  ## D. oligonucleotideTransitions()
  ## ---------------------------------------------------------------------
  ## Get nucleotide transition matrices for yeast1
  oligonucleotideTransitions(yeast1)
  oligonucleotideTransitions(yeast1, 2, as.prob=TRUE)

  ## ---------------------------------------------------------------------
  ## E. ADVANCED *Frequency() EXAMPLES
  ## ---------------------------------------------------------------------
  ## Note that when dropping the dimensions of the 'tri' array, elements
  ## in the resulting vector are ordered as if they were obtained with
  ## 'fast.moving.side="left"':
  triL <- trinucleotideFrequency(yeast1, fast.moving.side="left")
  all(as.vector(tri) == triL) # TRUE

  ## Convert the trinucleotide frequency into the amino acid frequency
  ## based on translation:
  tri1 <- trinucleotideFrequency(yeast1)
  names(tri1) <- GENETIC_CODE[names(tri1)]
  sapply(split(tri1, names(tri1)), sum) # 12512 occurrences of the stop codon

  ## When the returned vector is very long (e.g. width >= 10), using
  ## 'with.labels=FALSE' can improve performance significantly.
  ## Here for example, the observed speed up is between 25x and 500x:
  f12 <- oligonucleotideFrequency(yeast1, 12, with.labels=FALSE) # very fast!

  ## Spome related functions:
  dict1 <- mkAllStrings(LETTERS[1:3], 4)
  dict2 <- mkAllStrings(LETTERS[1:3], 4, fast.moving.side="left")
  stopifnot(identical(reverse(dict1), dict2)) 
}

\keyword{methods}
\keyword{manip}
