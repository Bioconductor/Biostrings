\name{reverseComplement}

\alias{reverse}
\alias{reverse,BString-method}
\alias{reverse,BStringSet-method}
\alias{reverse,BStringViews-method}
\alias{complement}
\alias{complement,DNAString-method}
\alias{complement,RNAString-method}
\alias{complement,DNAStringSet-method}
\alias{complement,RNAStringSet-method}
\alias{complement,BStringViews-method}
\alias{reverseComplement}
\alias{reverseComplement,DNAString-method}
\alias{reverseComplement,RNAString-method}
\alias{reverseComplement,DNAStringSet-method}
\alias{reverseComplement,RNAStringSet-method}
\alias{reverseComplement,BStringViews-method}


\title{Sequence reversing and complementing}

\description{
  These functions can reverse a \link{BString}, \link{DNAString}
  or \link{RNAString} object and complement each base
  of a \link{DNAString} object.
}

\usage{
  reverse(x, ...)
  complement(x, ...)
  reverseComplement(x, ...)
}

\arguments{
  \item{x}{
    A \link{BString} (or derived) object,
    a \link{BStringSet} (or derived) object,
    or a \link{BStringViews} object for \code{reverse}.

    A \link{DNAString}, \link{RNAString},
    \link{DNAStringSet} or \link{RNAStringSet} object,
    or a \link{BStringViews} object with a \link{DNAString}
    or \link{RNAString} subject for \code{complement}
    and \code{reverseComplement}.
  }
  \item{...}{
    Additional arguments to be passed to or from methods.
  }
}

\details{
  Given an object \code{x} of class \link{BString}, \link{DNAString}
  or \link{RNAString}, \code{reverse(x)} returns an object of the same class
  where letters in \code{x} are reordered in the reverse order.

  If \code{x} is a \link{DNAString} or \link{RNAString} object,
  \code{complement(x)} returns an object where each base in \code{x}
  is "complemented" i.e. A, C, G, T in a \link{DNAString} object are replaced
  by T, G, C, A respectively and A, C, G, U in a \link{RNAString} object
  are replaced by U, G, C, A respectively.

  Letters belonging to the "IUPAC extended genetic alphabet"
  are also replaced by their complement (M <-> K, R <-> Y, S <-> S, V <-> B,
  W <-> W, H <-> D, N <-> N) and the gap letter (-) is unchanged.
  \code{reverseComplement(x)} is equivalent to \code{reverse(complement(x))}
  but is faster and more memory efficient.
}

\value{
  An object of the same class and length as the original object.
}

\seealso{
  \link{DNAString-class},
  \link{RNAString-class},
  \link{DNAStringSet-class},
  \link{RNAStringSet-class},
  \link{BStringViews-class},
  \code{\link{strrev}},
  \code{\link{findPalindromes}}
}

\examples{
  reverseComplement(DNAString("ACGT-YN-"))

  ## Applying reverseComplement() to the pattern before calling matchPattern()
  ## is the standard way to search hits on the reverse strand of a chromosome:
  library(BSgenome.Dmelanogaster.FlyBase.r51)
  chrX <- Dmelanogaster[["X"]]
  pattern <- DNAString("GAACGGTGTCT")
  matchPattern(pattern, chrX) # 1 hit on strand +
  m0 <- matchPattern(reverseComplement(pattern), chrX) # 2 hits on strand -

  ## Applying reverseComplement() to the subject instead of the pattern is not
  ## a good idea for 2 reasons:
  ## (1) Chromosome sequences are generally huge so it's going to be a lot of
  ##     work and require a lot of memory to compute reverseComplement(subject).
  ## (2) Chromosome locations are generally given relatively to the positive
  ##     strand, even for features located in the negative strand, so after
  ##     doing this:
  m1 <- matchPattern(pattern, reverseComplement(chrX))
  ##     the start/end of the matches are now relative to the negative strand.
  ##     You need to apply reverseComplement() again on the result if you want
  ##     them to be relative to the positive strand:
  m2 <- reverseComplement(m1)
  ##     and finally to apply rev() to sort the matches from left to right
  ##     (5'3' direction) like in m0:
  m3 <- rev(m2) # same as m0, finally!

  ## Don't try the above example on human chromosome 1 since your computer
  ## would need to allocate about 250Mb of memory for this:
  if (FALSE) {
    library(BSgenome.Hsapiens.UCSC.hg18)
    chr1 <- Hsapiens$chr1
    matchPattern(pattern, reverseComplement(chr1)) # DON'T DO THIS!
    matchPattern(reverseComplement(pattern), chr1) # DO THIS INSTEAD
  }
}

\keyword{methods}
\keyword{manip}

