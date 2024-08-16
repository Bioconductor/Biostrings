\name{chartr}

\alias{chartr}
\alias{chartr,ANY,ANY,XString-method}
\alias{chartr,ANY,ANY,XStringSet-method}
\alias{chartr,ANY,ANY,XStringViews-method}
\alias{chartr,ANY,ANY,MaskedXString-method}

\alias{replaceAmbiguities}

\title{Replace letters in a sequence or set of sequences}

\description{
  Replace letters in a sequence or set of sequences.
}

\usage{
\S4method{chartr}{ANY,ANY,XString}(old, new, x)

replaceAmbiguities(x, new="N")
}

\arguments{
  \item{old}{
    A character string specifying the characters to be replaced.
  }
  \item{new}{
    A character string specifying the replacements. It must be a single
    letter for \code{replaceAmbiguities}.
  }
  \item{x}{
    The sequence or set of sequences to translate.
    If \code{x} is an \link{XString}, \link{XStringSet}, \link{XStringViews}
    or \link{MaskedXString} object, then the appropriate \code{chartr} method
    is called, otherwise the standard \code{\link[base]{chartr}} R function
    is called.
  }
}

\details{
  See \code{?\link[base]{chartr}} for the details.

  Note that, unlike the standard \code{\link[base]{chartr}} R function,
  the methods for \link{XString}, \link{XStringSet}, \link{XStringViews}
  and \link{MaskedXString} objects do NOT support character ranges in the
  specifications.

  \code{replaceAmbiguities()} is a simple wrapper around \code{chartr()}
  that replaces all IUPAC ambiguities with N for objects containing DNA
  or RNA sequence data.
}

\value{
  An object of the same class and length as the original object.
}

\seealso{
  \itemize{
    \item \code{\link[base]{chartr}} in the \pkg{base} package.

    \item The \code{\link{replaceAt}} function for extracting or replacing
          arbitrary subsequences from/in a sequence or set of sequences.

    \item The \code{\link{replaceLetterAt}} function for a DNA-specific
          single-letter replacement functions useful for SNP injections.

    \item \code{\link{IUPAC_CODE_MAP}} for the mapping between IUPAC
          nucleotide ambiguity codes and their meaning.

    \item \code{\link{alphabetFrequency}} (and \code{\link{uniqueLetters}})
          for tabulating letters in (and extracting the unique letters from)
          a sequence or set of sequences.

    \item The \link{XString}, \link{XStringSet}, \link{XStringViews},
          and \link{MaskedXString} classes.
  }
}

\examples{
## ---------------------------------------------------------------------
## A BASIC chartr() EXAMPLE
## ---------------------------------------------------------------------

x <- BString("MiXeD cAsE 123")
chartr("iXs", "why", x)

## ---------------------------------------------------------------------
## TRANSFORMING DNA WITH BISULFITE (AND SEARCHING IT...)
## ---------------------------------------------------------------------

library(BSgenome.Celegans.UCSC.ce2)
chrII <- Celegans[["chrII"]]
alphabetFrequency(chrII)
pattern <- DNAString("TGGGTGTATTTA")

## Transforming and searching the + strand
plus_strand <- chartr("C", "T", chrII)
alphabetFrequency(plus_strand)
matchPattern(pattern, plus_strand)
matchPattern(pattern, chrII)

## Transforming and searching the - strand
minus_strand <- chartr("G", "A", chrII)
alphabetFrequency(minus_strand)
matchPattern(reverseComplement(pattern), minus_strand)
matchPattern(reverseComplement(pattern), chrII)

## ---------------------------------------------------------------------
## replaceAmbiguities()
## ---------------------------------------------------------------------

dna <- DNAStringSet(c("TTTKYTT-GR", "", "NAASACVT"))
dna
replaceAmbiguities(dna)
}

\keyword{methods}
\keyword{manip}
