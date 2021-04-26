\name{getSeq-methods}

\alias{getSeq-methods}
\alias{getSeq,XStringSet-method}

\alias{[,XStringSet,GenomicRanges-method}
\alias{[,XStringSet,GRangesList-method}

\title{getSeq method for XStringSet objects}

\description{
   \code{\link[Biostrings]{getSeq}} methods for extracting a set of
   sequences (or subsequences) from a \link[Biostrings]{XStringSet}
   object. For XStringSets, there are also
   convenience methods on \code{[} that delegate to \code{getSeq}.
}

\usage{
\S4method{getSeq}{XStringSet}(x, names)
}

\arguments{
  \item{x}{
    A \link[Biostrings]{XStringSet} object.
  }
  \item{names}{
    \code{names} must be a character vector,
    \link[GenomicRanges]{GRanges} or
    \link[GenomicRanges]{GRangesList} object.
  }
}

\details{
  L, the number of sequences to extract, is determined as follow:
  \itemize{
  \item If \code{names} is a \link[GenomicRanges]{GRanges} or
        \link[IRanges]{IntegerRanges} object then L = \code{length(names)}.
  \item If \code{names} is a \link[GenomicRanges]{GRangesList} or
        \link[IRanges]{IntegerRangesList} object then
        L = \code{length(unlist(names))}.
  \item Otherwise, L is the length of the longest of \code{names},
        \code{start}, \code{end} and \code{width} and all these
        arguments are recycled to this length.
        \code{NA}s and negative values in these 3 arguments are
        solved according to the rules of the SEW (Start/End/Width)
        interface (see \code{?\link[IRanges]{solveUserSEW}} for
        the details).
  }

  If \code{names} is neither a \link[GenomicRanges]{GRanges} or
  \link[GenomicRanges]{GRangesList} object, then the \code{strand}
  argument is also recycled to length L.

  There are convenience methods for extracting sequences from
  \link[Biostrings]{XStringSet} objects using a
  \link[GenomicRanges]{GenomicRanges} or \link[GenomicRanges]{GRangesList}
  subscript (character subscripts are implicitly supported). Both methods
  are simple wrappers around \code{getSeq}, although the GRangesList method
  differs from the \code{getSeq} behavior in that the within-element results
  are concatenated and returned as an XStringSet, rather than an
  XStringSetList. See the examples.
}

\value{
  Normally a \link[Biostrings]{DNAStringSet} object (or character vector
  if \code{as.character=TRUE}).

  With the 2 following exceptions:
  \enumerate{
    \item A \link[Biostrings]{DNAStringSetList} object (or
          \link[IRanges]{CharacterList} object if \code{as.character=TRUE})
          of the same shape as \code{names} if \code{names} is a
          \link[GenomicRanges]{GRangesList} object.
    \item A \link[Biostrings]{DNAString} object (or single character string
          if \code{as.character=TRUE}) if L = 1 and \code{names}
          is not a \link[GenomicRanges]{GRanges},
          \link[GenomicRanges]{GRangesList}, \link[IRanges]{IntegerRangesList},
          or \link[IRanges]{IntegerRanges} object.
  }
}

\note{
  Be aware that using \code{as.character=TRUE} can be very inefficient
  when extracting a "big" amount of DNA sequences (e.g. millions of
  short sequences or a small number of very long sequences).

  Note that the masks in \code{x}, if any, are always ignored. In other
  words, masked regions in the genome are extracted in the same way as
  unmasked regions (this is achieved by dropping the masks before extraction).
  See \code{?`\link[Biostrings]{MaskedDNAString-class}`} for more
  information about masked DNA sequences.
}

\author{H. Pagès; improvements suggested by Matt Settles and others}

\seealso{
  \code{\link[Biostrings]{getSeq}},
  \link[Biostrings]{DNAString-class},
  \link[Biostrings]{DNAStringSet-class},
  \link[Biostrings]{MaskedDNAString-class},
  \link[GenomicRanges]{GRanges-class},
  \link[GenomicRanges]{GRangesList-class},
  \link[IRanges]{IntegerRangesList-class},
  \link[IRanges]{IntegerRanges-class},
  \code{\link[base]{grep}}
}

\examples{
## Note that you need to load the package(s) where the generic resides
## to have getMethod() show the function
library(BSgenome)
getMethod("getSeq", "XStringSet")

## ---------------------------------------------------------------------
## A. SIMPLE EXAMPLE
## ---------------------------------------------------------------------

library(drosophila2probe)
probes <- DNAStringSet(drosophila2probe)

prob <- probes[1:10]
names(prob) <- letters[1:10]

getSeq(prob, c("a", "b"))

## ---------------------------------------------------------------------
## B. USING A GRanges OBJECT
## ---------------------------------------------------------------------

gr1 <- GRanges(seqnames=c("a", "b", "c"),
               ranges=IRanges(start=1:3, width=12),
               strand="+")

getSeq(prob, gr1)

## ---------------------------------------------------------------------
## D. USING A GRangesList OBJECT
## ---------------------------------------------------------------------

gr2 <- shift(gr1, 5)
gr3 <- gr2
strand(gr3) <- "-"

grl <- GRangesList(gr1, gr2, gr3)
getSeq(prob, grl)
}

\keyword{manip}

