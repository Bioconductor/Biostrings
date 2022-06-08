\name{seqinfo-methods}

\alias{seqinfo}
\alias{seqinfo,DNAStringSet-method}
\alias{seqinfo<-,DNAStringSet-method}

\title{seqinfo() method for DNAStringSet objects}

\description{
    \code{\link[GenomeInfoDb]{seqinfo}} methods for extracting the
    sequence information stored in a \link[Biostrings]{DNAStringSet}
    object.
}

\usage{
\S4method{seqinfo}{DNAStringSet}(x)
}

\arguments{
  \item{x}{
    A \link[Biostrings]{DNAStringSet} object.
  }
}

\value{
    A Seqinfo object for the 'seqinfo' getter.

    A \link[Biostrings]{DNAStringSet} object containing sequence
    information for the 'seqinfo' setter.
}

\seealso{
  \code{\link[Biostrings]{getSeq}},
  \link[Biostrings]{DNAStringSet-class},
}

\examples{
## ---------------------------------------------------------------------
## A. SIMPLE EXAMPLE
## ---------------------------------------------------------------------

library(drosophila2probe)
probes <- DNAStringSet(drosophila2probe)

## Check metadata slot: empty
metadata(probes)

## Get generated seqinfo table
seqinfo(probes)

## Subsetting seqinfo table to 10 seqnames
probes10 <- probes[1:10]
seqinfo(probes10) <- seqinfo(probes)[as.character(1:10)]
## See result: 10 seqnames
seqinfo(probes10)

}
