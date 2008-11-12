\name{translate}

\alias{transcribe}
\alias{cDNA}
\alias{dna2rna}
\alias{rna2dna}

\alias{codons}
\alias{codons,DNAString-method}
\alias{codons,RNAString-method}
\alias{codons,MaskedDNAString-method}
\alias{codons,MaskedRNAString-method}

\alias{translate}
\alias{translate,DNAString-method}
\alias{translate,RNAString-method}
\alias{translate,DNAStringSet-method}
\alias{translate,RNAStringSet-method}
\alias{translate,MaskedDNAString-method}
\alias{translate,MaskedRNAString-method}


\title{DNA/RNA transcription and translation}

\description{
  Functions for transcription and/or translation of DNA or RNA
  sequences, and related utilities.
}

\usage{
  transcribe(x)
  cDNA(x)
  codons(x)
  translate(x)

  ## Related utilities
  dna2rna(x)
  rna2dna(x)
}

\arguments{
  \item{x}{
    A \link{DNAString} object for \code{transcribe} and \code{dna2rna}.

    An \link{RNAString} object for \code{cDNA} and \code{rna2dna}.

    A \link{DNAString}, \link{RNAString}, \link{MaskedDNAString} or
    \link{MaskedRNAString} object for \code{codons}.

    A \link{DNAString}, \link{RNAString}, \link{DNAStringSet},
    \link{RNAStringSet}, \link{MaskedDNAString} or \link{MaskedRNAString}
    object for \code{translate}.
  }
}

\value{
  An \link{RNAString} object for \code{transcribe} and \code{dna2rna}.

  A \link{DNAString} object for \code{cDNA} and \code{rna2dna}.

  An \link{XStringViews} object with 1 view per codon for \code{codons}.
  When \code{x} is a \link{MaskedDNAString} or \link{MaskedRNAString}
  object, its masked parts are interpreted as introns and filled with
  the + letter in the returned object. Therefore codons that span across
  masked regions are represented by views that have a width > 3 and
  contain the + letter. Note that each view is guaranteed to contain
  exactly 3 base letters.

  An \link{AAString} object for \code{translate}.
}

\seealso{
  \code{\link{reverseComplement}},
  \code{\link{GENETIC_CODE}},
  \link{DNAString-class},
  \link{RNAString-class},
  \link{AAString-class},
  \link{DNAStringSet-class},
  \link{RNAStringSet-class},
  \link{XStringViews-class}
}

\examples{
  file <- system.file("extdata", "someORF.fa", package="Biostrings")
  x <- read.DNAStringSet(file, "fasta")
  x

  ## The first and last 1000 nucleotides are not part of the ORFs:
  x <- DNAStringSet(x, start=1001, end=-1001)

  ## Before calling translate() on an ORF, we need to mask the introns
  ## if any. We can get this information fron the SGD database
  ## (http://www.yeastgenome.org/).
  ## According to SGD, the 1st ORF (YAL001C) has an intron at 71..160
  ## (see http://db.yeastgenome.org/cgi-bin/locus.pl?locus=YAL001C)
  y1 <- x[[1]]
  mask1 <- Mask(length(y1), start=71, end=160)
  masks(y1) <- mask1
  y1
  translate(y1)

  ## Codons
  codons(y1)
  which(width(codons(y1)) != 3)
  codons(y1)[20:28]
}

\keyword{methods}
\keyword{manip}
