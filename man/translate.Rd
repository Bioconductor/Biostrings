\name{translate}

\alias{translate}
\alias{translate,DNAString-method}
\alias{translate,RNAString-method}
\alias{translate,DNAStringSet-method}
\alias{translate,RNAStringSet-method}
\alias{translate,MaskedDNAString-method}
\alias{translate,MaskedRNAString-method}

\alias{codons}
\alias{codons,DNAString-method}
\alias{codons,RNAString-method}
\alias{codons,MaskedDNAString-method}
\alias{codons,MaskedRNAString-method}

% Old stuff (Deprecated or Defunct):
\alias{transcribe}
\alias{cDNA}
\alias{dna2rna}
\alias{rna2dna}


\title{Translating DNA/RNA sequences}

\description{
  Functions for translating DNA or RNA sequences into amino acid sequences.
}

\usage{
## Translating DNA/RNA:
translate(x)

## Extracting codons without translating them:
codons(x)
}

\arguments{
  \item{x}{
    A \link{DNAStringSet}, \link{RNAStringSet}, \link{DNAString},
    \link{RNAString}, \link{MaskedDNAString} or \link{MaskedRNAString}
    object for \code{translate}.

    A \link{DNAString}, \link{RNAString}, \link{MaskedDNAString} or
    \link{MaskedRNAString} object for \code{codons}.
  }
}

\details{
  \code{translate} reproduces the biological process of RNA
  translation that occurs in the cell.
  The input of the function can be either RNA or coding DNA.
  The Standard Genetic Code (see \code{?\link{GENETIC_CODE}}) is
  used to translate codons into amino acids.

  \code{codons} is a utility for extracting the codons involved
  in this translation without translating them. 
}

\value{
  An \link{AAStringSet} or \link{AAString} object for \code{translate}.

  An \link{XStringViews} object with 1 view per codon for \code{codons}.
  When \code{x} is a \link{MaskedDNAString} or \link{MaskedRNAString}
  object, its masked parts are interpreted as introns and filled with
  the + letter in the returned object. Therefore codons that span across
  masked regions are represented by views that have a width > 3 and
  contain the + letter. Note that each view is guaranteed to contain
  exactly 3 base letters.
}

\seealso{
  \itemize{
    \item The \link{DNAStringSet} and \link{AAStringSet} classes.

    \item The \link{XStringViews} and \link{MaskedXString} classes.

    \item \code{\link{GENETIC_CODE}} for The Standard Genetic Code.

    \item The \code{\link{reverseComplement}} function.

    \item The examples for
          \code{\link[GenomicFeatures]{extractTranscriptsFromGenome}}
          in the \pkg{GenomicFeatures} package for computing the
          full proteome of a given organism.
  }
}

\examples{
file <- system.file("extdata", "someORF.fa", package="Biostrings")
x <- readDNAStringSet(file)
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

## Translation on the '-' strand:
dna <- DNAStringSet(c("ATC", "GCTG", "CGACT"))
translate(reverseComplement(dna))

## Translate sequences on both '+' and '-' strand across all 
## possible reading frames (i.e., codon position 1, 2 or 3):
## First create a DNAStringSet of '+' and '-' strand sequences, 
## removing the nucleotides prior to the reading frame start position.
dnaSubseq <- lapply(1:3, function(pos, dna) 
    subseq(c(dna, reverseComplement(dna)), start=pos), dna)
## Translation of 'dnaSubseq' produces a list of length 3, each with 6 
## elements (3 '+' strand results followed by 3 '-' strand results).
lapply(dnaSubseq, translate)

## translate() throws a warning when the length of the sequence is not 
## divisible by 3. To avoid this warning wrap the function in 
## suppressWarnings(). 
}

\keyword{methods}
\keyword{manip}
