\name{translate}

\alias{translate}
\alias{translate,DNAStringSet-method}
\alias{translate,RNAStringSet-method}
\alias{translate,DNAString-method}
\alias{translate,RNAString-method}
\alias{translate,MaskedDNAString-method}
\alias{translate,MaskedRNAString-method}

\alias{codons}
\alias{codons,DNAString-method}
\alias{codons,RNAString-method}
\alias{codons,MaskedDNAString-method}
\alias{codons,MaskedRNAString-method}


\title{Translating DNA/RNA sequences}

\description{
  Functions for translating DNA or RNA sequences into amino acid sequences.
}

\usage{
## Translating DNA/RNA:
translate(x, genetic.code=GENETIC_CODE, if.fuzzy.codon="error")

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
  \item{genetic.code}{
    The genetic code to use for the translation of codons into Amino Acid
    letters. It must be represented as a named character vector of length
    64 similar to predefined constant \code{GENETIC_CODE} i.e. it must
    contain 1-letter strings in the Amino Acid alphabet and its names must
    be identical to \code{names(GENETIC_CODE)}.
    The default value for \code{genetic.code} is \code{GENETIC_CODE}
    which represents The Standard Genetic Code.
    See \code{?AA_ALPHABET}) for the Amino Acid alphabet and
    \code{?GENETIC_CODE} for The Standard Genetic Code and its known
    variants.
  }
  \item{if.fuzzy.codon}{
    How fuzzy codons (i.e codon with IUPAC ambiguities) should be handled.
    Accepted values are:
    \itemize{
      \item \code{"error"}: An error will be raised on the first occurence of
            a fuzzy codon. This is the default.
      \item \code{"solve"}: Fuzzy codons that can be translated non
            ambiguously to an amino acid or to * (stop codon) will be
            translated. Ambiguous fuzzy codons will be translated to X.
      \item \code{"error.if.X"}: Fuzzy codons that can be translated non
            ambiguously to an amino acid or to * (stop codon) will be
            translated. An error will be raised on the first occurence of
            an ambiguous fuzzy codon.
      \item \code{"X"}: All fuzzy codons (ambiguous and non-ambiguous) will
            be translated to X.
    }
    Alternatively \code{if.fuzzy.codon} can be specified as a character
    vector of length 2. The 1st string and 2nd strings specify how to handle
    non-ambiguous and ambiguous fuzzy codons, respectively.
    The accepted values for the 1st string are:
    \itemize{
      \item \code{"error"}: Any occurence of a non-ambiguous fuzzy codon
            will cause an error.
      \item \code{"solve"}: Non-ambiguous fuzzy codons will be translated to
            an amino acid or to *.
      \item \code{"X"}: Non-ambiguous fuzzy codons will be translated to X.
    }
    The accepted values for the 2nd string are:
    \itemize{
      \item \code{"error"}: Any occurence of an ambiguous fuzzy codon will
            cause an error.
      \item \code{"X"}: Ambiguous fuzzy codons will be translated to X.
    }
    All the 6 possible combinations of 1st and 2nd strings are supported.
    Note that \code{if.fuzzy.codon=c("error", "error")} is equivalent to
    \code{if.fuzzy.codon="error"},
    \code{if.fuzzy.codon=c("solve", "X")} is equivalent to
    \code{if.fuzzy.codon="solve"},
    \code{if.fuzzy.codon=c("solve", "error")} is equivalent to
    \code{if.fuzzy.codon="error.if.X"}, and
    \code{if.fuzzy.codon=c("X", "X")} is equivalent to
    \code{if.fuzzy.codon="X"}.
  }
}

\details{
  \code{translate} reproduces the biological process of RNA
  translation that occurs in the cell.
  The input of the function can be either RNA or coding DNA.
  By default The Standard Genetic Code (see \code{?\link{GENETIC_CODE}})
  is used to translate codons into amino acids but the user can
  supply a different genetic code via the \code{genetic.code} argument.

  \code{codons} is a utility for extracting the codons involved
  in this translation without translating them. 
}

\value{
  For \code{translate}: An \link{AAString} object when \code{x} is a
  \link{DNAString}, \link{RNAString}, \link{MaskedDNAString}, or
  \link{MaskedRNAString} object.
  An \link{AAStringSet} object \emph{parallel} to \code{x} (i.e. with 1
  amino acid sequence per DNA or RNA sequence in \code{x}) when \code{x}
  is a \link{DNAStringSet} or \link{RNAStringSet} object. If \code{x} has
  names on it, they're propagated to the returned object.

  For \code{codons}: An \link{XStringViews} object with 1 view per codon.
  When \code{x} is a \link{MaskedDNAString} or \link{MaskedRNAString} object,
  its masked parts are interpreted as introns and filled with the + letter
  in the returned object. Therefore codons that span across masked regions
  are represented by views that have a width > 3 and contain the + letter.
  Note that each view is guaranteed to contain exactly 3 base letters.
}

\seealso{
  \itemize{
    \item \code{\link{AA_ALPHABET}} for the Amino Acid alphabet.

    \item \code{\link{GENETIC_CODE}} for The Standard Genetic Code and
          its known variants.

    \item The examples for
          \code{\link[GenomicFeatures]{extractTranscriptsFromGenome}}
          in the \pkg{GenomicFeatures} package for computing the
          full proteome of a given organism.

    \item The \code{\link{reverseComplement}} function.

    \item The \link{DNAStringSet} and \link{AAStringSet} classes.

    \item The \link{XStringViews} and \link{MaskedXString} classes.
  }
}

\examples{
## ---------------------------------------------------------------------
## 1. BASIC EXAMPLES
## ---------------------------------------------------------------------
dna1 <- DNAString("TATAAATGGAGTAGATAA")
translate(dna1)

SGC1 <- getGeneticCode("SGC1")  # Vertebrate Mitochondrial code
translate(dna1, genetic.code=SGC1)

## All codons except 1st are fuzzy:
dna2 <- DNAString("TATANATGRAGYMGRTRA")

\dontrun{
  translate(dna2)  # error because of fuzzy codons
}
## Codons 4 to 6 are non-ambiguous and can be solved. 2nd and 3rd codons
## are ambiguous and are translated to X:
translate(dna2, if.fuzzy.codon="solve")

## Fuzzy codons that are non-ambiguous with a given genetic code can
## become ambiguous with another genetic code and vice versa:
translate(dna2, genetic.code=SGC1, if.fuzzy.codon="solve")

## ---------------------------------------------------------------------
## 2. TRANSLATING AN OPEN READING FRAME
## ---------------------------------------------------------------------
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

## Codons:
codons(y1)
which(width(codons(y1)) != 3)
codons(y1)[20:28]

## ---------------------------------------------------------------------
## 3. AN ADVANCED EXAMPLE
## ---------------------------------------------------------------------
## Translation on the '-' strand:
dna3 <- DNAStringSet(c("ATC", "GCTG", "CGACT"))
translate(reverseComplement(dna3))

## Translate sequences on both '+' and '-' strand across all 
## possible reading frames (i.e., codon position 1, 2 or 3):
## First create a DNAStringSet of '+' and '-' strand sequences, 
## removing the nucleotides prior to the reading frame start position.
dna3_subseqs <- lapply(1:3, function(pos) 
    subseq(c(dna3, reverseComplement(dna3)), start=pos))
## Translation of 'dna3_subseqs' produces a list of length 3, each with
## 6 elements (3 '+' strand results followed by 3 '-' strand results).
lapply(dna3_subseqs, translate)

## Note that translate() throws a warning when the length of the sequence
## is not divisible by 3. To avoid this warning wrap the function in 
## suppressWarnings().
}

\keyword{methods}
\keyword{manip}
