library("methods")

setClass("BioAlphabet",
         representation(letters="character",
                        mapping="integer"),
         contains="VIRTUAL")

setMethod("initialize",
          signature(.Object = "BioAlphabet"),
          function (.Object, letters, ...)
      {
          if (!is.character(letters) ||
              any(nchar(letters) != 1))
              stop("invalid alphabet letters")
          if (!any(letters == '-'))
              letters <- c('-', letters)
          .Object@letters <- paste(letters, collapse='')
          .Object@mapping <- sort(as.integer(2^(0:(length(letters)-1))))
          names(.Object@mapping) <- letters
          .Object
      })

setClass("NucleotideAlphabet",
         contains="BioAlphabet")

setMethod("initialize",
          signature(.Object = "NucleotideAlphabet"),
          function (.Object, letters, ...)
      {
          tmp <- sort(toupper(letters))
          if (length(tmp) != 4 &&
              (tmp[1] != '-' || length(tmp) != 5))
              stop("Nucleotide alphabets have 4 letters other than gap")
          if (length(tmp) == 5)
              tmp <- tmp[-1]
          if (tmp != c('A', 'C', 'G', 'T') &&
              tmp != c('A', 'C', 'G', 'U'))
              warning("given letters, (",
                      paste(letters, collapse=", "),
                      "), do not match usual RNA or DNA base encodings")
          callNextMethod(.Object, letters=letters, ...)
      })

setClass("AminoAcidAlphabet", # should also have the codons
         contains="BioAlphabet")

DNAAlphabet <- function()
    new("NucleotideAlphabet",
        letters=c('A', 'G', 'C', 'T'))
RNAAlphabet <- function()
    new("NucleotideAlphabet",
        letters=c('A', 'G', 'C', 'U'))
