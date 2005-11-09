#Copyright (C) 2003-2004 by Saikat DebRoy
library("methods")

setClass("BioAlphabet",
         representation(letters="character",
                        mapping="integer",
                        gap="character"),
         contains="VIRTUAL")

setGeneric("gapletter<-",
               function(x, value)
               standardGeneric("gapletter<-"))

setReplaceMethod("gapletter",
                 signature(x = "BioAlphabet", value = "character"),
                 function (x, value)
             {
                 if (length(value) != 1 || nchar(value) != 1)
                     stop("gap must be a single letter")
                 gapindex <- match(x@gap, names(x@mapping))
                 if (is.na(gapindex))
                     stop("gap character not in mapping")
                 names(x@mapping)[gapindex] <- value
                 substr(x@letters, gapindex, gapindex) <- value
                 x@gap <- value
                 x
             })

setMethod("initialize",
          signature(.Object = "BioAlphabet"),
          function (.Object, letters, gap, ...)
      {
          if (!is.character(letters) ||
              any(nchar(letters) != 1))
              stop("invalid alphabet letters")
          if (missing(gap))
              gap <- '-'
          else if (!is.character(gap) || length(gap) != 1 ||
                   nchar(gap) != 1)
              stop("gap must be a single character")
          if (any(letters == gap))
              letters <- c(letters[letters == gap], letters[letters != gap])
          else letters <- c(gap, letters)
          .Object@letters <- paste(letters, collapse='')
          .Object@mapping <- sort(as.integer(2^(0:(length(letters)-1))))
          names(.Object@mapping) <- letters
          .Object@gap <- gap
          .Object
      })

setClass("NucleotideAlphabet",
         contains="BioAlphabet")

setMethod("initialize",
          signature(.Object = "NucleotideAlphabet"),
          function (.Object, letters, gap, ...)
      {
          tmp <- sort(toupper(letters))
          if (missing(gap))
              gap <- '-'
          else if (is.character(gap) || length(gap) != 1 ||
                   nchar(gap) != 1)
              stop("gap must be a single character")
          if ((length(tmp) != 4 &&
               (tmp[1] != toupper(gap) || length(tmp) != 5)) ||
              tmp != unique(tmp))
              stop("Nucleotide alphabets have 4 distinct letters other than gap")
          if (length(tmp) == 5)
              tmp <- tmp[tmp != gap]
          if (tmp != c('A', 'C', 'G', 'T') &&
              tmp != c('A', 'C', 'G', 'U'))
              warning("given letters, (",
                      paste(letters, collapse=", "),
                      "), do not match usual RNA or DNA base encodings")
          callNextMethod(.Object, letters=letters, gap=gap, ...)
      })

setClass("AminoAcidAlphabet", # should also have the codons
         contains="BioAlphabet")

DNAAlphabet <- function()
    new("NucleotideAlphabet",
        letters=c('T', 'G', 'C', 'A'))
RNAAlphabet <- function()
    new("NucleotideAlphabet",
        letters=c('U', 'G', 'C', 'A'))

# StandardPeptideAlphabet <- function()
#     new("AminoAcidAlphabet",
#         letters= c('A', # alanine
#                    'B', # aspartatic acid or asparagine
#                    'C', # cysteine
#                    'D', # aspartic acid
#                    'E', # glutamic acid
#                    'F', # phenylalanine
#                    'G', # glycine
#                    'H', # histidine
#                    'I', # isoleucine
#                    'K', # lysine
#                    'L', # leucine
#                    'M', # methionine
#                    'N', # asparagine
#                    'P', # proline
#                    'Q', # glutamine
#                    'R', # arginine
#                    'S', # serine
#                    'T', # threonine
#                    'U', # selenocysteine
#                    'V', # valine
#                    'W', # tryptophan
#                    'Y', # tyrosine
#                    'Z', # glutamaic acid or glutamine
#                    '*'  # translation stop
#                    ))
