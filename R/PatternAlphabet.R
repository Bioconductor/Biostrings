#Copyright (C) 2003 by Saikat DebRoy
setClass("BioPatternAlphabet",
         representation(baseAlphabet="BioAlphabet"),
         contains="BioAlphabet")

setReplaceMethod("gapletter",
                 signature(x = "BioPatternAlphabet",
                           value = "character"),
                 function (x, value)
             {
                 gapletter(x@baseAlphabet) <- value
                 callNextMethod()
             })

setMethod("initialize",
          signature(.Object = "BioPatternAlphabet"),
          function (.Object, baseAlphabet, letters, ...)
      {
          .Object@baseAlphabet <- baseAlphabet
          if (length(letters) == 0) {
              .Object@mapping <- baseAlphabet@mapping
              .Object@letters <- baseAlphabet@letters
              return(.Object)
          }
          alphletters <- names(baseAlphabet@mapping)
          mapping <- letters
          letters <- names(letters)
          if (!is.character(letters) || any(nchar(letters) != 1))
              stop("pattern alphabet can only have single letters")
          if (!all(match(letters, alphletters, nomatch = 0) == 0))
              stop("must add new letters in pattern alphabet")
          mapping <- unlist(lapply(mapping,
                                   function(map)
                               {
                                   map <- strsplit(map, '')[[1]]
                                   if (!all(map %in%
                                            alphletters))
                                       stop("each pattern letter can only match letters in the original alphabet")
                                   .Call("IntegerBitOr", baseAlphabet@mapping[map])
                               }))
          names(mapping) <- letters
          .Object@mapping <- sort(c(baseAlphabet@mapping, mapping))
          .Object@letters <- paste(baseAlphabet@letters,
                                   paste(letters, collapse=''),
                                   sep='')
          .Object@gap <- baseAlphabet@gap
          .Object
      })

DNAPatternAlphabet <- function()
    new("BioPatternAlphabet",
        DNAAlphabet(), c(N="AGCT", # N matches everything but the gap character
                         B="CGT",
                         D="AGT",
                         H="ACT",
                         K="GT",
                         M="AC",
                         R="AG",
                         S="CG",
                         V="ACG",
                         W="AT",
                         Y="CT"))

RNAPatternAlphabet <- function()
    new("BioPatternAlphabet",
        RNAAlphabet(), c(N="AGCU", # N matches everything but the gap character
                         B="CGU",
                         D="AGU",
                         H="ACU",
                         K="GU",
                         M="AC",
                         R="AG",
                         S="CG",
                         V="ACG",
                         W="AU",
                         Y="CU"))
