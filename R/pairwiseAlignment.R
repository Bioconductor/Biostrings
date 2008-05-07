### =========================================================================
### The pairwiseAlignment() generic & related functions
### -------------------------------------------------------------------------
###
### The pairwiseAligment() function provides optimal pairwise alignment of
### the following types:
### - Global alignment
### - Local alignment
### - Overlap alignment
###
### -------------------------------------------------------------------------


XString.pairwiseAlignment <-
function(string1,
         string2,
         quality1 = 22L,
         quality2 = 22L,
         type = "global",
         substitutionMatrix = NULL,
         gapOpening = -10,
         gapExtension = -0.5,
         scoreOnly = FALSE)
{
  ## Check arguments
  if (class(string1) != class(string2))
    stop("'string1' and 'string2' must have the same class")
  type <- match.arg(tolower(type), c("global", "local", "overlap", "overlap1", "overlap2"))
  typeCode <-
    c("global" = 1L, "local" = 2L, "overlap" = 3L, "overlap1" = 4L, "overlap2" = 5L)[[type]]
  gapOpening <- as.double(- abs(gapOpening))
  if (is.na(gapOpening) || length(gapOpening) != 1)
    stop("'gapOpening' must be a non-positive numeric vector of length 1")
  gapExtension <- as.double(- abs(gapExtension))
  if (is.na(gapExtension) || length(gapExtension) != 1)
    stop("'gapExtension' must be a non-positive numeric vector of length 1")
  scoreOnly <- as.logical(scoreOnly)
  if (length(scoreOnly) != 1 || any(is.na(scoreOnly)))
    stop("'scoreOnly' must be a non-missing logical value")

  ## Process string information
  if (is.null(codec(string1))) {
    uniqueBases <-
      unique(c(unique(charToRaw(as.character(string1))),
               unique(charToRaw(as.character(string2)))))
    alphabetToCodes <- as.integer(uniqueBases)
    names(alphabetToCodes) <- rawToChar(uniqueBases, multiple = TRUE)
  } else {
    stringCodec <- codec(string1)
    alphabetToCodes <- stringCodec@codes
    names(alphabetToCodes) <- stringCodec@letters
  }

  ## Generate quality-based and constant substitution matrix information
  if (is.null(substitutionMatrix)) {
    if (is.numeric(quality1)) {
      if (any(is.na(quality1)) || any(quality1 < 0 || quality1 > 99))
        stop("integer 'quality1' values must be between 0 and 99")
      quality1 <- rawToChar(as.raw(33L + as.integer(quality1)))
    }
    if (!is(quality1, "XString"))
      quality1 <- BString(quality1)

    if (is.numeric(quality2)) {
      if (any(is.na(quality2)) || any(quality2 < 0 || quality2 > 99))
        stop("integer 'quality2' values must be between 0 and 99")
      quality2 <- rawToChar(as.raw(33L + as.integer(quality2)))
    }
    if (!is(quality2, "XString"))
      quality2 <- BString(quality2)

    nAlphabet <-
      switch(class(string1),
             DNAString =, RNAString = 4L,
             AAString = 20L,
             length(alphabetToCodes))

    errorProbs <- 10^seq(0, -9.9, by = -0.1)
    errorMatrix <-
      outer(errorProbs, errorProbs,
            function(e1,e2,n) e1 + e2 - (n/(n - 1)) * e1 * e2,
            n = nAlphabet)
    qualityLookupTable <- buildLookupTable(33:(99 + 33), 0:99)
    qualityMatchMatrix <- log2((1 - errorMatrix) * nAlphabet)
    qualityMismatchMatrix <- log2(errorMatrix * (nAlphabet / (nAlphabet - 1)))

    constantLookupTable <- integer(0)
    constantMatrix <- matrix(numeric(0), nrow = 0, ncol = 0)
  } else {
    quality1 <- BString("")
    quality2 <- BString("")
    qualityLookupTable <- integer(0)
    qualityMatchMatrix <- matrix(numeric(0), nrow = 0, ncol = 0)
    qualityMismatchMatrix <- matrix(numeric(0), nrow = 0, ncol = 0)
    if (is.character(substitutionMatrix)) {
      if (length(substitutionMatrix) != 1)
        stop("'substitutionMatrix' is a character vector of length != 1")
      tempMatrix <- substitutionMatrix
      substitutionMatrix <- try(getdata(tempMatrix), silent = TRUE)
      if (is(substitutionMatrix, "try-error"))
        stop("unknown scoring matrix \"", tempMatrix, "\"")
    }
    if (!is.matrix(substitutionMatrix) || !is.numeric(substitutionMatrix))
      stop("'substitutionMatrix' must be a numeric matrix")
    if (!identical(rownames(substitutionMatrix), colnames(substitutionMatrix)))
      stop("row and column names differ for matrix 'substitutionMatrix'")
    if (is.null(rownames(substitutionMatrix)))
      stop("matrix 'substitutionMatrix' must have row and column names")
    if (any(duplicated(rownames(substitutionMatrix))))
      stop("matrix 'substitutionMatrix' has duplicated row names")
    availableLetters <-
      intersect(names(alphabetToCodes), rownames(substitutionMatrix))
    constantMatrix <-
      matrix(as.double(substitutionMatrix[availableLetters, availableLetters]),
             nrow = length(availableLetters),
             ncol = length(availableLetters),
             dimnames = list(availableLetters, availableLetters))
    constantLookupTable <-
      buildLookupTable(alphabetToCodes[availableLetters],
                       0:(length(availableLetters) - 1))
  }
  answer <- .Call("align_pairwiseAlignment",
                  string1,
                  string2,
                  quality1,
                  quality2,
                  typeCode,
                  scoreOnly,
                  gapOpening,
                  gapExtension,
                  qualityLookupTable,
                  qualityMatchMatrix,
                  qualityMismatchMatrix,
                  dim(qualityMatchMatrix),
                  constantLookupTable,
                  constantMatrix,
                  dim(constantMatrix),
                  PACKAGE="Biostrings")
  if (scoreOnly) {
    output <- answer[["score"]]
  } else {
    output <- new("XStringAlign",
                  string1 = string1,
                  string2 = string2,
                  quality1 = quality1,
                  quality2 = quality2,
                  match1 =
                  IRanges(start = answer[["startMatch1"]],
                          width = answer[["widthMatch1"]]),
                  match2 =
                  IRanges(start = answer[["startMatch2"]],
                          width = answer[["widthMatch2"]]),
                  inserts1 =
                  IRanges(start = answer[["startInserts1"]],
                          width = answer[["widthInserts1"]]),
                  inserts2 =
                  IRanges(start = answer[["startInserts2"]],
                          width = answer[["widthInserts2"]]),
                  profile1 = answer[["profile1"]],
                  profile2 = answer[["profile2"]],
                  score = answer[["score"]],
                  type = type,
                  constantMatrix = constantMatrix,
                  gapOpening = gapOpening,
                  gapExtension = gapExtension)
  }
  return(output)
}


setGeneric("pairwiseAlignment", signature = c("string1", "string2"),
           function(string1, string2, quality1 = 22L, quality2 = 22L,
                    type = "global", substitutionMatrix = NULL,
                    gapOpening = -10, gapExtension = -0.5,
                    scoreOnly = FALSE)
           standardGeneric("pairwiseAlignment"))

setMethod("pairwiseAlignment",
          signature(string1 = "character", string2 = "character"),
          function(string1, string2, quality1 = 22L, quality2 = 22L,
                   type = "global", substitutionMatrix = NULL,
                   gapOpening = -10, gapExtension = -0.5,
                   scoreOnly = FALSE)
          XString.pairwiseAlignment(BString(string1), BString(string2),
                                    quality1 = quality1,
                                    quality2 = quality2,
                                    type = type,
                                    substitutionMatrix = substitutionMatrix,
                                    gapExtension = gapExtension,
                                    gapOpening = gapOpening,
                                    scoreOnly = scoreOnly))

setMethod("pairwiseAlignment",
          signature(string1 = "character", string2 = "XString"),
          function(string1, string2, quality1 = 22L, quality2 = 22L,
                   type = "global", substitutionMatrix = NULL,
                   gapOpening = -10, gapExtension = -0.5,
                   scoreOnly = FALSE)
          XString.pairwiseAlignment(XString(class(string2), string1), string2,
                                    quality1 = quality1,
                                    quality2 = quality2,
                                    type = type,
                                    substitutionMatrix = substitutionMatrix,
                                    gapExtension = gapExtension,
                                    gapOpening = gapOpening,
                                    scoreOnly = scoreOnly))

setMethod("pairwiseAlignment",
          signature(string1 = "XString", string2 = "character"),
          function(string1, string2, quality1 = 22L, quality2 = 22L,
                   type = "global", substitutionMatrix = NULL,
                   gapOpening = -10, gapExtension = -0.5,
                   scoreOnly = FALSE)
          XString.pairwiseAlignment(string1, XString(class(string1), string2),
                                    quality1 = quality1,
                                    quality2 = quality2,
                                    type = type,
                                    substitutionMatrix = substitutionMatrix,
                                    gapExtension = gapExtension,
                                    gapOpening = gapOpening,
                                    scoreOnly = scoreOnly))

setMethod("pairwiseAlignment",
          signature(string1 = "XString", string2 = "XString"),
          function(string1, string2, quality1 = 22L, quality2 = 22L,
                   type = "global", substitutionMatrix = NULL,
                   gapOpening = -10, gapExtension = -0.5,
                   scoreOnly = FALSE)
          XString.pairwiseAlignment(string1, string2,
                                    quality1 = quality1,
                                    quality2 = quality2,
                                    type = type,
                                    substitutionMatrix = substitutionMatrix,
                                    gapExtension = gapExtension,
                                    gapOpening = gapOpening,
                                    scoreOnly = scoreOnly))
