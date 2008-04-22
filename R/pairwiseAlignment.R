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
         quality1 = 1,
         quality2 = 1,
         substitutionMatrix,
         gapOpening = -5,
         gapExtension = -2,
         type = "global",
         scoreOnly = FALSE)
{
  if (class(string1) != class(string2))
    stop("'string1' and 'string2' must have the same class")
  if (!(length(quality1) %in% c(1, nchar(string1))))
    stop("length(quality1) must be 1 or nchar(string1)")
  if (!(length(quality2) %in% c(1, nchar(string2))))
    stop("length(quality2) must be 1 or nchar(string2)")
  if (is.character(quality1)) {
    quality1 <- as.numeric(charToRaw(paste(quality1, collapse = ""))) - 33
    quality1 <- quality1 / max(quality1)
  }
  if (is.character(quality2)) {
    quality2 <- as.numeric(charToRaw(paste(quality2, collapse = ""))) - 33
    quality2 <- quality2 / max(quality2)
  }
  if (any(is.na(quality1)) || any(quality1 < 0) || any(quality1 > 1))
    stop("all elements in 'quality1' must result in values between 0 and 1")
  if (any(is.na(quality2)) || any(quality2 < 0) || any(quality2 > 1))
    stop("all elements in 'quality2' must result in values between 0 and 1")
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
  substitutionMatrix <-
    matrix(as.double(substitutionMatrix),
           nrow = nrow(substitutionMatrix), ncol = ncol(substitutionMatrix),
           dimnames = dimnames(substitutionMatrix))
  gapOpening <- as.double(- abs(gapOpening))
  if (is.na(gapOpening) || length(gapOpening) != 1)
    stop("'gapOpening' must be a non-positive numeric vector of length 1")
  gapExtension <- as.double(- abs(gapExtension))
  if (is.na(gapExtension) || length(gapExtension) != 1)
    stop("'gapExtension' must be a non-positive numeric vector of length 1")
  type <- match.arg(tolower(type), c("global", "local", "overlap"))
  typeCode <- c("global" = 1L, "local" = 2L, "overlap" = 3L)[[type]]
  scoreOnly <- as.logical(scoreOnly)
  if (length(scoreOnly) != 1 || any(is.na(scoreOnly)))
    stop("'scoreOnly' must be a non-missing logical value")
  if (is.null(codec(string1))) {
    codes <-
      as.integer(charToRaw(paste(rownames(substitutionMatrix), collapse="")))
    gapCode <- charToRaw("-")
  } else {
    if (!all(rownames(substitutionMatrix) %in% alphabet(string1)))
      stop("matrix 'substitutionMatrix' is incompatible with 'string1' alphabet")
    lettersToCodes <- codec(string1)@codes
    names(lettersToCodes) <- codec(string1)@letters
    codes <- lettersToCodes[rownames(substitutionMatrix)]
    gapCode <- as.raw(lettersToCodes[["-"]])
  }
  lookupTable <- buildLookupTable(codes, 0:(nrow(substitutionMatrix) - 1))
  answer <- .Call("align_pairwiseAlignment",
                  string1,
                  string2,
                  quality1,
                  quality2,
                  gapCode,
                  typeCode,
                  scoreOnly,
                  lookupTable,
                  substitutionMatrix,
                  dim(substitutionMatrix),
                  gapOpening,
                  gapExtension,
                  PACKAGE="Biostrings")
  if (scoreOnly) {
    output <- answer[["score"]]
  } else {
    align1 <-
      new(class(string1),
          xdata = answer[["align1"]],
          length = length(answer[["align1"]]))
    align2 <-
      new(class(string2),
          xdata = answer[["align2"]],
          length = length(answer[["align2"]]))
    output <- new("XStringAlign",
                  align1 = align1,
                  align2 = align2,
                  quality1 = quality1,
                  quality2 = quality2,
                  type = type,
                  substitutionMatrix = substitutionMatrix,
                  gapOpening = gapOpening,
                  gapExtension = gapExtension,
                  score = answer[["score"]])
  }
  return(output)
}


setGeneric("pairwiseAlignment", signature = c("string1", "string2"),
           function(string1, string2, quality1 = 1, quality2 = 1,
                    substitutionMatrix, gapOpening = -5, gapExtension = -2,
                    type = "global", scoreOnly = FALSE)
           standardGeneric("pairwiseAlignment"))

setMethod("pairwiseAlignment",
          signature(string1 = "character", string2 = "character"),
          function(string1, string2, quality1 = 1, quality2 = 1,
                   substitutionMatrix, gapOpening = -5, gapExtension = -2,
                   type = "global", scoreOnly = FALSE)
          XString.pairwiseAlignment(BString(string1), BString(string2),
                                    quality1 = quality1,
                                    quality2 = quality2,
                                    substitutionMatrix = substitutionMatrix,
                                    gapExtension = gapExtension,
                                    gapOpening = gapOpening,
                                    type = type,
                                    scoreOnly = scoreOnly))

setMethod("pairwiseAlignment",
          signature(string1 = "character", string2 = "XString"),
          function(string1, string2, quality1 = 1, quality2 = 1,
                   substitutionMatrix, gapOpening = -5, gapExtension = -2,
                   type = "global", scoreOnly = FALSE)
          XString.pairwiseAlignment(XString(class(string2), string1), string2,
                                    quality1 = quality1,
                                    quality2 = quality2,
                                    substitutionMatrix = substitutionMatrix,
                                    gapExtension = gapExtension,
                                    gapOpening = gapOpening,
                                    type = type,
                                    scoreOnly = scoreOnly))

setMethod("pairwiseAlignment",
          signature(string1 = "XString", string2 = "character"),
          function(string1, string2, quality1 = 1, quality2 = 1,
                   substitutionMatrix, gapOpening = -5, gapExtension = -2,
                   type = "global", scoreOnly = FALSE)
          XString.pairwiseAlignment(string1, XString(class(string1), string2),
                                    quality1 = quality1,
                                    quality2 = quality2,
                                    substitutionMatrix = substitutionMatrix,
                                    gapExtension = gapExtension,
                                    gapOpening = gapOpening,
                                    type = type,
                                    scoreOnly = scoreOnly))

setMethod("pairwiseAlignment",
          signature(string1 = "XString", string2 = "XString"),
          function(string1, string2, quality1 = 1, quality2 = 1,
                   substitutionMatrix, gapOpening = -5, gapExtension = -2,
                   type = "global", scoreOnly = FALSE)
          XString.pairwiseAlignment(string1, string2,
                                    quality1 = quality1,
                                    quality2 = quality2,
                                    substitutionMatrix = substitutionMatrix,
                                    gapExtension = gapExtension,
                                    gapOpening = gapOpening,
                                    type = type,
                                    scoreOnly = scoreOnly))
