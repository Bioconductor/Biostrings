### =========================================================================
### The stringDist() generic
### -------------------------------------------------------------------------

XStringSet.stringDist <-
function(x,
         method = "levenshtein",
         ignoreCase = FALSE,
         diag = FALSE,
         upper = FALSE,
         type = "global",
         quality = PhredQuality(22L),
         substitutionMatrix = NULL,
         gapOpening = 0,
         gapExtension = -1)
{
  ## Check arguments
  method <- match.arg(method, c("levenshtein", "quality", "substitutionMatrix"))
  type <- match.arg(type, c("global", "local", "overlap"))
  typeCode <- c("global" = 1L, "local" = 2L, "overlap" = 3L)[[type]]
  gapOpening <- as.double(- abs(gapOpening))
  if (length(gapOpening) != 1 || is.na(gapOpening))
    stop("'gapOpening' must be a non-positive numeric vector of length 1")
  gapExtension <- as.double(- abs(gapExtension))
  if (length(gapExtension) != 1 || is.na(gapExtension))
    stop("'gapExtension' must be a non-positive numeric vector of length 1")

  ## Process string information
  if (is.null(codec(x))) {
    uniqueBases <- unique(unique(charToRaw(as.character(super(x)))))
    alphabetToCodes <- as.integer(uniqueBases)
    names(alphabetToCodes) <- rawToChar(uniqueBases, multiple = TRUE)
  } else {
    stringCodec <- codec(x)
    alphabetToCodes <- stringCodec@codes
    names(alphabetToCodes) <- stringCodec@letters
  }

  ## Set parameters when method == "levenshtein"
  if (method == "levenshtein") {
    type <- "global"
    typeCode <- 1L
    gapOpening <- 0
    gapExtension <- -1
    if (ignoreCase)
      caseAdjustedAlphabet <- tolower(names(alphabetToCodes))
    else
      caseAdjustedAlphabet <- names(alphabetToCodes)
    substitutionMatrix <-
      outer(caseAdjustedAlphabet, caseAdjustedAlphabet, function(x,y) -as.numeric(x!=y))
    dimnames(substitutionMatrix) <- list(names(alphabetToCodes), names(alphabetToCodes))
  }

  ## Generate quality-based and constant substitution matrix information
  if (method == "quality") {
    useQuality <- TRUE
    if (!is(quality, "XStringQuality"))
        stop("'quality' must be of class 'XStringSet'")
    if (!all(nchar(quality) == 1 | nchar(quality) == nchar(x)))
      stop("'quality' must either be constant or have the same length as 'x'")
	
    alphabetLength <-
      switch(class(x),
             DNAStringSet =, RNAStringSet = 4L,
             AAStringSet = 20L,
             256L)

    if (is(quality, "PhredQuality"))
      qualityLookupTable <- buildLookupTable(33:(33 + 99), 0:99)
    else
      qualityLookupTable <- buildLookupTable(59:(59 + 104), 0:104)
    qualityMatrices <-
      qualitySubstitutionMatrices(alphabetLength = alphabetLength,
                                  qualityClass = class(quality))

    constantLookupTable <- integer(0)
    constantMatrix <- matrix(numeric(0), nrow = 0, ncol = 0)
  } else {
    useQuality <- FALSE
    quality <- BStringSet("")
    qualityLookupTable <- integer(0)
    qualityMatrices <-
      list(match = matrix(numeric(0), nrow = 0, ncol = 0),
           mismatch = matrix(numeric(0), nrow = 0, ncol = 0))
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
    if (!isSymmetric(substitutionMatrix))
      stop("'substitutionMatrix' must be a symmetric matrix")
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
  answer <- .Call("XStringSet_align_distance",
                  x,
                  quality,
                  type,
                  typeCode,
                  gapOpening,
                  gapExtension,
                  useQuality,
                  qualityLookupTable,
                  qualityMatrices[["match"]],
                  qualityMatrices[["mismatch"]],
                  dim(qualityMatrices[["match"]]),
                  constantLookupTable,
                  constantMatrix,
                  dim(constantMatrix),
                  PACKAGE="Biostrings")
  if (method == "levenshtein")
    answer <- -answer
  attr(answer, "Size") <- length(x)
  attr(answer, "Labels") <- names(x)
  attr(answer, "Diag") <- diag
  attr(answer, "Upper") <- upper
  attr(answer, "method") <- method
  class(answer) <- "dist"
  return(answer)
}


setGeneric("stringDist", signature = "x",
           function(x, method = "levenshtein", ignoreCase = FALSE, diag = FALSE, upper = FALSE, ...)
           standardGeneric("stringDist"))

setMethod("stringDist",
          signature(x = "character"),
          function(x, method = "levenshtein", ignoreCase = FALSE, diag = FALSE, upper = FALSE,
                   type = "global", quality = PhredQuality(22L), substitutionMatrix = NULL,
                   gapOpening = 0, gapExtension = -1)
          XStringSet.stringDist(x = BStringSet(x),
                                method = method,
                                ignoreCase = ignoreCase,
                                diag = diag,
                                upper = upper,
                                type = type,
                                quality = quality,
                                substitutionMatrix = substitutionMatrix,
                                gapExtension = gapExtension,
                                gapOpening = gapOpening))

setMethod("stringDist",
          signature(x = "XStringSet"),
          function(x, method = "levenshtein", ignoreCase = FALSE, diag = FALSE, upper = FALSE,
                   type = "global", quality = PhredQuality(22L), substitutionMatrix = NULL,
                   gapOpening = 0, gapExtension = -1)
          XStringSet.stringDist(x = x,
                                method = method,
								ignoreCase = ignoreCase,
								diag = diag,
								upper = upper,
								type = type,
								quality = quality,
								substitutionMatrix = substitutionMatrix,
								gapExtension = gapExtension,
								gapOpening = gapOpening))
