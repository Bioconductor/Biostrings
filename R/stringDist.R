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
    unique_letters <- uniqueLetters(x)
    #Even if safeLettersToInt() will deal properly with embedded nuls, I
    #suspect bad things will happen downstream in case there are any.
    alphabetToCodes <- safeLettersToInt(unique_letters, letters.as.names=TRUE)
  } else {
    alphabetToCodes <- codes(x)
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

  useQuality <- FALSE
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

  substitutionMatrix <-
    matrix(as.double(substitutionMatrix[availableLetters, availableLetters]),
           nrow = length(availableLetters),
           ncol = length(availableLetters),
           dimnames = list(availableLetters, availableLetters))
  substitutionArray <-
    array(unlist(substitutionMatrix, substitutionMatrix), dim = c(dim(substitutionMatrix), 2),
          dimnames = list(availableLetters, availableLetters, c("0", "1")))
  substitutionLookupTable <-
    buildLookupTable(alphabetToCodes[availableLetters],
                     0:(length(availableLetters) - 1))
  ambiguityMatrix <-
    matrix(0L, length(availableLetters), length(availableLetters),
           dimnames = list(availableLetters, availableLetters))
  diag(ambiguityMatrix) <- 1L
  ambiguityLookupTable <-
    buildLookupTable(alphabetToCodes[availableLetters], 0:(length(availableLetters) - 1))

  answer <- .Call("XStringSet_align_distance",
                  x,
                  type,
                  typeCode,
                  gapOpening,
                  gapExtension,
                  useQuality,
                  substitutionArray,
                  dim(substitutionArray),
                  substitutionLookupTable,
                  ambiguityMatrix,
                  dim(ambiguityMatrix),
                  ambiguityLookupTable,
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


QualityScaledXStringSet.stringDist <-
function(x,
         ignoreCase = FALSE,
         diag = FALSE,
         upper = FALSE,
         type = "global",
         ambiguityMatrix = NULL,
         gapOpening = 0,
         gapExtension = -1)
{
  ## Check arguments
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
    unique_letters <- uniqueLetters(x)
    #Even if safeLettersToInt() will deal properly with embedded nuls, I
    #suspect bad things will happen downstream in case there are any.
    alphabetToCodes <- safeLettersToInt(unique_letters, letters.as.names=TRUE)
  } else {
    alphabetToCodes <- codes(x)
  }

  useQuality <- TRUE
  if (is.null(ambiguityMatrix)) {
    ambiguityMatrix <- diag(length(alphabetToCodes))
    dimnames(ambiguityMatrix) <- list(names(alphabetToCodes), names(alphabetToCodes))
  } else {
    if (!is.matrix(ambiguityMatrix) || !is.numeric(ambiguityMatrix) ||
        any(is.na(ambiguityMatrix)) || any(ambiguityMatrix < 0) || any(ambiguityMatrix > 1))
      stop("'ambiguityMatrix' must be a numeric matrix with values between 0 and 1 inclusive")
    if (!identical(rownames(ambiguityMatrix), colnames(ambiguityMatrix)))
      stop("row and column names differ for matrix 'ambiguityMatrix'")
    if (is.null(rownames(ambiguityMatrix)))
      stop("matrix 'ambiguityMatrix' must have row and column names")
    if (any(duplicated(rownames(ambiguityMatrix))))
      stop("matrix 'ambiguityMatrix' has duplicated row names")    
  }
  availableLetters <-
    intersect(names(alphabetToCodes), rownames(ambiguityMatrix))
  
  ambiguityMatrix <- ambiguityMatrix[availableLetters, availableLetters]
  uniqueAmbiguities <- sort(unique(ambiguityMatrix))
  ambiguityReferenceMatrix <-
    matrix(match(ambiguityMatrix, uniqueAmbiguities) - 1L,
           nrow = nrow(ambiguityMatrix), ncol = ncol(ambiguityMatrix),
           dimnames = dimnames(ambiguityMatrix))
  ambiguityLookupTable <-
    buildLookupTable(alphabetToCodes[availableLetters], 0:(length(availableLetters) - 1))
  
  alphabetLength <-
    switch(class(x),
           QualityScaledDNAStringSet =, QualityScaledRNAStringSet = 4L,
           QualityScaledAAStringSet = 20L,
           256L)

  substitutionArray <-
    qualitySubstitutionMatrices(matchAmbiguity = uniqueAmbiguities,
                                alphabetLength = alphabetLength,
                                qualityClass = class(quality(x)))
  substitutionLookupTable <-
    buildLookupTable((minQuality(quality(x)) + offset(quality(x))):
                     (maxQuality(quality(x)) + offset(quality(x))),
                     0:(maxQuality(quality(x)) - minQuality(quality(x))))

  answer <- .Call("XStringSet_align_distance",
                  x,
                  type,
                  typeCode,
                  gapOpening,
                  gapExtension,
                  useQuality,
                  substitutionArray,
                  dim(substitutionArray),
                  substitutionLookupTable,
                  ambiguityReferenceMatrix,
                  dim(ambiguityReferenceMatrix),
                  ambiguityLookupTable,
                  PACKAGE="Biostrings")
  attr(answer, "Size") <- length(x)
  attr(answer, "Labels") <- names(x)
  attr(answer, "Diag") <- diag
  attr(answer, "Upper") <- upper
  attr(answer, "method") <- "quality"
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
                   ambiguityMatrix = NULL, gapOpening = 0, gapExtension = -1) {
            if (method != "quality") {
              XStringSet.stringDist(x = BStringSet(x),
                                    method = method,
                                    ignoreCase = ignoreCase,
                                    diag = diag,
                                    upper = upper,
                                    type = type,
                                    substitutionMatrix = substitutionMatrix,
                                    gapExtension = gapExtension,
                                    gapOpening = gapOpening)
            } else {
              QualityScaledXStringSet.stringDist(x = QualityScaledBStringSet(x, quality),
                                                 ignoreCase = ignoreCase,
                                                 diag = diag,
                                                 upper = upper,
                                                 type = type,
                                                 ambiguityMatrix = ambiguityMatrix,
                                                 gapExtension = gapExtension,
                                                 gapOpening = gapOpening)
          }})

setMethod("stringDist",
          signature(x = "XStringSet"),
          function(x, method = "levenshtein", ignoreCase = FALSE, diag = FALSE, upper = FALSE,
                   type = "global", quality = PhredQuality(22L), substitutionMatrix = NULL,
                   ambiguityMatrix = NULL, gapOpening = 0, gapExtension = -1) {
            if (method != "quality") {
              XStringSet.stringDist(x = x,
                                    method = method,
                                    ignoreCase = ignoreCase,
                                    diag = diag,
                                    upper = upper,
                                    type = type,
                                    substitutionMatrix = substitutionMatrix,
                                    gapExtension = gapExtension,
                                    gapOpening = gapOpening)
             } else {
               QualityScaledXStringSet.stringDist(x = QualityScaledXStringSet(x, quality),
                                                  ignoreCase = ignoreCase,
                                                  diag = diag,
                                                  upper = upper,
                                                  type = type,
                                                  ambiguityMatrix = ambiguityMatrix,
                                                  gapExtension = gapExtension,
                                                  gapOpening = gapOpening)
          }})

setMethod("stringDist",
          signature(x = "QualityScaledXStringSet"),
          function(x, method = "quality", ignoreCase = FALSE, diag = FALSE, upper = FALSE,
                   type = "global", substitutionMatrix = NULL, ambiguityMatrix = NULL,
                   gapOpening = 0, gapExtension = -1) {
            if (method != "quality") {
              XStringSet.stringDist(x = as(x, "XStringSet"),
                                   method = method,
                                   ignoreCase = ignoreCase,
                                   diag = diag,
                                   upper = upper,
                                   type = type,
                                   substitutionMatrix = substitutionMatrix,
                                   gapExtension = gapExtension,
                                   gapOpening = gapOpening)
            } else {
              QualityScaledXStringSet.stringDist(x = x,
                                                 ignoreCase = ignoreCase,
                                                 diag = diag,
                                                 upper = upper,
                                                 type = type,
                                                 ambiguityMatrix = ambiguityMatrix,
                                                 gapExtension = gapExtension,
                                                 gapOpening = gapOpening)
            }})
