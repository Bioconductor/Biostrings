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

qualitySubstitutionMatrices <-
function(alphabetLength = 4L, bitScale = 1/2) {
  errorProbs <- 10^seq(0, -9.9, by = -0.1)
  errorMatrix <-
    outer(errorProbs, errorProbs,
          function(e1,e2,n) e1 + e2 - (n/(n - 1)) * e1 * e2,
          n = alphabetLength)
  list(matchMatrix = bitScale * log2((1 - errorMatrix) * alphabetLength),
       mismatchMatrix =
         bitScale * log2(errorMatrix * (alphabetLength / (alphabetLength - 1))))
}


XString.pairwiseAlignment <-
function(pattern,
         subject,
         patternQuality = 22L,
         subjectQuality = 22L,
         type = "global",
         substitutionMatrix = NULL,
         gapOpening = -10,
         gapExtension = -2,
         scoreOnly = FALSE)
{
  ## Check arguments
  if (class(pattern) != class(subject))
    stop("'pattern' and 'subject' must have the same class")
  type <-
    match.arg(tolower(type),
              c("global", "local", "overlap",
                "patternOverlap", "subjectOverlap"))
  typeCode <-
    c("global" = 1L, "local" = 2L, "overlap" = 3L,
      "patternOverlap" = 4L, "subjectOverlap" = 5L)[[type]]
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
  if (is.null(codec(pattern))) {
    uniqueBases <-
      unique(c(unique(charToRaw(as.character(pattern))),
               unique(charToRaw(as.character(subject)))))
    alphabetToCodes <- as.integer(uniqueBases)
    names(alphabetToCodes) <- rawToChar(uniqueBases, multiple = TRUE)
  } else {
    stringCodec <- codec(pattern)
    alphabetToCodes <- stringCodec@codes
    names(alphabetToCodes) <- stringCodec@letters
  }

  ## Generate quality-based and constant substitution matrix information
  if (is.null(substitutionMatrix)) {
    if (is.numeric(patternQuality)) {
      if (any(is.na(patternQuality)) ||
          any(patternQuality < 0 || patternQuality > 99))
        stop("integer 'patternQuality' values must be between 0 and 99")
      patternQuality <- rawToChar(as.raw(33L + as.integer(patternQuality)))
    }
    if (!is(patternQuality, "XString"))
      patternQuality <- BString(patternQuality)
    if (!(length(patternQuality) %in% c(1, length(pattern))))
      stop(paste("'patternQuality' must either be constant or",
                 "have the same length as 'pattern'"))

    if (is.numeric(subjectQuality)) {
      if (any(is.na(subjectQuality)) ||
          any(subjectQuality < 0 || subjectQuality > 99))
        stop("integer 'subjectQuality' values must be between 0 and 99")
      subjectQuality <- rawToChar(as.raw(33L + as.integer(subjectQuality)))
    }
    if (!is(subjectQuality, "XString"))
      subjectQuality <- BString(subjectQuality)
    if (!(length(subjectQuality) %in% c(1, length(subject))))
      stop(paste("'subjectQuality' must either be constant or",
                 "have the same length as 'subject'"))

    alphabetLength <-
      switch(class(pattern),
             DNAString =, RNAString = 4L,
             AAString = 20L,
             length(alphabetToCodes))

    qualityLookupTable <- buildLookupTable(33:(99 + 33), 0:99)
    qualityMatrices <- qualitySubstitutionMatrices(alphabetLength = alphabetLength)

    constantLookupTable <- integer(0)
    constantMatrix <- matrix(numeric(0), nrow = 0, ncol = 0)
  } else {
    patternQuality <- BString("")
    subjectQuality <- BString("")
    qualityLookupTable <- integer(0)
    qualityMatrices <-
      list(matchMatrix = matrix(numeric(0), nrow = 0, ncol = 0),
           mismatchMatrix = matrix(numeric(0), nrow = 0, ncol = 0))
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
  patternProfile <- XNumeric(ifelse(scoreOnly, 0L, nchar(pattern)))
  subjectProfile <- XNumeric(ifelse(scoreOnly, 0L, nchar(subject)))
  answer <- .Call("XString_align_pairwiseAlignment",
                  pattern,
                  subject,
                  patternQuality,
                  subjectQuality,
                  patternProfile@xp,
                  subjectProfile@xp,
                  typeCode,
                  scoreOnly,
                  gapOpening,
                  gapExtension,
                  qualityLookupTable,
                  qualityMatrices[["matchMatrix"]],
                  qualityMatrices[["mismatchMatrix"]],
                  dim(qualityMatrices[["matchMatrix"]]),
                  constantLookupTable,
                  constantMatrix,
                  dim(constantMatrix),
                  PACKAGE="Biostrings")
  if (scoreOnly) {
    output <- answer
  } else {
    output <- new("PairwiseAlignment",
                  pattern =
                  new("AlignedXString",
                      unaligned = pattern,
                      quality = patternQuality,
                      range =
                      IRanges(start = answer[["startPatternRange"]],
                              width = answer[["widthPatternRange"]]),
                      inserts =
                      IRanges(start = answer[["startPatternInserts"]],
                              width = answer[["widthPatternInserts"]]),
                      profile = patternProfile),
                  subject =
                  new("AlignedXString",
                      unaligned = subject,
                      quality = subjectQuality,
                      range =
                      IRanges(start = answer[["startSubjectRange"]],
                              width = answer[["widthSubjectRange"]]),
                      inserts =
                      IRanges(start = answer[["startSubjectInserts"]],
                              width = answer[["widthSubjectInserts"]]),
                      profile = subjectProfile),
                  score = answer[["score"]],
                  type = type,
                  constantMatrix = constantMatrix,
                  gapOpening = gapOpening,
                  gapExtension = gapExtension)
  }
  return(output)
}


XStringSet.pairwiseAlignment <-
function(pattern,
         subject,
         patternQuality = 22L,
         subjectQuality = 22L,
         type = "global",
         substitutionMatrix = NULL,
         gapOpening = -10,
         gapExtension = -2,
         scoreOnly = TRUE)
{
  ## Check arguments
  if (class(super(pattern)) != class(subject))
    stop("'pattern' and 'subject' must store the same underlying string type")
  type <-
    match.arg(tolower(type),
              c("global", "local", "overlap",
                "patternOverlap", "subjectOverlap"))
  typeCode <-
    c("global" = 1L, "local" = 2L, "overlap" = 3L,
      "patternOverlap" = 4L, "subjectOverlap" = 5L)[[type]]
  gapOpening <- as.double(- abs(gapOpening))
  if (is.na(gapOpening) || length(gapOpening) != 1)
    stop("'gapOpening' must be a non-positive numeric vector of length 1")
  gapExtension <- as.double(- abs(gapExtension))
  if (is.na(gapExtension) || length(gapExtension) != 1)
    stop("'gapExtension' must be a non-positive numeric vector of length 1")
  scoreOnly <- as.logical(scoreOnly)
  if (length(scoreOnly) != 1 || any(is.na(scoreOnly)))
    stop("'scoreOnly' must be a non-missing logical value")
  if (!scoreOnly)
    stop("'scoreOnly' must be TRUE when 'pattern' is an 'XStringSet'")

  ## Process string information
  if (is.null(codec(pattern))) {
    uniqueBases <-
      unique(c(unique(charToRaw(as.character(super(pattern)))),
               unique(charToRaw(as.character(subject)))))
    alphabetToCodes <- as.integer(uniqueBases)
    names(alphabetToCodes) <- rawToChar(uniqueBases, multiple = TRUE)
  } else {
    stringCodec <- codec(pattern)
    alphabetToCodes <- stringCodec@codes
    names(alphabetToCodes) <- stringCodec@letters
  }

  ## Generate quality-based and constant substitution matrix information
  if (is.null(substitutionMatrix)) {
    if (is.numeric(patternQuality)) {
      if (any(is.na(patternQuality)) ||
          any(patternQuality < 0 || patternQuality > 99))
        stop("integer 'patternQuality' values must be between 0 and 99")
      patternQuality <- rawToChar(as.raw(33L + as.integer(patternQuality)))
    }
    if (!is(patternQuality, "XStringSet"))
      patternQuality <- BStringSet(patternQuality)
    if (!(length(super(patternQuality)) %in% c(1, length(super(pattern)))))
      stop(paste("'patternQuality' must either be constant or",
                 "have the same length as 'pattern'"))

    if (is.numeric(subjectQuality)) {
      if (any(is.na(subjectQuality)) ||
          any(subjectQuality < 0 || subjectQuality > 99))
        stop("integer 'subjectQuality' values must be between 0 and 99")
      subjectQuality <- rawToChar(as.raw(33L + as.integer(subjectQuality)))
    }
    if (!is(subjectQuality, "XString"))
      subjectQuality <- BString(subjectQuality)
    if (!(length(subjectQuality) %in% c(1, length(subject))))
      stop(paste("'subjectQuality' must either be constant or",
                 "have the same length as 'subject'"))
	
    alphabetLength <-
      switch(class(pattern),
             DNAString =, RNAString = 4L,
             AAString = 20L,
             length(alphabetToCodes))

    qualityLookupTable <- buildLookupTable(33:(99 + 33), 0:99)
    qualityMatrices <- qualitySubstitutionMatrices(alphabetLength = alphabetLength)

    constantLookupTable <- integer(0)
    constantMatrix <- matrix(numeric(0), nrow = 0, ncol = 0)
  } else {
    patternQuality <- BStringSet("")
    subjectQuality <- BString("")
    qualityLookupTable <- integer(0)
    qualityMatrices <-
      list(matchMatrix = matrix(numeric(0), nrow = 0, ncol = 0),
           mismatchMatrix = matrix(numeric(0), nrow = 0, ncol = 0))
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
  answer <- .Call("XStringSet_align_pairwiseAlignment",
                  pattern,
                  subject,
                  patternQuality,
                  subjectQuality,
                  typeCode,
                  gapOpening,
                  gapExtension,
                  qualityLookupTable,
                  qualityMatrices[["matchMatrix"]],
                  qualityMatrices[["mismatchMatrix"]],
                  dim(qualityMatrices[["matchMatrix"]]),
                  constantLookupTable,
                  constantMatrix,
                  dim(constantMatrix),
                  PACKAGE="Biostrings")
  return(answer)
}


setGeneric("pairwiseAlignment", signature = c("pattern", "subject"),
           function(pattern, subject,
                    patternQuality = 22L, subjectQuality = 22L,
                    type = "global", substitutionMatrix = NULL,
                    gapOpening = -10, gapExtension = -2,
                    scoreOnly = is(pattern, "XStringSet"))
           standardGeneric("pairwiseAlignment"))

setMethod("pairwiseAlignment",
          signature(pattern = "character", subject = "character"),
          function(pattern, subject, patternQuality = 22L, subjectQuality = 22L,
                   type = "global", substitutionMatrix = NULL,
                   gapOpening = -10, gapExtension = -2,
                   scoreOnly = FALSE)
          XString.pairwiseAlignment(pattern = BString(pattern),
                                    subject = BString(subject),
                                    patternQuality = patternQuality,
                                    subjectQuality = subjectQuality,
                                    type = type,
                                    substitutionMatrix = substitutionMatrix,
                                    gapExtension = gapExtension,
                                    gapOpening = gapOpening,
                                    scoreOnly = scoreOnly))

setMethod("pairwiseAlignment",
          signature(pattern = "character", subject = "XString"),
          function(pattern, subject, patternQuality = 22L, subjectQuality = 22L,
                   type = "global", substitutionMatrix = NULL,
                   gapOpening = -10, gapExtension = -2,
                   scoreOnly = FALSE)
          XString.pairwiseAlignment(pattern = XString(class(subject), pattern),
                                    subject = subject,
                                    patternQuality = patternQuality,
                                    subjectQuality = subjectQuality,
                                    type = type,
                                    substitutionMatrix = substitutionMatrix,
                                    gapExtension = gapExtension,
                                    gapOpening = gapOpening,
                                    scoreOnly = scoreOnly))

setMethod("pairwiseAlignment",
          signature(pattern = "XString", subject = "character"),
          function(pattern, subject, patternQuality = 22L, subjectQuality = 22L,
                   type = "global", substitutionMatrix = NULL,
                   gapOpening = -10, gapExtension = -2,
                   scoreOnly = FALSE)
          XString.pairwiseAlignment(pattern = pattern,
                                    subject = XString(class(pattern), subject),
                                    patternQuality = patternQuality,
                                    subjectQuality = subjectQuality,
                                    type = type,
                                    substitutionMatrix = substitutionMatrix,
                                    gapExtension = gapExtension,
                                    gapOpening = gapOpening,
                                    scoreOnly = scoreOnly))

setMethod("pairwiseAlignment",
          signature(pattern = "XString", subject = "XString"),
          function(pattern, subject, patternQuality = 22L, subjectQuality = 22L,
                   type = "global", substitutionMatrix = NULL,
                   gapOpening = -10, gapExtension = -2,
                   scoreOnly = FALSE)
          XString.pairwiseAlignment(pattern = pattern,
                                    subject = subject,
                                    patternQuality = patternQuality,
                                    subjectQuality = subjectQuality,
                                    type = type,
                                    substitutionMatrix = substitutionMatrix,
                                    gapExtension = gapExtension,
                                    gapOpening = gapOpening,
                                    scoreOnly = scoreOnly))

setMethod("pairwiseAlignment",
          signature(pattern = "XStringSet", subject = "character"),
          function(pattern, subject, patternQuality = 22L, subjectQuality = 22L,
                   type = "global", substitutionMatrix = NULL,
                   gapOpening = -10, gapExtension = -2,
                   scoreOnly = TRUE)
          XStringSet.pairwiseAlignment(pattern = pattern,
                                       subject =
                                       XString(class(super(pattern)), subject),
                                       patternQuality = patternQuality,
                                       subjectQuality = subjectQuality,
                                       type = type,
                                       substitutionMatrix = substitutionMatrix,
                                       gapExtension = gapExtension,
                                       gapOpening = gapOpening,
                                       scoreOnly = scoreOnly))

setMethod("pairwiseAlignment",
          signature(pattern = "XStringSet", subject = "XString"),
          function(pattern, subject, patternQuality = 22L, subjectQuality = 22L,
                   type = "global", substitutionMatrix = NULL,
                   gapOpening = -10, gapExtension = -2,
                   scoreOnly = TRUE)
          XStringSet.pairwiseAlignment(pattern = pattern,
                                       subject = subject,
                                       patternQuality = patternQuality,
                                       subjectQuality = subjectQuality,
                                       type = type,
                                       substitutionMatrix = substitutionMatrix,
                                       gapExtension = gapExtension,
                                       gapOpening = gapOpening,
                                       scoreOnly = scoreOnly))
