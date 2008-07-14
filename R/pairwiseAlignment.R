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
function(alphabetLength = 4L, qualityType = "Phred", bitScale = 1) {
  qualityType <- match.arg(qualityType, c("Phred", "Solexa"))
  if (qualityType == "Phred") {
    errorProbs <- 10^seq(0, -9.9, by = -0.1)
    qualityLabels <- as.character(0:99)
  } else {
    errorProbs <- 1 - 1/(1 + 10^seq(0.5, -9.9, by = -0.1))
    qualityLabels <- as.character(-5:99)
  }
  errorMatrix <-
    outer(errorProbs, errorProbs,
          function(e1,e2,n) e1 + e2 - (n/(n - 1)) * e1 * e2,
          n = alphabetLength)
  output <- 
    lapply(list(matchMatrix = bitScale * log2((1 - errorMatrix) * alphabetLength),
                mismatchMatrix =
                bitScale * log2(errorMatrix * (alphabetLength / (alphabetLength - 1)))),
           function(x) {dimnames(x) <- list(qualityLabels, qualityLabels); x})
  output
}


XStringSet.pairwiseAlignment <-
function(pattern,
         subject,
         patternQuality = 22L,
         subjectQuality = 22L,
         type = "global",
         qualityType = "Phred",
         substitutionMatrix = NULL,
         gapOpening = -10,
         gapExtension = -4,
         scoreOnly = FALSE)
{
  ## Check arguments
  if (class(super(pattern)) != class(subject))
    stop("'pattern' and 'subject' must store the same underlying string type")
  type <-
    match.arg(type,
              c("global", "local", "overlap",
                "patternOverlap", "subjectOverlap"))
  typeCode <-
    c("global" = 1L, "local" = 2L, "overlap" = 3L,
      "patternOverlap" = 4L, "subjectOverlap" = 5L)[[type]]
  qualityType <- match.arg(qualityType, c("Phred", "Solexa"))
  gapOpening <- as.double(- abs(gapOpening))
  if (length(gapOpening) != 1 || is.na(gapOpening))
    stop("'gapOpening' must be a non-positive numeric vector of length 1")
  gapExtension <- as.double(- abs(gapExtension))
  if (length(gapExtension) != 1 || is.na(gapExtension))
    stop("'gapExtension' must be a non-positive numeric vector of length 1")
  scoreOnly <- as.logical(scoreOnly)
  if (length(scoreOnly) != 1 || any(is.na(scoreOnly)))
    stop("'scoreOnly' must be a non-missing logical value")

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
    useQuality <- TRUE
    if (is.numeric(patternQuality)) {
      if (qualityType == "Phred") {
        if (any(is.na(patternQuality)) ||
            any(patternQuality < 0 || patternQuality > 99))
          stop("integer 'patternQuality' values must be between 0 and 99 for qualityType 'Phred'")
        patternQuality <- rawToChar(as.raw(33L + as.integer(patternQuality)))
      } else {
        if (any(is.na(patternQuality)) ||
            any(patternQuality < -5 || patternQuality > 99))
          stop("integer 'patternQuality' values must be between -5 and 99 for qualityType 'Solexa'")
        patternQuality <- rawToChar(as.raw(64L + as.integer(patternQuality)))
      }
    }
    if (!is(patternQuality, "XStringSet"))
      patternQuality <- BStringSet(patternQuality)
    if (!(length(super(patternQuality)) %in% c(1, length(super(pattern)))))
      stop(paste("'patternQuality' must either be constant or",
                 "have the same length as 'pattern'"))

    if (is.numeric(subjectQuality)) {
      if (qualityType == "Phred") {
        if (any(is.na(subjectQuality)) ||
            any(subjectQuality < 0 || subjectQuality > 99))
          stop("integer 'subjectQuality' values must be between 0 and 99 for qualityType 'Phred'")
        subjectQuality <- rawToChar(as.raw(33L + as.integer(subjectQuality)))
      } else {
        if (any(is.na(subjectQuality)) ||
            any(subjectQuality < -5 || subjectQuality > 99))
          stop("integer 'subjectQuality' values must be between -5 and 99 for qualityType 'Solexa'")
        subjectQuality <- rawToChar(as.raw(64L + as.integer(subjectQuality)))
      }
    }
    if (!is(subjectQuality, "XString"))
      subjectQuality <- BString(subjectQuality)
    if (!(length(subjectQuality) %in% c(1, length(subject))))
      stop(paste("'subjectQuality' must either be constant or",
                 "have the same length as 'subject'"))
	
    alphabetLength <-
      switch(class(pattern),
             DNAStringSet =, RNAStringSet = 4L,
             AAStringSet = 20L,
             length(alphabetToCodes))

    if (qualityType == "Phred")
      qualityLookupTable <- buildLookupTable(33:(33 + 99), 0:99)
    else
      qualityLookupTable <- buildLookupTable(59:(59 + 104), 0:104)
    qualityMatrices <-
      qualitySubstitutionMatrices(alphabetLength = alphabetLength,
                                  qualityType = qualityType)

    constantLookupTable <- integer(0)
    constantMatrix <- matrix(numeric(0), nrow = 0, ncol = 0)
  } else {
    useQuality <- FALSE
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
                  type,
                  typeCode,
                  qualityType,
                  scoreOnly,
                  gapOpening,
                  gapExtension,
                  useQuality,
                  qualityLookupTable,
                  qualityMatrices[["matchMatrix"]],
                  qualityMatrices[["mismatchMatrix"]],
                  dim(qualityMatrices[["matchMatrix"]]),
                  constantLookupTable,
                  constantMatrix,
                  dim(constantMatrix),
                  PACKAGE="Biostrings")
  if (!scoreOnly) {
    answer@subject@unaligned <- XStringSet(class(answer@subject@unaligned), answer@subject@unaligned)
    if (useQuality)
      answer@subject@quality <- XStringSet(class(answer@subject@quality), answer@subject@quality)
  }
  return(answer)
}


setGeneric("pairwiseAlignment", signature = c("pattern", "subject"),
           function(pattern, subject, ...)
           standardGeneric("pairwiseAlignment"))

setMethod("pairwiseAlignment",
          signature(pattern = "character", subject = "character"),
          function(pattern, subject, patternQuality = 22L, subjectQuality = 22L,
                   type = "global", qualityType = "Phred", substitutionMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE)
          XStringSet.pairwiseAlignment(pattern = BStringSet(pattern),
                                       subject = BString(subject),
                                       patternQuality = patternQuality,
                                       subjectQuality = subjectQuality,
                                       type = type,
                                       qualityType = qualityType,
                                       substitutionMatrix = substitutionMatrix,
                                       gapExtension = gapExtension,
                                       gapOpening = gapOpening,
                                       scoreOnly = scoreOnly))

setMethod("pairwiseAlignment",
          signature(pattern = "character", subject = "XString"),
          function(pattern, subject, patternQuality = 22L, subjectQuality = 22L,
                   type = "global", qualityType = "Phred", substitutionMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE)
          XStringSet.pairwiseAlignment(pattern = XStringSet(class(subject), pattern),
                                       subject = subject,
                                       patternQuality = patternQuality,
                                       subjectQuality = subjectQuality,
                                       type = type,
                                       qualityType = qualityType,
                                       substitutionMatrix = substitutionMatrix,
                                       gapExtension = gapExtension,
                                       gapOpening = gapOpening,
                                       scoreOnly = scoreOnly))

setMethod("pairwiseAlignment",
          signature(pattern = "XString", subject = "character"),
          function(pattern, subject, patternQuality = 22L, subjectQuality = 22L,
                   type = "global", qualityType = "Phred", substitutionMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE)
          XStringSet.pairwiseAlignment(pattern = XStringSet(class(pattern), pattern),
                                       subject = XString(class(pattern), subject),
                                       patternQuality = patternQuality,
                                       subjectQuality = subjectQuality,
                                       type = type,
                                       qualityType = qualityType,
                                       substitutionMatrix = substitutionMatrix,
                                       gapExtension = gapExtension,
                                       gapOpening = gapOpening,
                                       scoreOnly = scoreOnly))

setMethod("pairwiseAlignment",
          signature(pattern = "XString", subject = "XString"),
          function(pattern, subject, patternQuality = 22L, subjectQuality = 22L,
                   type = "global", qualityType = "Phred", substitutionMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE)
          XStringSet.pairwiseAlignment(pattern = XStringSet(class(pattern), pattern),
                                       subject = subject,
                                       patternQuality = patternQuality,
                                       subjectQuality = subjectQuality,
                                       type = type,
                                       qualityType = qualityType,
                                       substitutionMatrix = substitutionMatrix,
                                       gapExtension = gapExtension,
                                       gapOpening = gapOpening,
                                       scoreOnly = scoreOnly))

setMethod("pairwiseAlignment",
          signature(pattern = "XStringSet", subject = "character"),
          function(pattern, subject, patternQuality = 22L, subjectQuality = 22L,
                   type = "global", qualityType = "Phred", substitutionMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE)
          XStringSet.pairwiseAlignment(pattern = pattern,
                                       subject =
                                       XString(class(super(pattern)), subject),
                                       patternQuality = patternQuality,
                                       subjectQuality = subjectQuality,
                                       type = type,
                                       qualityType = qualityType,
                                       substitutionMatrix = substitutionMatrix,
                                       gapExtension = gapExtension,
                                       gapOpening = gapOpening,
                                       scoreOnly = scoreOnly))

setMethod("pairwiseAlignment",
          signature(pattern = "XStringSet", subject = "XString"),
          function(pattern, subject, patternQuality = 22L, subjectQuality = 22L,
                   type = "global", qualityType = "Phred", substitutionMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE)
          XStringSet.pairwiseAlignment(pattern = pattern,
                                       subject = subject,
                                       patternQuality = patternQuality,
                                       subjectQuality = subjectQuality,
                                       type = type,
                                       qualityType = qualityType,
                                       substitutionMatrix = substitutionMatrix,
                                       gapExtension = gapExtension,
                                       gapOpening = gapOpening,
                                       scoreOnly = scoreOnly))
