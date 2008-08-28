### =========================================================================
### The pairwiseAlignment() generic & related functions
### -------------------------------------------------------------------------
###
### The pairwiseAligment() function provides optimal pairwise alignment of
### the following types:
### - Global alignment
### - Local alignment
### - Overlap alignment
### - Pattern Overlap alignment
### - Subject Overlap alignment
###
### -------------------------------------------------------------------------


nucleotideSubstitutionMatrix <- function(match = 1, mismatch = 0, baseOnly = FALSE, type = "DNA")
{
  "%safemult%" <- function(x, y) ifelse(is.infinite(x) & y == 0, 0, x * y)
  type <- match.arg(type, c("DNA", "RNA"))
  if (!isSingleNumber(match) || !isSingleNumber(mismatch))
    stop("'match' and 'mismatch' must be non-missing numbers")
  if (baseOnly)
    letters <- IUPAC_CODE_MAP[c("A", "C", "G", "T")]
  else
    letters <- IUPAC_CODE_MAP
  if (type == "RNA")
    names(letters) <- chartr("T", "U", names(letters))
  nLetters <- length(letters)
  splitLetters <- strsplit(letters,split="")
  submat <- matrix(0, nrow = nLetters, ncol = nLetters, dimnames = list(names(letters), names(letters)))
  for(i in 1:nLetters)
    for(j in i:nLetters)
      submat[i,j] <- submat[j,i] <- mean(outer(splitLetters[[i]], splitLetters[[j]], "=="))
  abs(match) * submat - abs(mismatch) %safemult% (1 - submat)
}

errorSubstitutionMatrices <-
function(prob, alphabetLength = 4L, bitScale = 1) {
  errorMatrix <-
    outer(prob, prob,
          function(e1,e2,n) e1 + e2 - (n/(n - 1)) * e1 * e2,
          n = alphabetLength)
  output <- 
    lapply(list(match = bitScale * log2((1 - errorMatrix) * alphabetLength),
                mismatch =
                bitScale * log2(errorMatrix * (alphabetLength / (alphabetLength - 1)))),
                function(x) {dimnames(x) <- list(names(prob), names(prob)); x})
  output
}

qualitySubstitutionMatrices <-
function(alphabetLength = 4L, qualityClass = "PhredQuality", bitScale = 1) {
  qualityClass <- match.arg(qualityClass, c("PhredQuality", "SolexaQuality"))
  if (qualityClass == "PhredQuality") {
    prob <- 10^seq(0, -9.9, by = -0.1)
    names(prob) <- as.character(0:99)
  } else {
    prob <- 1 - 1/(1 + 10^seq(0.5, -9.9, by = -0.1))
    names(prob) <- as.character(-5:99)
  }
  errorSubstitutionMatrices(prob, alphabetLength = alphabetLength, bitScale = bitScale)
}


XStringSet.pairwiseAlignment <-
function(pattern,
         subject,
         patternQuality = PhredQuality(22L),
         subjectQuality = PhredQuality(22L),
         type = "global",
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
    if (class(patternQuality) != class(subjectQuality))
      stop("'patternQuality' and 'subjectQuality' must be of the same class")

    if (class(patternQuality) %in% c("integer", "numeric", "BString", "BStringSet"))
        patternQuality <- PhredQuality(patternQuality)
    if (!is(patternQuality, "XStringQuality"))
      stop("'patternQuality' must be of class 'XStringSet'")
    if (!all(nchar(patternQuality) == 1 | nchar(patternQuality) == nchar(pattern)))
      stop(paste("'patternQuality' must either be constant or",
                 "have the same length as 'pattern'"))

    if (class(subjectQuality) %in% c("integer", "numeric", "BString", "BStringSet"))
      subjectQuality <- PhredQuality(subjectQuality)
    if (!is(subjectQuality, "XStringQuality"))
      stop("'subjectQuality' must be of class 'XStringSet'")
    if (!(length(subjectQuality) %in% c(1, length(subject))))
      stop(paste("'subjectQuality' must either be constant or",
                 "have the same length as 'subject'"))

    alphabetLength <-
      switch(class(pattern),
             DNAStringSet =, RNAStringSet = 4L,
             AAStringSet = 20L,
             256L)

    if (is(patternQuality, "PhredQuality"))
      qualityLookupTable <- buildLookupTable(33:(33 + 99), 0:99)
    else
      qualityLookupTable <- buildLookupTable(59:(59 + 104), 0:104)
    qualityMatrices <-
      qualitySubstitutionMatrices(alphabetLength = alphabetLength,
                                  qualityClass = class(patternQuality))

    constantLookupTable <- integer(0)
    constantMatrix <- matrix(numeric(0), nrow = 0, ncol = 0)
  } else {
    useQuality <- FALSE
    patternQuality <- PhredQuality("")
    subjectQuality <- PhredQuality("")
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
                  scoreOnly,
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
  if (!scoreOnly) {
    answer@subject@unaligned <- XStringSet(class(answer@subject@unaligned), answer@subject@unaligned)
  }
  return(answer)
}


setGeneric("pairwiseAlignment", signature = c("pattern", "subject"),
           function(pattern, subject, ...)
           standardGeneric("pairwiseAlignment"))

setMethod("pairwiseAlignment",
          signature(pattern = "character", subject = "character"),
          function(pattern, subject, patternQuality = PhredQuality(22L),
                   subjectQuality = PhredQuality(22L), type = "global",
                   substitutionMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE)
          XStringSet.pairwiseAlignment(pattern = BStringSet(pattern),
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
          function(pattern, subject, patternQuality = PhredQuality(22L),
                   subjectQuality = PhredQuality(22L), type = "global",
                   substitutionMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE)
          XStringSet.pairwiseAlignment(pattern = XStringSet(class(subject), pattern),
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
          function(pattern, subject, patternQuality = PhredQuality(22L),
                   subjectQuality = PhredQuality(22L), type = "global",
                   substitutionMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE)
          XStringSet.pairwiseAlignment(pattern = XStringSet(class(pattern), pattern),
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
          function(pattern, subject, patternQuality = PhredQuality(22L),
                   subjectQuality = PhredQuality(22L), type = "global",
                   substitutionMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE)
          XStringSet.pairwiseAlignment(pattern = XStringSet(class(pattern), pattern),
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
          function(pattern, subject, patternQuality = PhredQuality(22L),
                   subjectQuality = PhredQuality(22L), type = "global",
                   substitutionMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE)
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
          function(pattern, subject, patternQuality = PhredQuality(22L),
                   subjectQuality = PhredQuality(22L), type = "global",
                   substitutionMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE)
          XStringSet.pairwiseAlignment(pattern = pattern,
                                       subject = subject,
                                       patternQuality = patternQuality,
                                       subjectQuality = subjectQuality,
                                       type = type,
                                       substitutionMatrix = substitutionMatrix,
                                       gapExtension = gapExtension,
                                       gapOpening = gapOpening,
                                       scoreOnly = scoreOnly))
