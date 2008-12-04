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
function(errorProbability, fuzzyMatch = c(0, 1), alphabetLength = 4L, bitScale = 1) {
  if (!is.numeric(errorProbability) || !all(!is.na(errorProbability) & errorProbability >= 0 & errorProbability <= 1))
    stop("'errorProbability' must be a numeric vector with values between 0 and 1 inclusive")
  if (!is.numeric(fuzzyMatch) || !all(!is.na(fuzzyMatch) & fuzzyMatch >= 0 & fuzzyMatch <= 1))
    stop("'fuzzyMatch' must be a numeric vector with values between 0 and 1 inclusive")
  errorMatrix <-
    outer(errorProbability, errorProbability,
          function(e1,e2,n) e1 + e2 - (n/(n - 1)) * e1 * e2,
          n = alphabetLength)
  adjMatchProbs <-
    lapply(list(match = (1 - errorMatrix) * alphabetLength,
                mismatch = errorMatrix * (alphabetLength / (alphabetLength - 1))),
                function(x) {dimnames(x) <- list(names(errorProbability), names(errorProbability)); x})
  output <-
    array(NA_real_, dim = c(length(errorProbability), length(errorProbability), length(fuzzyMatch)),
          dimnames = list(names(errorProbability), names(errorProbability), as.character(fuzzyMatch)))
  for (i in seq_len(length(fuzzyMatch))) {
    output[,,i] <-
      bitScale *
        log2(fuzzyMatch[i] * adjMatchProbs[["match"]] + (1 - fuzzyMatch[i]) * adjMatchProbs[["mismatch"]])
  }
  output
}

qualitySubstitutionMatrices <-
function(fuzzyMatch = c(0, 1), alphabetLength = 4L, qualityClass = "PhredQuality", bitScale = 1) {
  if (!is.numeric(fuzzyMatch) || !all(!is.na(fuzzyMatch) & fuzzyMatch >= 0 & fuzzyMatch <= 1))
    stop("'fuzzyMatch' must be a numeric vector with values between 0 and 1 inclusive")
  if (!is(new(qualityClass), "XStringQuality"))
    stop("'qualityClass' must be one of the 'XStringQuality' classes")
  qualityIntegers <- minQuality(new(qualityClass)):maxQuality(new(qualityClass))
  errorProbability <- qualityConverter(qualityIntegers, qualityClass, "numeric")
  names(errorProbability) <- as.character(qualityIntegers)
  errorSubstitutionMatrices(errorProbability, fuzzyMatch, alphabetLength = alphabetLength, bitScale = bitScale)
}


XStringSet.pairwiseAlignment <-
function(pattern,
         subject,
         type = "global",
         substitutionMatrix = NULL,
         gapOpening = -10,
         gapExtension = -4,
         scoreOnly = FALSE)
{
  ## Check arguments
  if (baseXStringSubtype(pattern) != baseXStringSubtype(subject))
    stop("'pattern' and 'subject' must have the same XString base subtype")
  if (length(subject) != 1)
    stop("'subject' must be of length 1")
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
    unique_letters <- unique(c(uniqueLetters(pattern), uniqueLetters(subject)))
    #Even if safeLettersToInt() will deal properly with embedded nuls, I
    #suspect bad things will happen downstream in case there are any.
    alphabetToCodes <- safeLettersToInt(unique_letters, letters.as.names=TRUE)
  } else {
    alphabetToCodes <- codes(pattern)
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
  fuzzyMatrix <-
    matrix(0L, length(availableLetters), length(availableLetters),
           dimnames = list(availableLetters, availableLetters))
  diag(fuzzyMatrix) <- 1L
  fuzzyLookupTable <-
    buildLookupTable(alphabetToCodes[availableLetters], 0:(length(availableLetters) - 1))

  .Call("XStringSet_align_pairwiseAlignment",
        pattern,
        subject,
        type,
        typeCode,
        scoreOnly,
        gapOpening,
        gapExtension,
        useQuality,
        substitutionArray,
        dim(substitutionArray),
        substitutionLookupTable,
        fuzzyMatrix,
        dim(fuzzyMatrix),
        fuzzyLookupTable,
        PACKAGE="Biostrings")
}


QualityScaledXStringSet.pairwiseAlignment <-
function(pattern,
         subject,
         type = "global",
         fuzzyMatrix = NULL,
         gapOpening = -10,
         gapExtension = -4,
         scoreOnly = FALSE)
{
  ## Check arguments
  if (class(pattern) != class(subject))
    stop("'pattern' and 'subject' must be of the same class")
  if (length(subject) != 1)
    stop("'subject' must be of length 1")
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
  if (class(quality(pattern)) != class(quality(subject)))
    stop("'quality(pattern)' and 'quality(subject)' must be of the same class")

  ## Process string information
  if (is.null(codec(pattern))) {
    unique_letters <- unique(c(uniqueLetters(pattern), uniqueLetters(subject)))
    #Even if safeLettersToInt() will deal properly with embedded nuls, I
    #suspect bad things will happen downstream in case there are any.
    alphabetToCodes <- safeLettersToInt(unique_letters, letters.as.names=TRUE)
  } else {
    alphabetToCodes <- codes(pattern)
  }

  useQuality <- TRUE
  if (is.null(fuzzyMatrix)) {
      fuzzyMatrix <- diag(length(alphabetToCodes))
      dimnames(fuzzyMatrix) <- list(names(alphabetToCodes), names(alphabetToCodes))
  } else {
    if (!is.matrix(fuzzyMatrix) || !is.numeric(fuzzyMatrix) ||
        any(is.na(fuzzyMatrix)) || any(fuzzyMatrix < 0) || any(fuzzyMatrix > 1))
      stop("'fuzzyMatrix' must be a numeric matrix with values between 0 and 1 inclusive")
    if (!identical(rownames(fuzzyMatrix), colnames(fuzzyMatrix)))
      stop("row and column names differ for matrix 'fuzzyMatrix'")
    if (is.null(rownames(fuzzyMatrix)))
      stop("matrix 'fuzzyMatrix' must have row and column names")
    if (any(duplicated(rownames(fuzzyMatrix))))
      stop("matrix 'fuzzyMatrix' has duplicated row names")    
  }
  availableLetters <-
    intersect(names(alphabetToCodes), rownames(fuzzyMatrix))

  fuzzyMatrix <- fuzzyMatrix[availableLetters, availableLetters, drop = FALSE]
  uniqueFuzzyValues <- sort(unique(fuzzyMatrix))
  fuzzyReferenceMatrix <-
    matrix(match(fuzzyMatrix, uniqueFuzzyValues) - 1L,
           nrow = nrow(fuzzyMatrix), ncol = ncol(fuzzyMatrix),
           dimnames = dimnames(fuzzyMatrix))
  fuzzyLookupTable <-
    buildLookupTable(alphabetToCodes[availableLetters], 0:(length(availableLetters) - 1))

  alphabetLength <-
    switch(class(pattern),
           QualityScaledDNAStringSet =, QualityScaledRNAStringSet = 4L,
           QualityScaledAAStringSet = 20L,
           256L)

  substitutionArray <-
    qualitySubstitutionMatrices(fuzzyMatch = uniqueFuzzyValues,
                                alphabetLength = alphabetLength,
                                qualityClass = class(quality(pattern)))
  substitutionLookupTable <-
    buildLookupTable((minQuality(quality(pattern)) + offset(quality(pattern))):
                     (maxQuality(quality(pattern)) + offset(quality(pattern))),
                     0:(maxQuality(quality(pattern)) - minQuality(quality(pattern))))

  .Call("XStringSet_align_pairwiseAlignment",
        pattern,
        subject,
        type,
        typeCode,
        scoreOnly,
        gapOpening,
        gapExtension,
        useQuality,
        substitutionArray,
        dim(substitutionArray),
        substitutionLookupTable,
        fuzzyReferenceMatrix,
        dim(fuzzyReferenceMatrix),
        fuzzyLookupTable,
        PACKAGE="Biostrings")
}

mpi.collate.pairwiseAlignment <-
function(mpiOutput, pattern, subject) {
  value <- mpiOutput[[1]]
  value@score <- unlist(lapply(mpiOutput, score))

  value@pattern@unaligned <- pattern
  value@pattern@range <- do.call("c", lapply(mpiOutput, function(x) x@pattern@range))
  value@pattern@mismatch <- do.call("c", lapply(mpiOutput, function(x) x@pattern@mismatch))
  value@pattern@indel <- do.call("c", lapply(mpiOutput, function(x) x@pattern@indel))

  value@subject@unaligned <- subject
  value@subject@range <- do.call("c", lapply(mpiOutput, function(x) x@subject@range))
  value@subject@mismatch <- do.call("c", lapply(mpiOutput, function(x) x@subject@mismatch))
  value@subject@indel <- do.call("c", lapply(mpiOutput, function(x) x@subject@indel))

  value
}

mpi.XStringSet.pairwiseAlignment <-
function(pattern,
         subject,
         type = "global",
         substitutionMatrix = NULL,
         gapOpening = -10,
         gapExtension = -4,
         scoreOnly = FALSE)
{
  n <- length(pattern)
  if (n > 1 && is.loaded("mpi_comm_size")) {
      ## 'get()' are to quieten R CMD check, and for no other reason
      mpi.comm.size <- get("mpi.comm.size", mode="function")
      mpi.remote.exec <- get("mpi.remote.exec", mode="function")
      mpi.parLapply <- get("mpi.parLapply", mode="function")

      k <- min(mpi.comm.size() - 1, n)
      useMpi <- (k > 1)
  } else {
      useMpi <- FALSE
  }
  if (useMpi) {
    perNode <- n %/% k
    indices <- vector("list", k)
    for (i in seq_len(k-1)) {
        indices[[i]] <- ((i-1)*perNode+1):(i*perNode)
    }
    indices[[k]] <- ((k-1)*perNode+1):n

    mpi.remote.exec(library(Biostrings), ret = FALSE)
    mpiOutput <-
      mpi.parLapply(indices,
          function(i, pattern, subject,
                   type = "global",
                   substitutionMatrix = NULL,
                   gapOpening = -10,
                   gapExtension = -4,
                   scoreOnly = FALSE)
          Biostrings:::XStringSet.pairwiseAlignment(pattern = pattern[i],
                                       subject = subject,
                                       type = type,
                                       substitutionMatrix = substitutionMatrix,
                                       gapOpening = gapOpening,
                                       gapExtension = gapExtension,
                                       scoreOnly = scoreOnly),
          pattern = pattern, subject = subject,
          type = type,
          substitutionMatrix = substitutionMatrix,
          gapOpening = gapOpening,
          gapExtension = gapExtension,
          scoreOnly = scoreOnly)
    if (scoreOnly) {
      value <- unlist(mpiOutput)
    } else {
      value <- mpi.collate.pairwiseAlignment(mpiOutput, pattern, subject)
    }
  } else {
    value <- 
      XStringSet.pairwiseAlignment(pattern = pattern,
                                   subject = subject,
                                   type = type,
                                   substitutionMatrix = substitutionMatrix,
                                   gapOpening = gapOpening,
                                   gapExtension = gapExtension,
                                   scoreOnly = scoreOnly)
  }
  value
}

mpi.QualityScaledXStringSet.pairwiseAlignment <-
function(pattern,
         subject,
         type = "global",
         fuzzyMatrix = NULL,
         gapOpening = -10,
         gapExtension = -4,
         scoreOnly = FALSE)
{
  n <- length(pattern)
  if (n > 1 && is.loaded("mpi_comm_size")) {
    ## 'get()' are to quieten R CMD check, and for no other reason
    mpi.comm.size <- get("mpi.comm.size", mode="function")
    mpi.remote.exec <- get("mpi.remote.exec", mode="function")
    mpi.parLapply <- get("mpi.parLapply", mode="function")

    k <- min(mpi.comm.size() - 1, n)
    useMpi <- (k > 1)
  } else {
    useMpi <- FALSE
  }
  if (useMpi) {
    perNode <- n %/% k
    indices <- vector("list", k)
    for (i in seq_len(k-1)) {
      indices[[i]] <- ((i-1)*perNode+1):(i*perNode)
    }
    indices[[k]] <- ((k-1)*perNode+1):n

    mpi.remote.exec(library(Biostrings), ret = FALSE)
    mpiOutput <-
      mpi.parLapply(indices,
                    function(i, pattern, subject,
                             type = "global",
                             fuzzyMatrix = NULL,
                             gapOpening = -10,
                             gapExtension = -4,
                             scoreOnly = FALSE)
                    Biostrings:::QualityScaledXStringSet.pairwiseAlignment(pattern = pattern[i],
                                       subject = subject,
                                       type = type,
                                       fuzzyMatrix = fuzzyMatrix,
                                       gapOpening = gapOpening,
                                       gapExtension = gapExtension,
                                       scoreOnly = scoreOnly),
                    pattern = pattern, subject = subject,
                    type = type,
                    fuzzyMatrix = fuzzyMatrix,
                    gapOpening = gapOpening,
                    gapExtension = gapExtension,
                    scoreOnly = scoreOnly)
    if (scoreOnly) {
      value <- unlist(mpiOutput)
    } else {
      value <- mpi.collate.pairwiseAlignment(mpiOutput, pattern, subject)
    }
  } else {
    value <- 
      QualityScaledXStringSet.pairwiseAlignment(pattern = pattern,
                                                subject = subject,
                                                type = type,
                                                fuzzyMatrix = fuzzyMatrix,
                                                gapOpening = gapOpening,
                                                gapExtension = gapExtension,
                                                scoreOnly = scoreOnly)
  }
  value
}

setGeneric("pairwiseAlignment", signature = c("pattern", "subject"),
           function(pattern, subject, ...)
           standardGeneric("pairwiseAlignment"))

setMethod("pairwiseAlignment",
          signature(pattern = "character", subject = "character"),
          function(pattern, subject, patternQuality = PhredQuality(22L),
                   subjectQuality = PhredQuality(22L), type = "global",
                   substitutionMatrix = NULL, fuzzyMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE) {
            if (!is.null(substitutionMatrix)) {
              mpi.XStringSet.pairwiseAlignment(pattern = BStringSet(pattern),
                                               subject = BStringSet(subject),
                                               type = type,
                                               substitutionMatrix = substitutionMatrix,
                                               gapExtension = gapExtension,
                                               gapOpening = gapOpening,
                                               scoreOnly = scoreOnly)
            } else {
              mpi.QualityScaledXStringSet.pairwiseAlignment(
                  pattern = QualityScaledBStringSet(BStringSet(pattern), patternQuality),
                  subject = QualityScaledBStringSet(BStringSet(subject), subjectQuality),
                  type = type,
                  fuzzyMatrix = fuzzyMatrix,
                  gapExtension = gapExtension,
                  gapOpening = gapOpening,
                  scoreOnly = scoreOnly)
          }})

setMethod("pairwiseAlignment",
          signature(pattern = "character", subject = "XString"),
          function(pattern, subject, patternQuality = PhredQuality(22L),
                   subjectQuality = PhredQuality(22L), type = "global",
                   substitutionMatrix = NULL, fuzzyMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE) {
            if (!is.null(substitutionMatrix)) {
              mpi.XStringSet.pairwiseAlignment(pattern = XStringSet(class(subject), pattern),
                                               subject = XStringSet(class(subject), subject),
                                               type = type,
                                               substitutionMatrix = substitutionMatrix,
                                               gapExtension = gapExtension,
                                               gapOpening = gapOpening,
                                               scoreOnly = scoreOnly)
            } else {
              mpi.QualityScaledXStringSet.pairwiseAlignment(
                  pattern = QualityScaledXStringSet(XStringSet(class(subject), pattern), patternQuality),
                  subject = QualityScaledXStringSet(XStringSet(class(subject), subject), subjectQuality),
                  type = type,
                  fuzzyMatrix = fuzzyMatrix,
                  gapExtension = gapExtension,
                  gapOpening = gapOpening,
                  scoreOnly = scoreOnly)
            }})

setMethod("pairwiseAlignment",
          signature(pattern = "character", subject = "XStringSet"),
          function(pattern, subject, patternQuality = PhredQuality(22L),
                   subjectQuality = PhredQuality(22L), type = "global",
                   substitutionMatrix = NULL, fuzzyMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE) {
            if (!is.null(substitutionMatrix)) {
              mpi.XStringSet.pairwiseAlignment(pattern = XStringSet(baseXStringSubtype(subject), pattern),
                                               subject = subject,
                                               type = type,
                                               substitutionMatrix = substitutionMatrix,
                                               gapExtension = gapExtension,
                                               gapOpening = gapOpening,
                                               scoreOnly = scoreOnly)
            } else {
              mpi.QualityScaledXStringSet.pairwiseAlignment(
                  pattern = QualityScaledXStringSet(XStringSet(baseXStringSubtype(subject), pattern), patternQuality),
                  subject = QualityScaledXStringSet(subject, subjectQuality),
                  type = type,
                  fuzzyMatrix = fuzzyMatrix,
                  gapExtension = gapExtension,
                  gapOpening = gapOpening,
                  scoreOnly = scoreOnly)
          }})

setMethod("pairwiseAlignment",
          signature(pattern = "character", subject = "QualityScaledXStringSet"),
          function(pattern, subject, patternQuality = PhredQuality(22L), type = "global",
                   substitutionMatrix = NULL, fuzzyMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE) {
            if (!is.null(substitutionMatrix)) {
              mpi.XStringSet.pairwiseAlignment(pattern = XStringSet(baseXStringSubtype(subject), pattern),
                                               subject = as(subject, "XStringSet"),
                                               type = type,
                                               substitutionMatrix = substitutionMatrix,
                                               gapExtension = gapExtension,
                                               gapOpening = gapOpening,
                                               scoreOnly = scoreOnly)
            } else {
              mpi.QualityScaledXStringSet.pairwiseAlignment(
                  pattern = QualityScaledXStringSet(XStringSet(baseXStringSubtype(subject), pattern), patternQuality),
                  subject = subject,
                  type = type,
                  fuzzyMatrix = fuzzyMatrix,
                  gapExtension = gapExtension,
                  gapOpening = gapOpening,
                  scoreOnly = scoreOnly)
            }})

setMethod("pairwiseAlignment",
          signature(pattern = "XString", subject = "character"),
          function(pattern, subject, patternQuality = PhredQuality(22L),
                   subjectQuality = PhredQuality(22L), type = "global",
                   substitutionMatrix = NULL, fuzzyMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE) {
            if (!is.null(substitutionMatrix)) {
              mpi.XStringSet.pairwiseAlignment(pattern = XStringSet(class(pattern), pattern),
                                               subject = XStringSet(class(pattern), subject),
                                               type = type,
                                               substitutionMatrix = substitutionMatrix,
                                               gapExtension = gapExtension,
                                               gapOpening = gapOpening,
                                               scoreOnly = scoreOnly)
            } else {
              mpi.QualityScaledXStringSet.pairwiseAlignment(
                  pattern = QualityScaledXStringSet(XStringSet(class(pattern), pattern), patternQuality),
                  subject = QualityScaledXStringSet(XStringSet(class(pattern), subject), subjectQuality),
                  type = type,
                  fuzzyMatrix = fuzzyMatrix,
                  gapExtension = gapExtension,
                  gapOpening = gapOpening,
                  scoreOnly = scoreOnly)
          }})

setMethod("pairwiseAlignment",
          signature(pattern = "XString", subject = "XString"),
          function(pattern, subject, patternQuality = PhredQuality(22L),
                   subjectQuality = PhredQuality(22L), type = "global",
                   substitutionMatrix = NULL, fuzzyMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE) {
            if (!is.null(substitutionMatrix)) {
              mpi.XStringSet.pairwiseAlignment(pattern = XStringSet(class(pattern), pattern),
                                               subject = XStringSet(class(subject), subject),
                                               type = type,
                                               substitutionMatrix = substitutionMatrix,
                                               gapExtension = gapExtension,
                                               gapOpening = gapOpening,
                                               scoreOnly = scoreOnly)
            } else {
              mpi.QualityScaledXStringSet.pairwiseAlignment(
                  pattern = QualityScaledXStringSet(XStringSet(class(pattern), pattern), patternQuality),
                  subject = QualityScaledXStringSet(XStringSet(class(subject), subject), subjectQuality),
                  type = type,
                  fuzzyMatrix = fuzzyMatrix,
                  gapExtension = gapExtension,
                  gapOpening = gapOpening,
                  scoreOnly = scoreOnly)
            }})

setMethod("pairwiseAlignment",
          signature(pattern = "XString", subject = "XStringSet"),
          function(pattern, subject, patternQuality = PhredQuality(22L),
                   subjectQuality = PhredQuality(22L), type = "global",
                   substitutionMatrix = NULL, fuzzyMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE) {
            if (!is.null(substitutionMatrix)) {
              mpi.XStringSet.pairwiseAlignment(pattern = XStringSet(class(pattern), pattern),
                                               subject = subject,
                                               type = type,
                                               substitutionMatrix = substitutionMatrix,
                                               gapExtension = gapExtension,
                                               gapOpening = gapOpening,
                                               scoreOnly = scoreOnly)
            } else {
              mpi.QualityScaledXStringSet.pairwiseAlignment(
                  pattern = QualityScaledXStringSet(XStringSet(class(pattern), pattern), patternQuality),
                  subject = QualityScaledXStringSet(subject, subjectQuality),
                  type = type,
                  fuzzyMatrix = fuzzyMatrix,
                  gapExtension = gapExtension,
                  gapOpening = gapOpening,
                  scoreOnly = scoreOnly)
            }})

setMethod("pairwiseAlignment",
          signature(pattern = "XString", subject = "QualityScaledXStringSet"),
          function(pattern, subject, patternQuality = PhredQuality(22L), type = "global",
                   substitutionMatrix = NULL, fuzzyMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE) {
            if (!is.null(substitutionMatrix)) {
              mpi.XStringSet.pairwiseAlignment(pattern = XStringSet(class(pattern), pattern),
                                               subject = as(subject, "XStringSet"),
                                               type = type,
                                               substitutionMatrix = substitutionMatrix,
                                               gapExtension = gapExtension,
                                               gapOpening = gapOpening,
                                               scoreOnly = scoreOnly)
            } else {
              mpi.QualityScaledXStringSet.pairwiseAlignment(
                  pattern = QualityScaledXStringSet(XStringSet(class(pattern), pattern), patternQuality),
                  subject = subject,
                  type = type,
                  fuzzyMatrix = fuzzyMatrix,
                  gapExtension = gapExtension,
                  gapOpening = gapOpening,
                  scoreOnly = scoreOnly)
            }})

setMethod("pairwiseAlignment",
          signature(pattern = "XStringSet", subject = "character"),
          function(pattern, subject, patternQuality = PhredQuality(22L),
                   subjectQuality = PhredQuality(22L), type = "global",
                   substitutionMatrix = NULL, fuzzyMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE) {
          if (!is.null(substitutionMatrix)) {
            mpi.XStringSet.pairwiseAlignment(pattern = pattern,
                                             subject = XStringSet(baseXStringSubtype(pattern), subject),
                                             type = type,
                                             substitutionMatrix = substitutionMatrix,
                                             gapExtension = gapExtension,
                                             gapOpening = gapOpening,
                                             scoreOnly = scoreOnly)
          } else {
            mpi.QualityScaledXStringSet.pairwiseAlignment(
                pattern = QualityScaledXStringSet(pattern, patternQuality),
                subject = QualityScaledXStringSet(XStringSet(baseXStringSubtype(pattern), subject), subjectQuality),
                type = type,
                fuzzyMatrix = fuzzyMatrix,
                gapExtension = gapExtension,
                gapOpening = gapOpening,
                scoreOnly = scoreOnly)
          }})

setMethod("pairwiseAlignment",
          signature(pattern = "XStringSet", subject = "XString"),
          function(pattern, subject, patternQuality = PhredQuality(22L),
                   subjectQuality = PhredQuality(22L), type = "global",
                   substitutionMatrix = NULL, fuzzyMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE) {
            if (!is.null(substitutionMatrix)) {
              mpi.XStringSet.pairwiseAlignment(pattern = pattern,
                                               subject = XStringSet(class(subject), subject),
                                               type = type,
                                               substitutionMatrix = substitutionMatrix,
                                               gapExtension = gapExtension,
                                               gapOpening = gapOpening,
                                               scoreOnly = scoreOnly)
            } else {
              mpi.QualityScaledXStringSet.pairwiseAlignment(
                  pattern = QualityScaledXStringSet(pattern, patternQuality),
                  subject = QualityScaledXStringSet(XStringSet(class(subject), subject), subjectQuality),
                  type = type,
                  fuzzyMatrix = fuzzyMatrix,
                  gapExtension = gapExtension,
                  gapOpening = gapOpening,
                  scoreOnly = scoreOnly)
            }})

setMethod("pairwiseAlignment",
          signature(pattern = "XStringSet", subject = "XStringSet"),
          function(pattern, subject, patternQuality = PhredQuality(22L),
                   subjectQuality = PhredQuality(22L), type = "global",
                   substitutionMatrix = NULL, fuzzyMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE) {
            if (!is.null(substitutionMatrix)) {
              mpi.XStringSet.pairwiseAlignment(pattern = pattern,
                                               subject = subject,
                                               type = type,
                                               substitutionMatrix = substitutionMatrix,
                                               gapExtension = gapExtension,
                                               gapOpening = gapOpening,
                                               scoreOnly = scoreOnly)
            } else {
              mpi.QualityScaledXStringSet.pairwiseAlignment(
                  pattern = QualityScaledXStringSet(pattern, patternQuality),
                  subject = QualityScaledXStringSet(subject, subjectQuality),
                  type = type,
                  fuzzyMatrix = fuzzyMatrix,
                  gapExtension = gapExtension,
                  gapOpening = gapOpening,
                  scoreOnly = scoreOnly)
            }})

setMethod("pairwiseAlignment",
          signature(pattern = "XStringSet", subject = "QualityScaledXStringSet"),
          function(pattern, subject, patternQuality = PhredQuality(22L), type = "global",
                   substitutionMatrix = NULL, fuzzyMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE) {
            if (!is.null(substitutionMatrix)) {
              mpi.XStringSet.pairwiseAlignment(pattern = pattern,
                                               subject = as(subject, "XStringSet"),
                                               type = type,
                                               substitutionMatrix = substitutionMatrix,
                                               gapExtension = gapExtension,
                                               gapOpening = gapOpening,
                                               scoreOnly = scoreOnly)
            } else {
              mpi.QualityScaledXStringSet.pairwiseAlignment(
                  pattern = QualityScaledXStringSet(pattern, patternQuality),
                  subject = subject,
                  type = type,
                  fuzzyMatrix = fuzzyMatrix,
                  gapExtension = gapExtension,
                  gapOpening = gapOpening,
                  scoreOnly = scoreOnly)
            }})

setMethod("pairwiseAlignment",
          signature(pattern = "QualityScaledXStringSet", subject = "character"),
          function(pattern, subject, subjectQuality = PhredQuality(22L), type = "global",
                   substitutionMatrix = NULL, fuzzyMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE) {
            if (!is.null(substitutionMatrix)) {
              mpi.XStringSet.pairwiseAlignment(pattern = as(pattern, "XStringSet"),
                                               subject = XStringSet(baseXStringSubtype(pattern), subject),
                                               type = type,
                                               substitutionMatrix = substitutionMatrix,
                                               gapExtension = gapExtension,
                                               gapOpening = gapOpening,
                                               scoreOnly = scoreOnly)
            } else {
              mpi.QualityScaledXStringSet.pairwiseAlignment(
                  pattern = pattern,
                  subject = QualityScaledXStringSet(XStringSet(baseXStringSubtype(pattern), subject), subjectQuality),
                  type = type,
                  fuzzyMatrix = fuzzyMatrix,
                  gapExtension = gapExtension,
                  gapOpening = gapOpening,
                  scoreOnly = scoreOnly)
            }})

setMethod("pairwiseAlignment",
          signature(pattern = "QualityScaledXStringSet", subject = "XString"),
          function(pattern, subject, subjectQuality = PhredQuality(22L), type = "global",
                   substitutionMatrix = NULL, fuzzyMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE) {
            if (!is.null(substitutionMatrix)) {
              mpi.XStringSet.pairwiseAlignment(pattern = as(pattern, "XStringSet"),
                                               subject = XStringSet(class(subject), subject),
                                               type = type,
                                               substitutionMatrix = substitutionMatrix,
                                               gapExtension = gapExtension,
                                               gapOpening = gapOpening,
                                               scoreOnly = scoreOnly)
            } else {
              mpi.QualityScaledXStringSet.pairwiseAlignment(
                  pattern = pattern,
                  subject = QualityScaledXStringSet(XStringSet(class(subject), subject), subjectQuality),
                  type = type,
                  fuzzyMatrix = fuzzyMatrix,
                  gapExtension = gapExtension,
                  gapOpening = gapOpening,
                  scoreOnly = scoreOnly)
            }})

setMethod("pairwiseAlignment",
          signature(pattern = "QualityScaledXStringSet", subject = "XStringSet"),
          function(pattern, subject, subjectQuality = PhredQuality(22L), type = "global",
                   substitutionMatrix = NULL, fuzzyMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE) {
            if (!is.null(substitutionMatrix)) {
              mpi.XStringSet.pairwiseAlignment(pattern = as(pattern, "XStringSet"),
                                               subject = subject,
                                               type = type,
                                               substitutionMatrix = substitutionMatrix,
                                               gapExtension = gapExtension,
                                               gapOpening = gapOpening,
                                               scoreOnly = scoreOnly)
            } else {
              mpi.QualityScaledXStringSet.pairwiseAlignment(
                  pattern = pattern,
                  subject = QualityScaledXStringSet(subject, subjectQuality),
                  type = type,
                  fuzzyMatrix = fuzzyMatrix,
                  gapExtension = gapExtension,
                  gapOpening = gapOpening,
                  scoreOnly = scoreOnly)
            }})

setMethod("pairwiseAlignment",
          signature(pattern = "QualityScaledXStringSet", subject = "QualityScaledXStringSet"),
          function(pattern, subject, type = "global",
                   substitutionMatrix = NULL, fuzzyMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE) {
            if (!is.null(substitutionMatrix)) {
              mpi.XStringSet.pairwiseAlignment(pattern = as(pattern, "XStringSet"),
                                               subject = as(subject, "XStringSet"),
                                               type = type,
                                               substitutionMatrix = substitutionMatrix,
                                               gapExtension = gapExtension,
                                               gapOpening = gapOpening,
                                               scoreOnly = scoreOnly)
            } else {
              mpi.QualityScaledXStringSet.pairwiseAlignment(
                  pattern = pattern,
                  subject = subject,
                  type = type,
                  fuzzyMatrix = fuzzyMatrix,
                  gapExtension = gapExtension,
                  gapOpening = gapOpening,
                  scoreOnly = scoreOnly)
            }})
