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
    letters <- IUPAC_CODE_MAP[DNA_BASES]
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
         gapOpening = 10,
         gapExtension = 4,
         scoreOnly = FALSE)
{
  ## Check arguments
  if (seqtype(pattern) != seqtype(subject))
    stop("'pattern' and 'subject' must contain ",
         "sequences of the same type")
  if (!(length(subject) %in% c(1, length(pattern))))
    stop("'length(subject)' must equal 1 or 'length(pattern)'")
  type <-
    match.arg(type,
              c("global", "local", "overlap", "global-local", "local-global",
                "subjectOverlap", "patternOverlap"))
  if (type == "subjectOverlap") {
    warning("type = 'subjectOverlap' has been renamed type = 'global-local'")
    type <- "global-local"
  }
  if (type == "patternOverlap") {
    warning("type = 'patternOverlap' has been renamed type = 'local-global'")
    type <- "local-global"
  }
  typeCode <-
    c("global" = 1L, "local" = 2L, "overlap" = 3L, "global-local" = 4L,
      "local-global" = 5L)[[type]]
  gapOpening <- as.double(abs(gapOpening))
  if (length(gapOpening) != 1 || is.na(gapOpening))
    stop("'gapOpening' must be a non-negative numeric vector of length 1")
  gapExtension <- as.double(abs(gapExtension))
  if (length(gapExtension) != 1 || is.na(gapExtension))
    stop("'gapExtension' must be a non-negative numeric vector of length 1")
  scoreOnly <- as.logical(scoreOnly)
  if (length(scoreOnly) != 1 || any(is.na(scoreOnly)))
    stop("'scoreOnly' must be a non-missing logical value")

  ## Process string information
  if (is.null(xscodec(pattern))) {
    unique_letters <- unique(c(uniqueLetters(pattern), uniqueLetters(subject)))
    #Even if safeLettersToInt() will deal properly with embedded nuls, I
    #suspect bad things would probably happen downstream in case there are any.
    alphabetToCodes <- safeLettersToInt(unique_letters, letters.as.names=TRUE)
  } else {
    alphabetToCodes <- xscodes(pattern)
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
    buildLookupTable(alphabetToCodes[availableLetters],
                     0:(length(availableLetters) - 1))

  .Call2("XStringSet_align_pairwiseAlignment",
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

.normargFuzzyMatrix <- function(fuzzyMatrix, rownames)
{
    if (is.null(fuzzyMatrix)) {
        fuzzyMatrix <- diag(length(rownames))
        dimnames(fuzzyMatrix) <- list(rownames, rownames)
        return(fuzzyMatrix)
    }
    if (!is.matrix(fuzzyMatrix) || !is.numeric(fuzzyMatrix) ||
        any(is.na(fuzzyMatrix)) || any(fuzzyMatrix < 0) ||
        any(fuzzyMatrix > 1))
        stop("'fuzzyMatrix' must be a numeric matrix with values ",
             "between 0 and 1 inclusive")
    if (!identical(rownames(fuzzyMatrix), colnames(fuzzyMatrix)))
        stop("row and column names differ for matrix 'fuzzyMatrix'")
    if (is.null(rownames(fuzzyMatrix)))
        stop("matrix 'fuzzyMatrix' must have row and column names")
    if (any(duplicated(rownames(fuzzyMatrix))))
        stop("matrix 'fuzzyMatrix' has duplicated row names")    
    availableLetters <- intersect(rownames, rownames(fuzzyMatrix))
    fuzzyMatrix[availableLetters, availableLetters, drop = FALSE]
}

.makeSubstitutionLookupTable <- function(qpattern)
{
    keys <- (minQuality(qpattern) + offset(qpattern)):
                     (maxQuality(qpattern) + offset(qpattern))
    vals <- 0:(maxQuality(qpattern) - minQuality(qpattern))
    buildLookupTable(keys, vals)
}

QualityScaledXStringSet.pairwiseAlignment <- function(pattern, subject,
                                                      type = "global",
                                                      fuzzyMatrix = NULL,
                                                      gapOpening = 10,
                                                      gapExtension = 4,
                                                      scoreOnly = FALSE)
{
    ## Check arguments
    if (class(pattern) != class(subject))
        stop("'pattern' and 'subject' must be of the same class")
    if (!(length(subject) %in% c(1L, length(pattern))))
        stop("'length(subject)' must equal 1 or 'length(pattern)'")
    type <- match.arg(type, c("global", "local", "overlap",
                              "global-local", "local-global",
                              "subjectOverlap", "patternOverlap"))
    if (type == "subjectOverlap") {
        warning("type \"subjectOverlap\" has been renamed \"global-local\"")
        type <- "global-local"
    }
    if (type == "patternOverlap") {
        warning("type \"patternOverlap\" has been renamed \"local-global\"")
        type <- "local-global"
    }
    typeCode <- c("global" = 1L, "local" = 2L, "overlap" = 3L,
                  "global-local" = 4L, "local-global" = 5L)[[type]]
    gapOpening <- as.double(abs(gapOpening))
    if (length(gapOpening) != 1L || is.na(gapOpening))
        stop("'gapOpening' must be a non-negative numeric vector of length 1")
    gapExtension <- as.double(abs(gapExtension))
    if (length(gapExtension) != 1L || is.na(gapExtension))
        stop("'gapExtension' must be a non-negative numeric vector of length 1")
    scoreOnly <- as.logical(scoreOnly)
    if (length(scoreOnly) != 1L || any(is.na(scoreOnly)))
        stop("'scoreOnly' must be a non-missing logical value")
    if (class(quality(pattern)) != class(quality(subject)))
        stop("'quality(pattern)' and 'quality(subject)' must be ",
             "of the same class")

    ## Process string information
    if (is.null(xscodec(pattern))) {
        unique_letters <- unique(c(uniqueLetters(pattern),
                                   uniqueLetters(subject)))
        #Even if safeLettersToInt() will deal properly with embedded nuls, I
        #suspect bad things will happen downstream in case there are any.
        alphabetToCodes <- safeLettersToInt(unique_letters,
                                            letters.as.names=TRUE)
    } else {
        alphabetToCodes <- xscodes(pattern)
    }

    useQuality <- TRUE
    fuzzyMatrix <- .normargFuzzyMatrix(fuzzyMatrix, names(alphabetToCodes))
    uniqueFuzzyValues <- sort(unique(fuzzyMatrix))
    fuzzyReferenceMatrix <- matrix(match(fuzzyMatrix, uniqueFuzzyValues) - 1L,
                                   nrow = nrow(fuzzyMatrix),
                                   ncol = ncol(fuzzyMatrix),
                                   dimnames = dimnames(fuzzyMatrix))
    fuzzyLookupTable <-
      buildLookupTable(alphabetToCodes[rownames(fuzzyMatrix)],
                       seq_len(nrow(fuzzyMatrix)) - 1L)
    alphabetLength <-
      switch(class(pattern),
             QualityScaledDNAStringSet =, QualityScaledRNAStringSet = 4L,
             QualityScaledAAStringSet = 20L,
             length(alphabetToCodes))
    substitutionArray <-
      qualitySubstitutionMatrices(fuzzyMatch = uniqueFuzzyValues,
                                  alphabetLength = alphabetLength,
                                  qualityClass = class(quality(pattern)))
    substitutionLookupTable <- .makeSubstitutionLookupTable(quality(pattern))

    .Call2("XStringSet_align_pairwiseAlignment",
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
  value@pattern@range <- do.call(c, lapply(mpiOutput, function(x) x@pattern@range))
  value@pattern@mismatch <- do.call(c, lapply(mpiOutput, function(x) x@pattern@mismatch))
  value@pattern@indel <- do.call(c, lapply(mpiOutput, function(x) x@pattern@indel))

  value@subject@unaligned <- subject
  value@subject@range <- do.call(c, lapply(mpiOutput, function(x) x@subject@range))
  value@subject@mismatch <- do.call(c, lapply(mpiOutput, function(x) x@subject@mismatch))
  value@subject@indel <- do.call(c, lapply(mpiOutput, function(x) x@subject@indel))

  value
}

mpi.XStringSet.pairwiseAlignment <-
  function(pattern, subject,
           type = "global",
           substitutionMatrix = NULL,
           gapOpening = 10,
           gapExtension = 4,
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
    subsets <- vector("list", k)
    for (i in seq_len(k)) {
      indices <- ((i-1)*perNode+1):ifelse(i < k, i*perNode, n)
      if (length(subject) == 1) {
          subsets[[i]] <-
            list(pattern = XStringSet(seqtype(pattern), as.character(pattern[indices])),
                 subject = subject)
      } else {
          subsets[[i]] <-
            list(pattern = XStringSet(seqtype(pattern), as.character(pattern[indices])),
                 subject = XStringSet(seqtype(subject), as.character(subject[indices])))
      }
    }

    mpi.remote.exec(library(Biostrings), ret = FALSE)
    mpiOutput <-
      mpi.parLapply(subsets,
          function(x,
                   type = "global",
                   substitutionMatrix = NULL,
                   gapOpening = 10,
                   gapExtension = 4,
                   scoreOnly = FALSE) {
            output <-
              XStringSet.pairwiseAlignment(pattern = x$pattern,
                        subject = x$subject,
                        type = type,
                        substitutionMatrix = substitutionMatrix,
                        gapOpening = gapOpening,
                        gapExtension = gapExtension,
                        scoreOnly = scoreOnly)
            if (!scoreOnly) {
              output@pattern@unaligned <- BStringSet("")
              output@subject@unaligned <- BStringSet("")
            }
            output
          },
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
  function(pattern, subject,
           type = "global",
           fuzzyMatrix = NULL,
           gapOpening = 10,
           gapExtension = 4,
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
    subsets <- vector("list", k)
    for (i in seq_len(k)) {
      indices <- ((i-1)*perNode+1):ifelse(i < k, i*perNode, n)
      if (length(subject) == 1) {
        subsets[[i]] <-
          list(pattern =
               QualityScaledXStringSet(XStringSet(seqtype(pattern),
                 as.character(pattern[indices])),
                   do.call(class(quality(pattern)),
                     list(as.character(quality(pattern[indices]))))),
               subject = subject)
      } else {
        subsets[[i]] <-
          list(pattern =
               QualityScaledXStringSet(XStringSet(seqtype(pattern),
                 as.character(pattern[indices])),
                   do.call(class(quality(pattern)),
                     list(as.character(quality(pattern[indices]))))),
               subject =
               QualityScaledXStringSet(XStringSet(seqtype(subject),
                 as.character(subject[indices])),
                   do.call(class(quality(subject)),
                     list(as.character(quality(subject[indices]))))))
      }
    }

    mpi.remote.exec(library(Biostrings), ret = FALSE)
    mpiOutput <-
      mpi.parLapply(subsets,
                    function(x,
                             type = "global",
                             fuzzyMatrix = NULL,
                             gapOpening = 10,
                             gapExtension = 4,
                             scoreOnly = FALSE) {
                      output <-
                        QualityScaledXStringSet.pairwiseAlignment(pattern = x$pattern,
                                  subject = x$subject,
                                  type = type,
                                  fuzzyMatrix = fuzzyMatrix,
                                  gapOpening = gapOpening,
                                  gapExtension = gapExtension,
                                  scoreOnly = scoreOnly)
                      if (!scoreOnly) {
                        output@pattern@unaligned <- BStringSet("")
                        output@subject@unaligned <- BStringSet("")
                      }
                      output
                    },
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### pairwiseAlignment() generic and methods.
###
### We want pairwiseAlignment() to work when 'pattern' and 'subject' are any
### of the 4 following objects: character vector, XString, XStringSet, and
### QualityScaledXStringSet. With the 4 methods defined below, we cover the 16
### type pairs we want to support (and more).
### TODO: Maybe consider making pairwiseAlignment() just an ordinary function?
###

setGeneric("pairwiseAlignment",
    function(pattern, subject, ...) standardGeneric("pairwiseAlignment")
)

setMethod("pairwiseAlignment", c("ANY", "ANY"),
    function(pattern, subject,
             patternQuality=PhredQuality(22L),
             subjectQuality=PhredQuality(22L),
             type="global",
             substitutionMatrix=NULL, fuzzyMatrix=NULL,
             gapOpening=10, gapExtension=4,
             scoreOnly=FALSE)
    {
        ## Turn each of 'pattern' and 'subject' into an instance of one of
        ## the 4 direct concrete subclasses of the XStringSet virtual class.
        pattern_seqtype <- try(seqtype(pattern), silent=TRUE)
        if (is(pattern_seqtype, "try-error"))
            pattern_seqtype <- "B"
        subject_seqtype <- try(seqtype(subject), silent=TRUE)
        if (is(subject_seqtype, "try-error"))
            subject_seqtype <- "B"
        if (pattern_seqtype == "B")
            pattern_seqtype <- subject_seqtype
        if (subject_seqtype == "B")
            subject_seqtype <- pattern_seqtype
        pattern <- XStringSet(pattern_seqtype, pattern)
        subject <- XStringSet(subject_seqtype, subject)

        if (!is.null(substitutionMatrix)) {
            mpi.XStringSet.pairwiseAlignment(pattern, subject,
                                    type=type,
                                    substitutionMatrix=substitutionMatrix,
                                    gapOpening=gapOpening,
                                    gapExtension=gapExtension,
                                    scoreOnly=scoreOnly)
        } else {
            pattern <- QualityScaledXStringSet(pattern, patternQuality)
            subject <- QualityScaledXStringSet(subject, subjectQuality)
            mpi.QualityScaledXStringSet.pairwiseAlignment(pattern, subject,
                                    type=type,
                                    fuzzyMatrix=fuzzyMatrix,
                                    gapOpening=gapOpening,
                                    gapExtension=gapExtension,
                                    scoreOnly=scoreOnly)
        }
    }
)

setMethod("pairwiseAlignment", c("ANY", "QualityScaledXStringSet"),
    function(pattern, subject,
             patternQuality=PhredQuality(22L),
             type="global",
             substitutionMatrix=NULL, fuzzyMatrix=NULL,
             gapOpening=10, gapExtension=4,
             scoreOnly=FALSE)
    {
        if (is.character(pattern)) {
            pattern <- XStringSet(seqtype(subject), pattern)
        } else {
            pattern <- as(pattern, "XStringSet")
        }
        if (!is.null(substitutionMatrix)) {
            subject <- as(subject, "XStringSet")
            mpi.XStringSet.pairwiseAlignment(pattern, subject,
                                    type=type,
                                    substitutionMatrix=substitutionMatrix,
                                    gapOpening=gapOpening,
                                    gapExtension=gapExtension,
                                    scoreOnly=scoreOnly)
        } else {
            pattern <- QualityScaledXStringSet(pattern, patternQuality)
            mpi.QualityScaledXStringSet.pairwiseAlignment(pattern, subject,
                                    type=type,
                                    fuzzyMatrix=fuzzyMatrix,
                                    gapOpening=gapOpening,
                                    gapExtension=gapExtension,
                                    scoreOnly=scoreOnly)
        }
    }
)

setMethod("pairwiseAlignment", c("QualityScaledXStringSet", "ANY"),
    function(pattern, subject,
             subjectQuality=PhredQuality(22L),
             type="global",
             substitutionMatrix=NULL, fuzzyMatrix=NULL,
             gapOpening=10, gapExtension=4,
             scoreOnly=FALSE)
    {
        if (is.character(subject)) {
            subject <- XStringSet(seqtype(pattern), subject)
        } else {
            subject <- as(subject, "XStringSet")
        }
        if (!is.null(substitutionMatrix)) {
            pattern <- as(pattern, "XStringSet")
            mpi.XStringSet.pairwiseAlignment(pattern, subject,
                                    type=type,
                                    substitutionMatrix=substitutionMatrix,
                                    gapOpening=gapOpening,
                                    gapExtension=gapExtension,
                                    scoreOnly=scoreOnly)
        } else {
            subject <- QualityScaledXStringSet(subject, subjectQuality)
            mpi.QualityScaledXStringSet.pairwiseAlignment(pattern, subject,
                                    type=type,
                                    fuzzyMatrix=fuzzyMatrix,
                                    gapOpening=gapOpening,
                                    gapExtension=gapExtension,
                                    scoreOnly=scoreOnly)
        }
    }
)

setMethod("pairwiseAlignment", c("QualityScaledXStringSet",
                                 "QualityScaledXStringSet"),
    function(pattern, subject,
             type="global",
             substitutionMatrix=NULL, fuzzyMatrix=NULL,
             gapOpening=10, gapExtension=4,
             scoreOnly=FALSE)
    {
        if (!is.null(substitutionMatrix)) {
            pattern <- as(pattern, "XStringSet")
            subject <- as(subject, "XStringSet")
            mpi.XStringSet.pairwiseAlignment(pattern, subject,
                                    type=type,
                                    substitutionMatrix=substitutionMatrix,
                                    gapOpening=gapOpening,
                                    gapExtension=gapExtension,
                                    scoreOnly=scoreOnly)
        } else {
            mpi.QualityScaledXStringSet.pairwiseAlignment(pattern, subject,
                                    type=type,
                                    fuzzyMatrix=fuzzyMatrix,
                                    gapOpening=gapOpening,
                                    gapExtension=gapExtension,
                                    scoreOnly=scoreOnly)
        }
    }
)

