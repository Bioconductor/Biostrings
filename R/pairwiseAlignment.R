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
  if (!is(new(qualityClass), "XStringQuality"))
    stop("'qualityClass' must be one of the 'XStringQuality' classes")
  qualityIntegers <- minQuality(new(qualityClass)):maxQuality(new(qualityClass))
  prob <- qualityConverter(qualityIntegers, qualityClass, "numeric")
  names(prob) <- as.character(qualityIntegers)
  errorSubstitutionMatrices(prob, alphabetLength = alphabetLength, bitScale = bitScale)
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
  #Herve - Sept 30, 2008 - The test below is probably too stringent. It will
  #reject a pattern or a subject that both extend DNAStringSet for example.
  #Testing whether baseXStringSubtype(pattern) != baseXStringSubtype(subject)
  #should be enough.
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

  ## Process string information
  if (is.null(codec(pattern))) {
    #Herve - Sept 30, 2008 - Never assume that the elements of an XStringSet
    #object 'x' cover 'super(x)' entirely! It's not true in general. For example
    #when 'x' is the result of the coercion of an XStringViews object where the
    #views don't cover the subject entirely. Or when 'x' has is the result of a
    #narrowing operation with DNAStringSet(x0, start=4, end=-4).
    #It's always better to try to remain agnostic about the internals of an
    #XStringSet object e.g. you should avoid to use 'super(x)' (the internals
    #of the XStringSet container will change and the slot @super will go away).
    #For uniqueBases you could do:
    #  afp <- alphabetFrequency(pattern, collapse=TRUE)
    #  afs <- alphabetFrequency(subject, collapse=TRUE)
    #  unique_letters <- unique(c(names(afp)[afp != 0], names(afs)[afs != 0]))
    #  alphabetToCodes <- as.integer(charToRaw(paste(unique_letters, collapse="")))
    #  names(alphabetToCodes) <- unique_letters
    #It does the right thing, whether 'super(pattern)' is entirely covered or not.
    #It's up to 20x faster and much more memory efficient (calls to as.character()
    #should be avoided at all cost, there are almost always better ways).
    #It won't break when the internals of the XStringSet container change.
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

  useQuality <- FALSE
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

  .Call("XStringSet_align_pairwiseAlignment",
        pattern,
        subject,
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
}


QualityScaledXStringSet.pairwiseAlignment <-
function(pattern,
         subject,
         type = "global",
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

  useQuality <- TRUE
  alphabetLength <-
    switch(class(pattern),
           QualityScaledDNAStringSet =, QualityScaledRNAStringSet = 4L,
           QualityScaledAAStringSet = 20L,
           256L)

  qualityLookupTable <-
    buildLookupTable((minQuality(quality(pattern)) + offset(quality(pattern))):
                     (maxQuality(quality(pattern)) + offset(quality(pattern))),
                     0:(maxQuality(quality(pattern)) - minQuality(quality(pattern))))
  qualityMatrices <-
    qualitySubstitutionMatrices(alphabetLength = alphabetLength,
                                qualityClass = class(quality(pattern)))

  constantLookupTable <- integer(0)
  constantMatrix <- matrix(numeric(0), nrow = 0, ncol = 0)

  .Call("XStringSet_align_pairwiseAlignment",
        pattern,
        subject,
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
                   scoreOnly = FALSE) {
            if (!is.null(substitutionMatrix)) {
              XStringSet.pairwiseAlignment(pattern = BStringSet(pattern),
                                           subject = BStringSet(subject),
                                           type = type,
                                           substitutionMatrix = substitutionMatrix,
                                           gapExtension = gapExtension,
                                           gapOpening = gapOpening,
                                           scoreOnly = scoreOnly)
            } else {
              QualityScaledXStringSet.pairwiseAlignment(
                  pattern = QualityScaledBStringSet(BStringSet(pattern), patternQuality),
                  subject = QualityScaledBStringSet(BStringSet(subject), subjectQuality),
                  type = type,
                  gapExtension = gapExtension,
                  gapOpening = gapOpening,
                  scoreOnly = scoreOnly)
          }})

setMethod("pairwiseAlignment",
          signature(pattern = "character", subject = "XString"),
          function(pattern, subject, patternQuality = PhredQuality(22L),
                   subjectQuality = PhredQuality(22L), type = "global",
                   substitutionMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE) {
            if (!is.null(substitutionMatrix)) {
              XStringSet.pairwiseAlignment(pattern = XStringSet(class(subject), pattern),
                                           subject = XStringSet(class(subject), subject),
                                           type = type,
                                           substitutionMatrix = substitutionMatrix,
                                           gapExtension = gapExtension,
                                           gapOpening = gapOpening,
                                           scoreOnly = scoreOnly)
            } else {
              QualityScaledXStringSet.pairwiseAlignment(
                  pattern = QualityScaledXStringSet(XStringSet(class(subject), pattern), patternQuality),
                  subject = QualityScaledXStringSet(XStringSet(class(subject), subject), subjectQuality),
                  type = type,
                  gapExtension = gapExtension,
                  gapOpening = gapOpening,
                  scoreOnly = scoreOnly)
            }})

setMethod("pairwiseAlignment",
          signature(pattern = "character", subject = "XStringSet"),
          function(pattern, subject, patternQuality = PhredQuality(22L),
                   subjectQuality = PhredQuality(22L), type = "global",
                   substitutionMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE) {
            if (!is.null(substitutionMatrix)) {
              XStringSet.pairwiseAlignment(pattern = XStringSet(baseXStringSubtype(subject), pattern),
                                           subject = subject,
                                           type = type,
                                           substitutionMatrix = substitutionMatrix,
                                           gapExtension = gapExtension,
                                           gapOpening = gapOpening,
                                           scoreOnly = scoreOnly)
            } else {
              QualityScaledXStringSet.pairwiseAlignment(
                  pattern = QualityScaledXStringSet(XStringSet(baseXStringSubtype(subject), pattern), patternQuality),
                  subject = QualityScaledXStringSet(subject, subjectQuality),
                  type = type,
                  gapExtension = gapExtension,
                  gapOpening = gapOpening,
                  scoreOnly = scoreOnly)
          }})

setMethod("pairwiseAlignment",
          signature(pattern = "character", subject = "QualityScaledXStringSet"),
          function(pattern, subject, patternQuality = PhredQuality(22L), type = "global",
                   substitutionMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE) {
            if (!is.null(substitutionMatrix)) {
              XStringSet.pairwiseAlignment(pattern = XStringSet(baseXStringSubtype(subject), pattern),
                                           subject = as(subject, "XStringSet"),
                                           type = type,
                                           substitutionMatrix = substitutionMatrix,
                                           gapExtension = gapExtension,
                                           gapOpening = gapOpening,
                                           scoreOnly = scoreOnly)
            } else {
              QualityScaledXStringSet.pairwiseAlignment(
                  pattern = QualityScaledXStringSet(XStringSet(baseXStringSubtype(subject), pattern), patternQuality),
                  subject = subject,
                  type = type,
                  gapExtension = gapExtension,
                  gapOpening = gapOpening,
                  scoreOnly = scoreOnly)
            }})

setMethod("pairwiseAlignment",
          signature(pattern = "XString", subject = "character"),
          function(pattern, subject, patternQuality = PhredQuality(22L),
                   subjectQuality = PhredQuality(22L), type = "global",
                   substitutionMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE) {
            if (!is.null(substitutionMatrix)) {
              XStringSet.pairwiseAlignment(pattern = XStringSet(class(pattern), pattern),
                                           subject = XStringSet(class(pattern), subject),
                                           type = type,
                                           substitutionMatrix = substitutionMatrix,
                                           gapExtension = gapExtension,
                                           gapOpening = gapOpening,
                                           scoreOnly = scoreOnly)
            } else {
              QualityScaledXStringSet.pairwiseAlignment(
                  pattern = QualityScaledXStringSet(XStringSet(class(pattern), pattern), patternQuality),
                  subject = QualityScaledXStringSet(XStringSet(class(pattern), subject), subjectQuality),
                  type = type,
                  gapExtension = gapExtension,
                  gapOpening = gapOpening,
                  scoreOnly = scoreOnly)
          }})

setMethod("pairwiseAlignment",
          signature(pattern = "XString", subject = "XString"),
          function(pattern, subject, patternQuality = PhredQuality(22L),
                   subjectQuality = PhredQuality(22L), type = "global",
                   substitutionMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE) {
            if (!is.null(substitutionMatrix)) {
              XStringSet.pairwiseAlignment(pattern = XStringSet(class(pattern), pattern),
                                           subject = XStringSet(class(subject), subject),
                                           type = type,
                                           substitutionMatrix = substitutionMatrix,
                                           gapExtension = gapExtension,
                                           gapOpening = gapOpening,
                                           scoreOnly = scoreOnly)
            } else {
              QualityScaledXStringSet.pairwiseAlignment(
                  pattern = QualityScaledXStringSet(XStringSet(class(pattern), pattern), patternQuality),
                  subject = QualityScaledXStringSet(XStringSet(class(subject), subject), subjectQuality),
                  type = type,
                  gapExtension = gapExtension,
                  gapOpening = gapOpening,
                  scoreOnly = scoreOnly)
            }})

setMethod("pairwiseAlignment",
          signature(pattern = "XString", subject = "XStringSet"),
          function(pattern, subject, patternQuality = PhredQuality(22L),
                   subjectQuality = PhredQuality(22L), type = "global",
                   substitutionMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE) {
            if (!is.null(substitutionMatrix)) {
              XStringSet.pairwiseAlignment(pattern = XStringSet(class(pattern), pattern),
                                           subject = subject,
                                           type = type,
                                           substitutionMatrix = substitutionMatrix,
                                           gapExtension = gapExtension,
                                           gapOpening = gapOpening,
                                           scoreOnly = scoreOnly)
            } else {
              QualityScaledXStringSet.pairwiseAlignment(
                  pattern = QualityScaledXStringSet(XStringSet(class(pattern), pattern), patternQuality),
                  subject = QualityScaledXStringSet(subject, subjectQuality),
                  type = type,
                  gapExtension = gapExtension,
                  gapOpening = gapOpening,
                  scoreOnly = scoreOnly)
            }})

setMethod("pairwiseAlignment",
          signature(pattern = "XString", subject = "QualityScaledXStringSet"),
          function(pattern, subject, patternQuality = PhredQuality(22L), type = "global",
                   substitutionMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE) {
            if (!is.null(substitutionMatrix)) {
              XStringSet.pairwiseAlignment(pattern = XStringSet(class(pattern), pattern),
                                           subject = as(subject, "XStringSet"),
                                           type = type,
                                           substitutionMatrix = substitutionMatrix,
                                           gapExtension = gapExtension,
                                           gapOpening = gapOpening,
                                           scoreOnly = scoreOnly)
            } else {
              QualityScaledXStringSet.pairwiseAlignment(
                  pattern = QualityScaledXStringSet(XStringSet(class(pattern), pattern), patternQuality),
                  subject = subject,
                  type = type,
                  gapExtension = gapExtension,
                  gapOpening = gapOpening,
                  scoreOnly = scoreOnly)
            }})

setMethod("pairwiseAlignment",
          signature(pattern = "XStringSet", subject = "character"),
          function(pattern, subject, patternQuality = PhredQuality(22L),
                   subjectQuality = PhredQuality(22L), type = "global",
                   substitutionMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE) {
          if (!is.null(substitutionMatrix)) {
            XStringSet.pairwiseAlignment(pattern = pattern,
                                         subject = XStringSet(baseXStringSubtype(pattern), subject),
                                         type = type,
                                         substitutionMatrix = substitutionMatrix,
                                         gapExtension = gapExtension,
                                         gapOpening = gapOpening,
                                         scoreOnly = scoreOnly)
          } else {
            QualityScaledXStringSet.pairwiseAlignment(
                pattern = QualityScaledXStringSet(pattern, patternQuality),
                subject = QualityScaledXStringSet(XStringSet(baseXStringSubtype(pattern), subject), subjectQuality),
                type = type,
                gapExtension = gapExtension,
                gapOpening = gapOpening,
                scoreOnly = scoreOnly)
          }})

setMethod("pairwiseAlignment",
          signature(pattern = "XStringSet", subject = "XString"),
          function(pattern, subject, patternQuality = PhredQuality(22L),
                   subjectQuality = PhredQuality(22L), type = "global",
                   substitutionMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE) {
            if (!is.null(substitutionMatrix)) {
              XStringSet.pairwiseAlignment(pattern = pattern,
                                           subject = XStringSet(class(subject), subject),
                                           type = type,
                                           substitutionMatrix = substitutionMatrix,
                                           gapExtension = gapExtension,
                                           gapOpening = gapOpening,
                                           scoreOnly = scoreOnly)
            } else {
              QualityScaledXStringSet.pairwiseAlignment(
                  pattern = QualityScaledXStringSet(pattern, patternQuality),
                  subject = QualityScaledXStringSet(XStringSet(class(subject), subject), subjectQuality),
                  type = type,
                  gapExtension = gapExtension,
                  gapOpening = gapOpening,
                  scoreOnly = scoreOnly)
            }})

setMethod("pairwiseAlignment",
          signature(pattern = "XStringSet", subject = "XStringSet"),
          function(pattern, subject, patternQuality = PhredQuality(22L),
                   subjectQuality = PhredQuality(22L), type = "global",
                   substitutionMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE) {
            if (!is.null(substitutionMatrix)) {
              XStringSet.pairwiseAlignment(pattern = pattern,
                                           subject = subject,
                                           type = type,
                                           substitutionMatrix = substitutionMatrix,
                                           gapExtension = gapExtension,
                                           gapOpening = gapOpening,
                                           scoreOnly = scoreOnly)
            } else {
              QualityScaledXStringSet.pairwiseAlignment(
                  pattern = QualityScaledXStringSet(pattern, patternQuality),
                  subject = QualityScaledXStringSet(subject, subjectQuality),
                  type = type,
                  gapExtension = gapExtension,
                  gapOpening = gapOpening,
                  scoreOnly = scoreOnly)
            }})

setMethod("pairwiseAlignment",
          signature(pattern = "XStringSet", subject = "QualityScaledXStringSet"),
          function(pattern, subject, patternQuality = PhredQuality(22L), type = "global",
                   substitutionMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE) {
            if (!is.null(substitutionMatrix)) {
              XStringSet.pairwiseAlignment(pattern = pattern,
                                           subject = as(subject, "XStringSet"),
                                           type = type,
                                           substitutionMatrix = substitutionMatrix,
                                           gapExtension = gapExtension,
                                           gapOpening = gapOpening,
                                           scoreOnly = scoreOnly)
            } else {
              QualityScaledXStringSet.pairwiseAlignment(
                  pattern = QualityScaledXStringSet(pattern, patternQuality),
                  subject = subject,
                  type = type,
                  gapExtension = gapExtension,
                  gapOpening = gapOpening,
                  scoreOnly = scoreOnly)
            }})

setMethod("pairwiseAlignment",
          signature(pattern = "QualityScaledXStringSet", subject = "character"),
          function(pattern, subject, subjectQuality = PhredQuality(22L), type = "global",
                   substitutionMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE) {
            if (!is.null(substitutionMatrix)) {
              XStringSet.pairwiseAlignment(pattern = as(pattern, "XStringSet"),
                                           subject = XStringSet(baseXStringSubtype(pattern), subject),
                                           type = type,
                                           substitutionMatrix = substitutionMatrix,
                                           gapExtension = gapExtension,
                                           gapOpening = gapOpening,
                                           scoreOnly = scoreOnly)
            } else {
              QualityScaledXStringSet.pairwiseAlignment(
                  pattern = pattern,
                  subject = QualityScaledXStringSet(XStringSet(baseXStringSubtype(pattern), subject), subjectQuality),
                  type = type,
                  gapExtension = gapExtension,
                  gapOpening = gapOpening,
                  scoreOnly = scoreOnly)
            }})

setMethod("pairwiseAlignment",
          signature(pattern = "QualityScaledXStringSet", subject = "XString"),
          function(pattern, subject, subjectQuality = PhredQuality(22L), type = "global",
                   substitutionMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE) {
            if (!is.null(substitutionMatrix)) {
              XStringSet.pairwiseAlignment(pattern = as(pattern, "XStringSet"),
                                           subject = XStringSet(class(subject), subject),
                                           type = type,
                                           substitutionMatrix = substitutionMatrix,
                                           gapExtension = gapExtension,
                                           gapOpening = gapOpening,
                                           scoreOnly = scoreOnly)
            } else {
              QualityScaledXStringSet.pairwiseAlignment(
                  pattern = pattern,
                  subject = QualityScaledXStringSet(XStringSet(class(subject), subject), subjectQuality),
                  type = type,
                  gapExtension = gapExtension,
                  gapOpening = gapOpening,
                  scoreOnly = scoreOnly)
            }})

setMethod("pairwiseAlignment",
          signature(pattern = "QualityScaledXStringSet", subject = "XStringSet"),
          function(pattern, subject, subjectQuality = PhredQuality(22L), type = "global",
                   substitutionMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE) {
            if (!is.null(substitutionMatrix)) {
              XStringSet.pairwiseAlignment(pattern = as(pattern, "XStringSet"),
                                           subject = subject,
                                           type = type,
                                           substitutionMatrix = substitutionMatrix,
                                           gapExtension = gapExtension,
                                           gapOpening = gapOpening,
                                           scoreOnly = scoreOnly)
            } else {
              QualityScaledXStringSet.pairwiseAlignment(
                  pattern = pattern,
                  subject = QualityScaledXStringSet(subject, subjectQuality),
                  type = type,
                  gapExtension = gapExtension,
                  gapOpening = gapOpening,
                  scoreOnly = scoreOnly)
            }})

setMethod("pairwiseAlignment",
          signature(pattern = "QualityScaledXStringSet", subject = "QualityScaledXStringSet"),
          function(pattern, subject, type = "global",
                   substitutionMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE) {
            if (!is.null(substitutionMatrix)) {
              XStringSet.pairwiseAlignment(pattern = as(pattern, "XStringSet"),
                                           subject = as(subject, "XStringSet"),
                                           type = type,
                                           substitutionMatrix = substitutionMatrix,
                                           gapExtension = gapExtension,
                                           gapOpening = gapOpening,
                                           scoreOnly = scoreOnly)
            } else {
              QualityScaledXStringSet.pairwiseAlignment(
                  pattern = pattern,
                  subject = subject,
                  type = type,
                  gapExtension = gapExtension,
                  gapOpening = gapOpening,
                  scoreOnly = scoreOnly)
            }})
