### =========================================================================
### PairwiseAlignments objects
### -------------------------------------------------------------------------
###
### A PairwiseAlignments object contains two aligned XStringSet objects.
###

setClass("PairwiseAlignments",
    contains="Vector",
    representation(
        pattern="AlignedXStringSet0",  # of length N
        subject="AlignedXStringSet0",  # of length N
        score="numeric",               # of length N
        type="character",
        gapOpening="numeric",
        gapExtension="numeric"
    )
)

### Combine the new parallel slots with those of the parent class. Make sure
### to put the new parallel slots *first*.
setMethod("parallelSlotNames", "PairwiseAlignments",
    function(x) c("pattern", "subject", "score", callNextMethod())
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructors.
###

newPairwiseAlignments <-
function(pattern, subject, type = "global", substitutionMatrix = NULL,
         gapOpening = 0, gapExtension = 1,
         baseClass = "BString", pwaClass = "PairwiseAlignments") {
    seqtype <- substr(baseClass, 1, nchar(baseClass) - 6)  # remove "String" suffix
    getMismatches <- function(x) {
        whichMismatches <- which(x[["values"]] == "?")
        if (length(whichMismatches) == 0) {
            value <- integer(0)
        } else {
            start <- cumsum(x[["lengths"]])[whichMismatches]
            end <- start + (x[["lengths"]][whichMismatches] - 1L)
            value <-
              eval(parse(text =
                   paste("c(", paste(start, ":", end, sep = "", collapse = ", "), ")")))
        }
        IntegerList(value)
    }
    getRange <- function(x) {
        if (!(x[["values"]][1] %in% c("-", "+"))) {
            start <- 1L
        } else if (length(x[["values"]]) == 1) {
            start <- integer(0)
        } else {
            start <- x[["lengths"]][1] + 1L
        }
        if (!(x[["values"]][length(x[["values"]])] %in% c("-", "+"))) {
            end <- sum(x[["lengths"]])
        } else if (length(x[["values"]]) == 1) {
            end <- integer(0)
        } else {
            end <- sum(x[["lengths"]][-length(x[["lengths"]])])
        }
        IRanges(start = start, end = end)
    }
    getIndels <- function(x, indelChar) {
        if (x[["values"]][1] %in% c("-", "+")) {
            x[["values"]] <- x[["values"]][-1]
            x[["lengths"]] <- x[["lengths"]][-1]
        }
        if (x[["values"]][length(x[["values"]])] %in% c("-", "+")) {
            x[["values"]] <- x[["values"]][-length(x[["values"]])]
            x[["lengths"]] <- x[["lengths"]][-length(x[["lengths"]])]
        }
        isIndels <- (x[["values"]] == indelChar)
        if (!any(isIndels))
            IRangesList(IRanges(integer(0), integer(0)))
        else
            IRangesList(IRanges(
                            cumsum(c(1L, ifelse(isIndels, 0L, x[["lengths"]])[-length(x[["lengths"]])]))[isIndels],
                            width = x[["lengths"]][isIndels]))
    }
    if (length(pattern) != 1 || length(subject) != 1)
        stop("'pattern' and 'subject' must both be of length 1")
    if (nchar(pattern) != nchar(subject))
        stop("'pattern' and 'subject' must have the same number of characters")
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
    gapOpening <- as.double(abs(gapOpening))
    if (length(gapOpening) != 1 || is.na(gapOpening))
        stop("'gapOpening' must be a non-negative numeric vector of length 1")
    gapExtension <- as.double(abs(gapExtension))
    if (length(gapExtension) != 1 || is.na(gapExtension))
        stop("'gapExtension' must be a non-negative numeric vector of length 1")

    explodedPattern <- safeExplode(pattern)
    explodedSubject <- safeExplode(subject)
    degappedPattern <- explodedPattern[explodedPattern != "-"]
    degappedSubject <- explodedSubject[explodedSubject != "-"]
    availableLetters <-
      sort(unique(c(unique(degappedPattern), unique(degappedSubject))))
    if (is.null(substitutionMatrix)) {
        substitutionMatrix <- diag(length(availableLetters)) - 1
        dimnames(substitutionMatrix) <- list(availableLetters, availableLetters)
    } else if (is.character(substitutionMatrix)) {
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
      intersect(availableLetters, rownames(substitutionMatrix))

    substitutionMatrix <-
      matrix(as.double(substitutionMatrix[availableLetters, availableLetters]),
             nrow = length(availableLetters),
             ncol = length(availableLetters),
             dimnames = list(availableLetters, availableLetters))

    comparison <- rle(safeExplode(compareStrings(pattern, subject)))
    whichPattern <- which(comparison[["values"]] != "-")
    patternRle <-
      structure(list(lengths = comparison[["lengths"]][whichPattern],
                     values = comparison[["values"]][whichPattern]),
                class = "rle")
    whichSubject <- which(comparison[["values"]] != "+")
    subjectRle <-
      structure(list(lengths = comparison[["lengths"]][whichSubject],
                     values = comparison[["values"]][whichSubject]),
                class = "rle")
    substitutionIndices <- (explodedPattern != "-") & (explodedSubject != "-")
    new(pwaClass,
        pattern =
          new("AlignedXStringSet",
              unaligned = XStringSet(seqtype, paste(degappedPattern, collapse = "")),
              range = getRange(patternRle), mismatch = getMismatches(patternRle),
              indel = getIndels(comparison, "-")),
        subject =
          new("AlignedXStringSet",
              unaligned = XStringSet(seqtype, paste(degappedSubject, collapse = "")),
              range = getRange(subjectRle), mismatch = getMismatches(subjectRle),
              indel = getIndels(comparison, "+")),
        type = type,
        score =
          sum(substitutionMatrix[
                match(explodedPattern[substitutionIndices], availableLetters) +
                  length(availableLetters) *
                     (match(explodedSubject[substitutionIndices], availableLetters) - 1)]) +
            gapOpening * sum(comparison[["values"]] %in% c("+", "-")) +
              gapExtension * sum(comparison[["lengths"]][comparison[["values"]] %in% c("+", "-")]),
        gapOpening = gapOpening, gapExtension = gapExtension)
}

setGeneric("PairwiseAlignments",
    function(pattern, subject, type = "global", substitutionMatrix = NULL,
             gapOpening = 0, gapExtension = 1, ...)
        standardGeneric("PairwiseAlignments"))

setMethod("PairwiseAlignments", signature(pattern = "XString", subject = "XString"),
    function(pattern, subject, type = "global", substitutionMatrix = NULL,
             gapOpening = 0, gapExtension = 1) {
        if (seqtype(pattern) != seqtype(subject))
            stop("'pattern' and 'subject' must contain ",
                 "sequences of the same type")
        PairwiseAlignments(as.character(pattern), as.character(subject),
                                  type = type, substitutionMatrix = substitutionMatrix,
                                  gapOpening = gapOpening, gapExtension = gapExtension,
                                  baseClass = xsbaseclass(pattern))
    }
)

setMethod("PairwiseAlignments", signature(pattern = "XStringSet", subject = "missing"),
    function(pattern, subject, type = "global", substitutionMatrix = NULL,
             gapOpening = 0, gapExtension = 1) {
        if (length(pattern) != 2)
            stop("'pattern' must be of length 2 when 'subject' is missing")
        if (diff(nchar(pattern)) != 0)
            stop("'pattern' elements must have the same number of characters")
        PairwiseAlignments(as.character(pattern[1]), as.character(pattern[2]),
                                  type = type, substitutionMatrix = substitutionMatrix,
                                  gapOpening = gapOpening, gapExtension = gapExtension,
                                  baseClass = xsbaseclass(pattern))
    }
)

setMethod("PairwiseAlignments", signature(pattern = "character", subject = "missing"),
    function(pattern, subject, type = "global", substitutionMatrix = NULL,
             gapOpening = 0, gapExtension = 1, baseClass = "BString") {
        if (length(pattern) != 2)
            stop("'pattern' must be of length 2 when 'subject' is missing")
        if (diff(nchar(pattern)) != 0)
            stop("'pattern' elements must have the same number of characters")
        PairwiseAlignments(pattern[1], pattern[2],
                                  type = type, substitutionMatrix = substitutionMatrix,
                                  gapOpening = gapOpening, gapExtension = gapExtension,
                                  baseClass = baseClass)
    }
)

setMethod("PairwiseAlignments", signature(pattern = "character", subject = "character"),
    function(pattern, subject, type = "global", substitutionMatrix = NULL,
             gapOpening = 0, gapExtension = 1, baseClass = "BString") {
        newPairwiseAlignments(pattern = pattern, subject = subject, type = type,
                                     substitutionMatrix = substitutionMatrix,
                                     gapOpening = gapOpening, gapExtension = gapExtension,
                                     baseClass = baseClass, pwaClass = "PairwiseAlignments")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.PairwiseAlignments <- function(object)
{
    message <- NULL
    if (!identical(class(unaligned(pattern(object))), class(unaligned(subject(object)))))
        message <- c(message, "'unaligned(pattern)' and 'unaligned(subject)' must be XString objects of the same base type")
    if (length(object@type) != 1 || !(object@type %in% c("global", "local", "overlap", "global-local", "local-global")))
        message <- c(message, "'type' must be one of 'global', 'local', 'overlap', 'global-local', or 'local-global'")
    if (!isSingleNumber(object@gapOpening) || object@gapOpening < 0)
        message <- c(message, "'gapOpening' must be a non-negative numeric vector of length 1")
    if (!isSingleNumber(object@gapExtension) || object@gapExtension < 0)
        message <- c(message, "'gapExtension' must be a non-negative numeric vector of length 1")
    message
}

setValidity("PairwiseAlignments",
    function(object)
    {
       problems <- .valid.PairwiseAlignments(object)
        if (is.null(problems)) TRUE else problems
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods.
###

setMethod("pattern", "PairwiseAlignments", function(x) x@pattern)
setMethod("subject", "PairwiseAlignments", function(x) x@subject)
setGeneric("type", function(x) standardGeneric("type"))
setMethod("type", "PairwiseAlignments", function(x) x@type)
setMethod("score", "PairwiseAlignments", function(x) x@score)
setMethod("insertion", "PairwiseAlignments", function(x) indel(subject(x)))
setMethod("deletion", "PairwiseAlignments", function(x) indel(pattern(x)))
setMethod("indel", "PairwiseAlignments",
          function(x) new("InDel", insertion = insertion(x), deletion = deletion(x)))
setMethod("nindel", "PairwiseAlignments",
          function(x) new("InDel", insertion = nindel(subject(x)), deletion = nindel(pattern(x))))
setMethod("nchar", "PairwiseAlignments", function(x, type="chars", allowNA=FALSE) nchar(subject(x)))
setMethod("seqtype", "PairwiseAlignments", function(x) seqtype(subject(x)))
setGeneric("pid", signature="x", function(x, type="PID1") standardGeneric("pid"))
setMethod("pid", "PairwiseAlignments",
          function(x, type="PID1") {
              type <- match.arg(type, c("PID1", "PID2", "PID3", "PID4"))
              denom <-
                switch(type,
                       "PID1" = nchar(x),
                       "PID2" = nmatch(x) + nmismatch(x),
                       "PID3" = pmin(nchar(unaligned(pattern(x))), nchar(unaligned(subject(x)))),
                       "PID4" = (nchar(unaligned(pattern(x))) + nchar(unaligned(subject(x)))) / 2)
              100 * nmatch(x)/denom
		  })


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method.
###

### TODO: Make the "show" method to format the alignment in a SGD fashion
### i.e. split in 60-letter blocks and use the "|" character to highlight
### exact matches.
setMethod("show", "PairwiseAlignments",
    function(object) {
        if (length(object) == 0)
            cat("Empty Pairwise Alignment\n")
        else {
            cat(switch(type(object), "global" = "Global", "overlap" = "Overlap",
                       "local" = "Local", "global-local" = "Global-Local",
                       "local-global" = "Local-Global"),
                " ", class(object), " (1 of ", length(object), ")\n", sep = "")
            if (width(pattern(object))[1] == 0 || width(subject(object))[1] == 0) {
                patternSpaces <- 0
                subjectSpaces <- 0
            } else {
                patternSpaces <-
                  floor(log10(start(subject(object))[1])) - floor(log10(start(pattern(object))[1]))
		        subjectSpaces <- max(0, - patternSpaces)
		        patternSpaces <- max(0, patternSpaces)
            }
            cat(paste(c("pattern: ", rep(" ", patternSpaces)), collapse = ""))
            show(pattern(object)[1])
            cat(paste(c("subject: ", rep(" ", subjectSpaces)), collapse = ""))
            show(subject(object)[1])
            cat("score:", score(object)[1], "\n")
        }
    }
)

