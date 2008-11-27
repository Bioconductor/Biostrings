### ==========================================================================
### PairwiseAlignment objects
### --------------------------------------------------------------------------
### An PairwiseAlignment object contains the result of the alignment of
### two XString objects of the same subtype.


setClass("PairwiseAlignment",
    representation(
        pattern="AlignedXStringSet",
        subject="AlignedXStringSet",
        type="character",
        score="numeric",
        substitutionArray="array",
        gapOpening="numeric",
        gapExtension="numeric"
    )
)

setClass("PairwiseAlignmentSummary",
    representation(
        type="character",
        score="numeric",
        nmatch="numeric",
        nmismatch="numeric",
        ninsertion="matrix",
        ndeletion="matrix",
        coverage="XRleInteger",
        mismatchSummary="list"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructors.
###

setGeneric("PairwiseAlignment",
    function(pattern, subject, type = "global", substitutionMatrix = NULL,
             gapOpening = 0, gapExtension = -1, ...)
    standardGeneric("PairwiseAlignment"))
setMethod("PairwiseAlignment", signature(pattern = "XString", subject = "XString"),
    function(pattern, subject, type = "global", substitutionMatrix = NULL,
             gapOpening = 0, gapExtension = -1) {
        if (baseXStringSubtype(pattern) != baseXStringSubtype(subject))
            stop("'pattern' and 'subject' must have the same XString base subtype")
        PairwiseAlignment(as.character(pattern), as.character(subject),
                          type = type, substitutionMatrix = substitutionMatrix,
                          gapOpening = gapOpening, gapExtension = gapExtension,
                          baseClass = baseXStringSubtype(pattern))
    }
)

setMethod("PairwiseAlignment", signature(pattern = "XStringSet", subject = "missing"),
    function(pattern, subject, type = "global", substitutionMatrix = NULL,
             gapOpening = 0, gapExtension = -1) {
        if (length(pattern) != 2)
            stop("'pattern' must be of length 2 when 'subject' is missing")
        if (diff(nchar(pattern)) != 0)
            stop("'pattern' elements must have the same number of characters")
        PairwiseAlignment(as.character(pattern[1]), as.character(pattern[2]),
                          type = type, substitutionMatrix = substitutionMatrix,
                          gapOpening = gapOpening, gapExtension = gapExtension,
                          baseClass = baseXStringSubtype(pattern))
    }
)

setMethod("PairwiseAlignment", signature(pattern = "character", subject = "missing"),
    function(pattern, subject, type = "global", substitutionMatrix = NULL,
             gapOpening = 0, gapExtension = -1, baseClass = "BString") {
        if (length(pattern) != 2)
            stop("'pattern' must be of length 2 when 'subject' is missing")
        if (diff(nchar(pattern)) != 0)
            stop("'pattern' elements must have the same number of characters")
        PairwiseAlignment(pattern[1], pattern[2],
                          type = type, substitutionMatrix = substitutionMatrix,
                          gapOpening = gapOpening, gapExtension = gapExtension,
                          baseClass = baseClass)
    }
)

setMethod("PairwiseAlignment", signature(pattern = "character", subject = "character"),
    function(pattern, subject, type = "global", substitutionMatrix = NULL,
             gapOpening = 0, gapExtension = -1, baseClass = "BString") {
        getMismatches <- function(x) {
            whichMismatches <- which(x[["values"]] == "?")
            if (length(whichMismatches) == 0) {
                value <- numeric(0)
            } else {
                start <- cumsum(x[["lengths"]])[whichMismatches]
                end <- start + x[["lengths"]][whichMismatches]
                value <-
                  eval(parse(text =
                             paste("c(", paste(start, ":", end, sep = "", collapse = ", "), ")")
                             ))
            }
            list(value)
        }
        getRange <- function(x) {
            if (!(x[["values"]][1] %in% c("-", "+"))) {
                start <- 1
            } else if (length(x[["values"]]) == 1) {
                start <- numeric(0)
            } else {
                start <- x[["lengths"]][1] + 1
            }
            if (!(x[["values"]][length(x[["values"]])] %in% c("-", "+"))) {
                end <- sum(x[["lengths"]])
            } else if (length(x[["values"]]) == 1) {
                end <- numeric(0)
            } else {
                end <- sum(x[["lengths"]][-length(x[["lengths"]])])
            }
            IRanges(start = start, end = end)
        }
        getIndels <- function(x) {
            if (x[["values"]][1] %in% c("-", "+")) {
                x[["values"]] <- x[["values"]][-1]
                x[["lengths"]] <- x[["lengths"]][-1]
            }
            if (x[["values"]][length(x[["values"]])] %in% c("-", "+")) {
                x[["values"]] <- x[["values"]][-length(x[["values"]])]
                x[["lengths"]] <- x[["lengths"]][-length(x[["lengths"]])]
            }
            isIndels <- (x[["values"]] %in% c("-", "+"))
            if (!any(isIndels))
                IRangesList(IRanges(numeric(0), numeric(0)))
            else
                IRangesList(IRanges(
                  cumsum(c(1, ifelse(isIndels, 0, x[["lengths"]])[-length(x[["lengths"]])]))[isIndels],
                                    width = x[["lengths"]][isIndels]))
        }
        if (length(pattern) != 1 || length(subject) != 1)
            stop("'pattern' and 'subject' must both be of length 1")
        if (nchar(pattern) != nchar(subject))
            stop("'pattern' and 'subject' must have the same number of characters")
        type <-
          match.arg(type,
                    c("global", "local", "overlap", "patternOverlap", "subjectOverlap"))
        gapOpening <- as.double(- abs(gapOpening))
        if (length(gapOpening) != 1 || is.na(gapOpening))
            stop("'gapOpening' must be a non-positive numeric vector of length 1")
        gapExtension <- as.double(- abs(gapExtension))
        if (length(gapExtension) != 1 || is.na(gapExtension))
            stop("'gapExtension' must be a non-positive numeric vector of length 1")

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
        substitutionArray <-
          array(unlist(substitutionMatrix, substitutionMatrix), dim = c(dim(substitutionMatrix), 2),
                dimnames = list(availableLetters, availableLetters, c("0", "1")))

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
        new("PairwiseAlignment",
            pattern =
            new("AlignedXStringSet",
                unaligned = XStringSet(baseClass, paste(degappedPattern, collapse = "")),
                range = getRange(patternRle), mismatch = getMismatches(patternRle),
                indel = getIndels(subjectRle)),
            subject =
            new("AlignedXStringSet",
                unaligned = XStringSet(baseClass, paste(degappedSubject, collapse = "")),
                range = getRange(subjectRle), mismatch = getMismatches(subjectRle),
                indel = getIndels(patternRle)),
            type = type,
            score =
              sum(substitutionMatrix[
                    match(explodedPattern[substitutionIndices], availableLetters) +
                      length(availableLetters) *
                        (match(explodedSubject[substitutionIndices], availableLetters) - 1)]) +
              gapOpening * sum(comparison[["values"]] %in% c("+", "-")) +
              gapExtension * sum(comparison[["lengths"]][comparison[["values"]] %in% c("+", "-")]),
            substitutionArray = substitutionArray, gapOpening = gapOpening, gapExtension = gapExtension)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Initialization.
###

setMethod("initialize", "PairwiseAlignment",
    function(.Object, pattern, subject, type, score, substitutionArray,
             gapOpening, gapExtension, check = TRUE)
    {
        if (!identical(class(unaligned(pattern)), class(unaligned(subject))))
            stop("'unaligned(pattern)' and 'unaligned(subject)' must be XString objects of the same subtype")
        if (length(type) != 1 || !(type %in% c("global", "local", "overlap", "patternOverlap", "subjectOverlap")))
            stop("'type' must be one of 'global', 'local', 'overlap', 'patternOverlap', or 'subjectOverlap'")
        if (length(pattern) != length(subject))
            stop("'length(pattern)' must equal 'length(subject)'")
        gapOpening <- as.double(- abs(gapOpening))
        if (length(gapOpening) != 1 || is.na(gapOpening))
            stop("'gapOpening' must be a non-positive numeric vector of length 1")
        gapExtension <- as.double(- abs(gapExtension))
        if (length(gapExtension) != 1 || is.na(gapExtension))
            stop("'gapExtension' must be a non-positive numeric vector of length 1")
        slot(.Object, "pattern", check = check) <- pattern
        slot(.Object, "subject", check = check) <- subject
        slot(.Object, "type", check = check) <- type
        slot(.Object, "score", check = check) <- score
        slot(.Object, "substitutionArray", check = check) <- substitutionArray
        slot(.Object, "gapOpening", check = check) <- gapOpening
        slot(.Object, "gapExtension", check = check) <- gapExtension
        .Object
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.PairwiseAlignment <- function(object)
{
    message <- character(0)
    if (!identical(class(unaligned(pattern(object))), class(unaligned(subject(object)))))
        message <- c(message, "'unaligned(pattern)' and 'unaligned(subject)' must be XString objects of the same subtype")
    if (length(object@type) != 1 || !(object@type %in% c("global", "local", "overlap", "patternOverlap", "subjectOverlap")))
        message <- c(message, "'type' must be one of 'global', 'local', 'overlap', 'patternOverlap', or 'subjectOverlap'")
    if (length(pattern) != length(subject))
        message <- c(message, "'length(pattern)' must equal 'length(subject)'")
    if (length(object@gapOpening) != 1 || is.na(object@gapOpening))
        message <- c(message, "'gapOpening' must be a non-positive numeric vector of length 1")
    if (length(object@gapExtension) != 1 || is.na(object@gapExtension))
        message <- c(message, "'gapExtension' must be a non-positive numeric vector of length 1")
    if (length(message) == 0)
        message <- NULL
    message
}

setValidity("PairwiseAlignment",
    function(object)
    {
        problems <- .valid.PairwiseAlignment(object)
        if (is.null(problems)) TRUE else problems
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods.
###

setMethod("pattern", "PairwiseAlignment", function(x) x@pattern)
setMethod("subject", "PairwiseAlignment", function(x) x@subject)

setGeneric("type", function(x) standardGeneric("type"))
setMethod("type", "PairwiseAlignment", function(x) x@type)
setMethod("type", "PairwiseAlignmentSummary", function(x) x@type)

setGeneric("score", function(x) standardGeneric("score"))
setMethod("score", "PairwiseAlignment", function(x) x@score)
setMethod("score", "PairwiseAlignmentSummary", function(x) x@score)

setMethod("nindel", "PairwiseAlignment",
          function(x) new("InDel", insertion = nindel(subject(x)), deletion = nindel(pattern(x))))
setMethod("nindel", "PairwiseAlignmentSummary",
          function(x) new("InDel", insertion = x@ninsertion, deletion = x@ndeletion))

setMethod("length", "PairwiseAlignment", function(x) length(score(x)))
setMethod("length", "PairwiseAlignmentSummary", function(x) length(score(x)))

setMethod("nchar", "PairwiseAlignment", function(x, type="chars", allowNA=FALSE) nchar(subject(x)))
setMethod("nchar", "PairwiseAlignmentSummary", 
          function(x, type="chars", allowNA=FALSE)
          unname(nmatch(x) + nmismatch(x) + x@ninsertion[,"WidthSum"] + x@ndeletion[,"WidthSum"]))

setMethod("alphabet", "PairwiseAlignment", function(x) alphabet(subject(x)))
setMethod("codec", "PairwiseAlignment", function(x) codec(subject(x)))

setGeneric("pid", signature="x", function(x, type="PID1") standardGeneric("pid"))
setMethod("pid", "PairwiseAlignment",
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

setMethod("Views", signature = c(subject = "PairwiseAlignment"),
          function(subject, start=NA, end=NA, names=NULL)
          {
              if (all(is.na(start)))
                  start <- start(subject(subject))
              else if (!is.numeric(start) || length(start) > 1)
                  stop("'start' must be either NA or an integer vector of length 1")
              else
                  start <- as.integer(start) + start(subject(subject))
              if (all(is.na(end)))
                  end <- end(subject(subject))
              else if (!is.numeric(end) || length(end) > 1)
                  stop("'end' must be either NA or an integer vector of length 1")
              else
                  end <- as.integer(end) + start(subject(subject))
              Views(super(unaligned(subject(subject))), start=start, end=end, names=names)
          })


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method.
###

### TODO: Make the "show" method to format the alignment in a SGD fashion
### i.e. split in 60-letter blocks and use the "|" character to highlight
### exact matches.
setMethod("show", "PairwiseAlignment", function(object)
          {
              if (length(object) == 0)
                  cat("Empty Pairwise Alignment\n")
              else {
                  cat(switch(type(object), "global" = "Global", "overlap" = "Overlap",
                             "patternOverlap" = "Pattern Overlap", "subjectOverlap" = "Subject Overlap",
                             "local" = "Local"), " Pairwise Alignment (1 of ", length(object), ")\n", sep = "")
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
          })


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "summary" method.
###

setMethod("summary", "PairwiseAlignment", function(object, weight=1L, ...)
          {
              if (!is.numeric(weight) || !(length(weight) %in% c(1, length(object))))
                  stop("'weight' must be an integer vector with length 1 or 'length(object)'")
              if (!is.integer(weight))
                weight <- as.integer(weight)
              if (all(weight == 1))
                  new("PairwiseAlignmentSummary",
                      type = type(object),
                      score = score(object),
                      nmatch = nmatch(object),
                      nmismatch = nmismatch(object),
                      ninsertion = nindel(subject(object)),
                      ndeletion = nindel(pattern(object)),
                      coverage = coverage(object),
                      mismatchSummary = mismatchSummary(object))
              else
                  new("PairwiseAlignmentSummary",
                      type = type(object),
                      score = rep(score(object), weight),
                      nmatch = rep(nmatch(object), weight),
                      nmismatch = rep(nmismatch(object), weight),
                      ninsertion = nindel(subject(object))[rep(seq_len(length(object)), weight), , drop = FALSE],
                      ndeletion = nindel(pattern(object))[rep(seq_len(length(object)), weight), , drop = FALSE],
                      coverage = coverage(object, weight = weight),
                      mismatchSummary = mismatchSummary(object, weight = weight))
          })

setMethod("show", "PairwiseAlignmentSummary", function(object)
          {
              cat(switch(type(object), "global" = "Global", "overlap" = "Overlap",
                         "patternOverlap" = "Pattern Overlap", "subjectOverlap" = "Subject Overlap",
                         "local" = "Local"), " Pairwise Alignment\n", sep = "")
              cat("Number of Alignments:  ", length(score(object)), "\n", sep = "")
              cat("\nScores:\n")
              print(summary(score(object)))
              cat("\nNumber of matches:\n")
              print(summary(nmatch(object)))
              n <- min(nrow(mismatchSummary(object)[["subject"]]), 10)
              cat(paste("\nTop", n, "Mismatch Counts:\n"))
              print(mismatchSummary(object)[["subject"]][
                order(mismatchSummary(object)[["subject"]][["Count"]],
                      mismatchSummary(object)[["subject"]][["Probability"]],
                      decreasing = TRUE)[seq_len(n)],,drop=FALSE])
          })


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "as.character" method.
###

setMethod("aligned", "PairwiseAlignment",
          function(x, degap = FALSE) {
              if (degap) {
                  value <- aligned(pattern(x), degap = degap)
              } else {
                  codecX <- codec(x)
                  if (is.null(codecX)) {
                      gapCode <- charToRaw("-")
                  } else {
                      letters2codes <- codecX@codes
                      names(letters2codes) <- codecX@letters
                      gapCode <- as.raw(letters2codes[["-"]])
                  }
                  value <-
                    .Call("PairwiseAlignment_align_aligned", x, gapCode, PACKAGE="Biostrings")
              }
              value
          })

setMethod("as.character", "PairwiseAlignment",
          function(x)
          {
              as.character(aligned(x))
          })

setMethod("toString", "PairwiseAlignment", function(x, ...) toString(as.character(x), ...))

setMethod("as.matrix", "PairwiseAlignment",
          function(x) {
              as.matrix(aligned(x))
          })


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

setMethod("[", "PairwiseAlignment",
    function(x, i, j, ..., drop)
    {
        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")
        if (missing(i) || (is.logical(i) && all(i)))
            return(x)
        if (is.logical(i))
            i <- which(i)
        if (!is.numeric(i) || any(is.na(i)))
            stop("invalid subsetting")
        if (any(i < 1) || any(i > length(x)))
            stop("subscript out of bounds")
        new("PairwiseAlignment",
            pattern = x@pattern[i],
            subject = x@subject[i],
            type = x@type,
            score = x@score[i],
            substitutionArray = x@substitutionArray,
            gapOpening = x@gapOpening,
            gapExtension = x@gapExtension)
    }
)

setReplaceMethod("[", "PairwiseAlignment",
    function(x, i, j,..., value)
    {
        stop("attempt to modify the value of a ", class(x), " instance")
    }
)

setMethod("rep", "PairwiseAlignment",
    function(x, times)
		x[rep.int(seq_len(length(x)), times)]
)
