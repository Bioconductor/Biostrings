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
        constantMatrix="matrix",
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
        coverage="numeric",
        mismatchSummary="list"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Initialization.
###

setMethod("initialize", "PairwiseAlignment",
    function(.Object, pattern, subject, type, score, constantMatrix,
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
        slot(.Object, "constantMatrix", check = check) <- constantMatrix
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

setMethod("views", signature = c(subject = "PairwiseAlignment"),
          function(subject, start=NA, end=NA)
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
              views(super(unaligned(subject(subject))), start = start, end = end)
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

setMethod("as.character", "PairwiseAlignment",
    function(x)
    {
        rbind(pattern = as.character(pattern(x)), subject = as.character(subject(x)))
    }
)

setMethod("toString", "PairwiseAlignment", function(x, ...) as.character(x))


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
            constantMatrix = x@constantMatrix,
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
