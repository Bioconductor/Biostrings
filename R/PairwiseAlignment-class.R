### ==========================================================================
### XStringPairwiseAlignment objects
### --------------------------------------------------------------------------
### An XStringPairwiseAlignment object contains the result of the alignment of
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

setGeneric("score", function(x) standardGeneric("score"))
setMethod("score", "PairwiseAlignment", function(x) x@score)

setMethod("length", "PairwiseAlignment", function(x) length(subject(x)))
setMethod("nchar", "PairwiseAlignment", function(x, type="chars", allowNA=FALSE) nchar(subject(x)))

setMethod("alphabet", "PairwiseAlignment", function(x) alphabet(subject(x)))
setMethod("codec", "PairwiseAlignment", function(x) codec(subject(x)))

setGeneric("mismatchTable", signature = "x",
           function(x, shiftLeft = 0L, shiftRight = 0L)
           standardGeneric("mismatchTable"))
setMethod("mismatchTable", "PairwiseAlignment",
          function(x, shiftLeft = 0L, shiftRight = 0L)
          {
              if (!isSingleNumber(shiftLeft) || shiftLeft > 0)
                  stop("'shiftLeft' must be a non-positive integer")
              if (!isSingleNumber(shiftRight) || shiftRight < 0)
                  stop("'shiftRight' must be a non-negative integer")
              shiftLeft <- as.integer(shiftLeft)
              shiftRight <- as.integer(shiftRight)
              nMismatch <- nmismatch(x)
              patternNumber <- rep.int(1:length(nMismatch), nMismatch)
              patternSubset <- unaligned(pattern(x))[patternNumber]
              subjectSubset <- unaligned(subject(x))[[1]]
              patternPosition <- unlist(mismatch(pattern(x)))
              subjectPosition <- unlist(mismatch(subject(x)))
              if (shiftLeft == 0L) {
                  patternStart <- patternPosition
                  subjectStart <- subjectPosition
              } else {
                  patternStart <- pmax(patternPosition + shiftLeft, 1L)
                  subjectStart <- pmax(subjectPosition + shiftLeft, 1L)
              }
              if (shiftRight == 0L) {
                  patternEnd <- patternPosition
                  subjectEnd <- subjectPosition
              } else {
                  patternEnd <- pmin(patternPosition + shiftRight, width(patternSubset))
                  subjectEnd <- pmin(subjectPosition + shiftRight, width(subjectSubset))
              }
              patternSubstring <-
                narrow(patternSubset, start = patternStart, end = patternEnd)
              subjectSubstring <-
                views(subjectSubset, start = subjectStart, end = subjectEnd)
              if (length(x@constantMatrix) > 0) {
                  output <-
                    data.frame("PatternNumber" = patternNumber,
                               "PatternStart" = patternStart,
                               "PatternEnd" = patternEnd,
                               "SubjectStart" = subjectStart,
                               "SubjectEnd" = subjectEnd,
                               "PatternSubstring" = as.character(patternSubstring),
                               "SubjectSubstring" = as.character(subjectSubstring))
              } else {
                  if (length(pattern(x)@quality) == 1 && width(pattern(x)@quality) == 1)
                      patternQuality <-
                        views(pattern(x)@quality[[1]][rep.int(1L, max(patternEnd))],
                              start = patternStart, end = patternEnd)
                  else
                      patternQuality <-
                        narrow(pattern(x)@quality[patternNumber],
                               start = patternStart, end = patternEnd)
                  if (length(subject(x)@quality) == 1 && width(subject(x)@quality) == 1)
                      subjectQuality <-
                        views(subject(x)@quality[[1]][rep.int(1L, max(subjectEnd))],
                              start = subjectStart, end = subjectEnd)
				  else
                      subjectQuality <-
                        views(subject(x)@quality[[1]], start = subjectStart, end = subjectEnd)
                  output <-
                    data.frame("PatternNumber" = patternNumber,
                               "PatternStart" = patternStart,
                               "PatternEnd" = patternEnd,
                               "SubjectStart" = subjectStart,
                               "SubjectEnd" = subjectEnd,
                               "PatternSubstring" = as.character(patternSubstring),
                               "SubjectSubstring" = as.character(subjectSubstring),
                               "PatternQuality" = as.character(patternQuality),
                               "SubjectQuality" = as.character(subjectQuality))
              }
              output
          })

setMethod("nmismatch", c(pattern = "PairwiseAlignment", x = "missing"),
          function(pattern, x, fixed) nmismatch(pattern(pattern)))

setMethod("coverage", "PairwiseAlignment", function(x, start = NA, end = NA, weight = 1L)
          {
              if (any(is.na(start)))
                  start <- 1
              if (any(is.na(end)))
                  end <- nchar(super(unaligned(subject(x))))
              coverage(subject(x)@range, start = start, end = end, weight = weight)
          })

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
              show(pattern(object))
              cat(paste(c("subject: ", rep(" ", subjectSpaces)), collapse = ""))
              show(subject(object))
              cat("score:", score(object)[1], "\n")
          })


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "summary" method.
###

setMethod("summary", "PairwiseAlignment", function(object, ...)
          {
              totalIndel <- function(object) {
                  unlist(lapply(object, function(x) sum(width(x))))
              }
              cat(switch(type(object), "global" = "Global", "overlap" = "Overlap",
                         "patternOverlap" = "Pattern Overlap", "subjectOverlap" = "Subject Overlap",
                         "local" = "Local"), " Pairwise Alignment\n", sep = "")
              cat("Number of Alignments:  ", length(object), "\n", sep = "")
              cat("\nScores:\n")
              print(summary(score(object)))
              cat("\nNumber of matched characters:\n")
              print(summary(nchar(object) - nmismatch(object) - totalIndel(indel(subject(object))) - 
                            totalIndel(indel(pattern(object)))))
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
    {
        x[rep.int(1:length(x), times)]
    }
)
