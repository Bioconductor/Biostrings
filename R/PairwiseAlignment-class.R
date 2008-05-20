### ==========================================================================
### XStringPairwiseAlignment objects
### --------------------------------------------------------------------------
### An XStringPairwiseAlignment object contains the result of the alignment of
### two XString objects of the same subtype.


setClass("PairwiseAlignment",
    representation(
        pattern="AlignedXString",
        subject="AlignedXString",
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
        if (is.na(gapOpening) || length(gapOpening) != 1)
            stop("'gapOpening' must be a non-positive numeric vector of length 1")
        gapExtension <- as.double(- abs(gapExtension))
        if (is.na(gapExtension) || length(gapExtension) != 1)
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
    if (is.na(object@gapOpening) || length(object@gapOpening) != 1)
        message <- c(message, "'gapOpening' must be a non-positive numeric vector of length 1")
    if (is.na(object@gapExtension) || length(object@gapExtension) != 1)
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
setMethod("nchar", "PairwiseAlignment", function(x, type="chars", allowNA=FALSE) length(x))

setMethod("alphabet", "PairwiseAlignment", function(x) alphabet(subject(x)))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method.
###

### TODO: Make the "show" method to format the alignment in a SGD fashion
### i.e. split in 60-letter blocks and use the "|" character to highlight
### exact matches.
setMethod("show", "PairwiseAlignment",
    function(object)
    {
        cat(switch(type(object), "global" = "Global", "overlap" = "Overlap",
                   "patternOverlap" = "Pattern Overlap", "subjectOverlap" = "Subject Overlap",
                   "local" = "Local"), "Pairwise Alignment\n")
        patternSpaces <- floor(log10(start(subject(object)))) - floor(log10(start(pattern(object))))
		subjectSpaces <- max(0, - patternSpaces)
		patternSpaces <- max(0, patternSpaces)
		cat(paste("pattern: ", rep(" ", patternSpaces), sep = "", collapse = ""))
        show(pattern(object))
        cat(paste("subject: ", rep(" ", subjectSpaces), sep = "", collapse = ""))
        show(subject(object))
        cat("score:", score(object), "\n")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "as.character" method.
###

setMethod("as.character", "PairwiseAlignment",
    function(x)
    {
        c(pattern = as.character(pattern(x)), subject = as.character(subject(x)))
    }
)
