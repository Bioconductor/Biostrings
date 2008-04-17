### =========================================================================
### XStringAlign objects
### -------------------------------------------------------------------------
### An XStringAlign object contains the result of the alignment of 2 XString
### objects of the same subtype.


setClass("XStringAlign",
    representation(
        align1="XString",
        align2="XString",
        type="character",
        matchScoring="matrix",
        gapOpening="integer",
        gapExtension="integer",
        score="integer"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Initialization.
###

setMethod("initialize", "XStringAlign",
    function(.Object, align1, align2, type, matchScoring, gapOpening, gapExtension, score)
    {
        if (!identical(class(align1), class(align2)))
            stop("'align1' and 'align2' must be XString objects of the same subtype")
        if (length(align1) != length(align2))
            stop("'align1' and 'align2' must have the same length")
        if (length(type) != 1 || !(type %in% c("global", "local", "overlap")))
            stop("'type' must be one of 'global', 'local', or 'overlap'")
        if (!is.matrix(matchScoring) || !is.integer(matchScoring))
            stop("'matchScoring' must be a matrix of integers")
        if (!identical(rownames(matchScoring), colnames(matchScoring)))
            stop("row and column names differ for matrix 'matchScoring'")
        if (is.null(rownames(matchScoring)))
            stop("matrix 'matchScoring' must have row and column names")
        if (any(duplicated(rownames(matchScoring))))
            stop("matrix 'matchScoring' has duplicated row names")
        gapOpening <- as.integer(- abs(gapOpening))
        if (is.na(gapOpening) || length(gapOpening) != 1)
            stop("'gapOpening' must be a non-positive integer vector of length 1")
        gapExtension <- as.integer(- abs(gapExtension))
        if (is.na(gapExtension) || length(gapExtension) != 1)
            stop("'gapExtension' must be a non-positive integer vector of length 1")
        .Object@align1 <- align1
        .Object@align2 <- align2
        .Object@type <- type
        .Object@matchScoring <- matchScoring
        .Object@gapOpening <- gapOpening
        .Object@gapExtension <- gapExtension
        .Object@score <- score
        .Object
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods.
###

setGeneric("align1", function(x) standardGeneric("align1"))
setMethod("align1", "XStringAlign", function(x) x@align1)

setGeneric("align2", function(x) standardGeneric("align2"))
setMethod("align2", "XStringAlign", function(x) x@align2)

setGeneric("type", function(x) standardGeneric("type"))
setMethod("type", "XStringAlign", function(x) x@type)

setGeneric("score", function(x) standardGeneric("score"))
setMethod("score", "XStringAlign", function(x) x@score)

setMethod("length", "XStringAlign", function(x) length(align1(x)))
setMethod("nchar", "XStringAlign", function(x, type="chars", allowNA=FALSE) length(x))

setMethod("alphabet", "XStringAlign", function(x) alphabet(align1(x)))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method.
###

### TODO: Make the "show" method to format the alignment in a SGD fashion
### i.e. split in 60-letter blocks and use the "|" character to highlight
### exact matches.
setMethod("show", "XStringAlign",
    function(object)
    {
        cat(switch(type(object), "global" = "Global", "overlap" = "Overlap",
                   "local" = "Local"), "Pairwise Alignment\n")
        cat("1: ", toSeqSnippet(align1(object), getOption("width") - 8), "\n")
        cat("2: ", toSeqSnippet(align2(object), getOption("width") - 8), "\n")
        cat("Score: ", score(object), "\n")
        cat("Gap Opening: ", object@gapOpening, "\n")
        cat("Gap Extension: ", object@gapExtension, "\n")
        cat("Match/Mismatch Scoring:\n")
        print(object@matchScoring)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "as.character" method.
###

setMethod("as.character", "XStringAlign",
    function(x)
    {
        c(align1 = as.character(align1(x)), align2 = as.character(align2(x)))
    }
)

