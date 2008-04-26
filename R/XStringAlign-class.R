### =========================================================================
### XStringAlign objects
### -------------------------------------------------------------------------
### An XStringAlign object contains the result of the alignment of 2 XString
### objects of the same subtype.


setClass("XStringAlign",
    representation(
        align1="XString",
        align2="XString",
        quality1="BString",
        quality2="BString",
        type="character",
        score="numeric",
        constantMatrix="matrix",
        gapOpening="numeric",
        gapExtension="numeric",
        endGap="logical"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Initialization.
###

setMethod("initialize", "XStringAlign",
    function(.Object, align1, align2, quality1, quality2, type, score, constantMatrix, gapOpening, gapExtension, endGap)
    {
        if (!identical(class(align1), class(align2)))
            stop("'align1' and 'align2' must be XString objects of the same subtype")
        if (length(align1) != length(align2))
            stop("'align1' and 'align2' must have the same length")
        if (length(type) != 1 || !(type %in% c("global", "local", "overlap")))
            stop("'type' must be one of 'global', 'local', or 'overlap'")
        gapOpening <- as.double(- abs(gapOpening))
        if (is.na(gapOpening) || length(gapOpening) != 1)
            stop("'gapOpening' must be a non-positive numeric vector of length 1")
        gapExtension <- as.double(- abs(gapExtension))
        if (is.na(gapExtension) || length(gapExtension) != 1)
            stop("'gapExtension' must be a non-positive numeric vector of length 1")
        endGap <- as.logical(endGap)
        if (length(endGap) != 1 || any(is.na(endGap)))
            stop("'endGap' must be a non-missing logical value")
		.Object@align1 <- align1
        .Object@align2 <- align2
        .Object@quality1 <- quality1
        .Object@quality2 <- quality2
        .Object@type <- type
        .Object@score <- score
        .Object@constantMatrix <- constantMatrix
        .Object@gapOpening <- gapOpening
        .Object@gapExtension <- gapExtension
        .Object@endGap <- endGap
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

