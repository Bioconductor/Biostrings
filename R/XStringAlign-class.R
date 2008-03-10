### =========================================================================
### XStringAlign objects
### -------------------------------------------------------------------------
### An XStringAlign object contains the result of the alignment of 2 XString
### objects of the same subtype.


setClass("XStringAlign",
    representation(
        align1="XString",
        align2="XString",
        score="integer"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Initialization.
###

setMethod("initialize", "XStringAlign",
    function(.Object, align1, align2, score)
    {
        if (!identical(class(align1), class(align2)))
            stop("'align1' and 'align2' must be XString objects of the same subtype")
        if (length(align1) != length(align2))
            stop("'align1' and 'align2' must have the same length")
        if (length(score) != 1)
            stop("'score' must be a single integer")
        if (!is.integer(score))
            score <- as.integer(score)
        .Object@align1 <- align1
        .Object@align2 <- align2
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
        al1 <- align1(object)
        al2 <- align2(object)
        #l1 <- length(al1)
        #cat("  ", l1, "-letter \"", class(al1), "\" objects", sep="")
        cat("align1:", XString.get_snippet(al1, getOption("width") - 8))
        cat("\nalign2:", XString.get_snippet(al2, getOption("width") - 8))
        cat("\nscore:", score(object))
        cat("\n")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "as.character" method.
###

setMethod("as.character", "XStringAlign",
    function(x)
    {
        c(align1=as.character(x@align1), align2=as.character(x@align2))
    }
)

