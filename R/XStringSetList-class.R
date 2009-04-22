### =========================================================================
### XStringSetList objects
### -------------------------------------------------------------------------
###
### Example:
###   unlisted <- DNAStringSet(c("AAA", "AC", "GGATA"))
###   xsl <- new("DNAStringSetList", unlisted=unlisted, cum_eltlength=c(0L, 2L, 2L, 3L))
###

setClass("XStringSetList",
    contains="ListLike",
    representation(
        "VIRTUAL",
        unlisted="XStringSet",
        cum_eltlength="integer"
    )
)

setClass("BStringSetList",
    contains="XStringSetList",
    representation(
        unlisted="BStringSet"
    )
)
setClass("DNAStringSetList",
    contains="XStringSetList",
    representation(
        unlisted="DNAStringSet"
    )
)
setClass("RNAStringSetList",
    contains="XStringSetList",
    representation(
        unlisted="RNAStringSet"
    )
)
setClass("AAStringSetList",
    contains="XStringSetList",
    representation(
        unlisted="AAStringSet"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor-like methods.
###

setMethod("unlist", "XStringSetList",
    function(x, recursive=TRUE, use.names=TRUE) x@unlisted
)

setGeneric("cum_elementLengths", function(x) standardGeneric("cum_elementLengths"))
setMethod("cum_elementLengths", "XStringSetList", function(x) x@cum_eltlength)

setMethod("length", "XStringSetList", function(x) length(cum_elementLengths(x)))

setMethod("elementLengths", "XStringSetList",
    function(x) diff(c(0L, cum_elementLengths(x)))
)

setMethod("start", "XStringSetList",
    function(x)
    {
        if (length(x) == 0L)
            return(integer())
        c(1L, cum_elementLengths(x)[-length(x)] + 1L)
    }
)

setMethod("end", "XStringSetList", function(x) cum_elementLengths(x))

setMethod("width", "XStringSetList", function(x) elementLengths(x))

### XStringSetList objects don't support names for now.
setReplaceMethod("names", "XStringSetList",
    function(x, value)
    {
        if (!is.null(value))
            stop("attempt to put names on a ", class(x), " instance")
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Unsafe constructor (not exported). Use only when 'cum_eltlength' is
### guaranteed to describe valid adjacent ranges on 'unlisted' i.e. ranges
### that are within the limits of 'unlisted'.
###

unsafe.newXStringSetList <- function(class, unlisted, cum_eltlength)
{
    new2(class, unlisted=unlisted, cum_eltlength=cum_eltlength, check=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "xsbasetype" and "xsbasetype<-" methods.
###

setMethod("xsbasetype", "XStringSetList", function(x) xsbasetype(unlist(x)))

### Downgrades 'x' to a B/DNA/RNA/AAStringSetList instance!
setReplaceMethod("xsbasetype", "XStringSetList",
    function(x, value)
    {
        ## could be done with 'xsbasetype(unlisted(x)) <- value'
        ## if `unlisted<-` was available
        ans_class <- paste(value, "StringSetList", sep="")
        ans_unlisted <- unlist(x)
        xsbasetype(ans_unlisted) <- value
        unsafe.newXStringSetList(ans_class, ans_unlisted, x@cum_eltlength)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

setAs("XStringSetList", "IRanges",
    function(from)
        new2("IRanges", start=start(from), width=width(from), check=FALSE)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method.
###

setMethod("show", "XStringSetList",
    function(object)
    {
        cat("  A ", class(object), " instance of length ", length(object), "\n", sep="")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###
### No "[" method for now.
###

### Returns an XStringSet object of the same base type as 'x'.
setMethod("[[", "XStringSetList",
    function(x, i, j, ...)
    {
        i <- callNextMethod()
        ii <- seq_len(width(x)[i])
        if (i >= 2L)
            ii <- ii + end(x)[i - 1L]
        unlist(x)[ii]
    }
)

setReplaceMethod("[[", "XStringSetList",
    function(x, i, j,..., value)
    {
        stop("attempt to modify the value of a ", class(x), " instance")
    }
)

