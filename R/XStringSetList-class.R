### =========================================================================
### XStringSetList objects
### -------------------------------------------------------------------------
###
### Example:
###   unlisted <- DNAStringSet(c("AAA", "AC", "GGATA"))
###   xsl <- new("DNAStringSetList", unlisted=unlisted, end=c(0L, 2L, 2L, 3L))
###

setClass("XStringSetList",
    contains=c("IPartitioning", "ListLike", "VIRTUAL"),
    representation(
        unlisted="XStringSet"
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Unsafe constructor (not exported). Use only when 'end' is guaranteed to
### describe a valid Partitioning on 'unlisted'.
###

unsafe.newXStringSetList <- function(class, unlisted, end)
{
    new2(class, unlisted=unlisted, end=end, check=FALSE)
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
        unsafe.newXStringSetList(ans_class, ans_unlisted, x@end)
    }
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
        i <- callNextMethod()  # calls "[[" method for ListLike objects
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

