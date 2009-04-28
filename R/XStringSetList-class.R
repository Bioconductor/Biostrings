### =========================================================================
### XStringSetList objects
### -------------------------------------------------------------------------
###
### Example:
###   unlisted <- DNAStringSet(c("AAA", "AC", "GGATA"))
###   partitioning <- Partitioning(c(0, 2, 2, 3))
###   x <- new("DNAStringSetList", unlisted=unlisted, partitioning=partitioning)
###

setClass("XStringSetList",
    contains="ListLike",
    representation(
        "VIRTUAL",
        unlisted="XStringSet",
        partitioning="Partitioning"
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

setGeneric("partitioning", function(x) standardGeneric("partitioning"))

setMethod("partitioning", "XStringSetList", function(x) x@partitioning)

setMethod("length", "XStringSetList", function(x) length(x@partitioning))

setMethod("names", "XStringSetList", function(x) names(x@partitioning))

setReplaceMethod("names", "XStringSetList",
    function(x, value)
    {
        names(x@partitioning) <- value
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Unsafe constructor (not exported). Use only when 'partitioning' is
### guaranteed to describe a valid Partitioning on 'unlisted'.
###

unsafe.newXStringSetList <- function(class, unlisted, partitioning)
{
    new2(class, unlisted=unlisted, partitioning=partitioning, check=FALSE)
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
        unsafe.newXStringSetList(ans_class, ans_unlisted, x@partitioning)
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
        ii <- x@partitioning[[i]]
        unlist(x)[ii]
    }
)

setReplaceMethod("[[", "XStringSetList",
    function(x, i, j,..., value)
    {
        stop("attempt to modify the value of a ", class(x), " instance")
    }
)

