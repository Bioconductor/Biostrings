### =========================================================================
### XStringSetList objects
### -------------------------------------------------------------------------
###

setClass("XStringSetList",
    contains="CompressedList",
    representation(
        "VIRTUAL",
        unlistData="XStringSet"
    ),
    prototype(
        elementType="XStringSet"
    )
)

setClass("BStringSetList",
    contains="XStringSetList",
    representation(
        unlistData="BStringSet"
    ),
    prototype(
        elementType="BStringSet"
    )
)
setClass("DNAStringSetList",
    contains="XStringSetList",
    representation(
        unlistData="DNAStringSet"
    ),
    prototype(
        elementType="DNAStringSet"
    )
)
setClass("RNAStringSetList",
    contains="XStringSetList",
    representation(
        unlistData="RNAStringSet"
    ),
    prototype(
        elementType="RNAStringSet"
    )
)
setClass("AAStringSetList",
    contains="XStringSetList",
    representation(
        unlistData="AAStringSet"
    ),
    prototype(
        elementType="AAStringSet"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor-like methods.
###

### TODO: Move the partitioning() generic to IRanges and make the
### "partitioning" method for XStringSetList objects below the method for
### for CompressedList objects.
setGeneric("partitioning", function(x) standardGeneric("partitioning"))

setMethod("partitioning", "XStringSetList", function(x) x@partitioning)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Unsafe constructor (not exported). Use only when 'partitioning' is
### guaranteed to describe a valid Partitioning on 'unlistData'.
###

unsafe.newXStringSetList <- function(class, unlistData, partitioning)
{
    new2(class, unlistData=unlistData, partitioning=partitioning, check=FALSE)
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
        ans_unlistData <- unlist(x)
        xsbasetype(ans_unlistData) <- value
        unsafe.newXStringSetList(ans_class, ans_unlistData, x@partitioning)
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
        i <- IRanges:::checkAndTranslateDbleBracketSubscript(x, i)
        ii <- x@partitioning[[i]]
        unlist(x)[ii]
    }
)

