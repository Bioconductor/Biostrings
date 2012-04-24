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

setGeneric("XStringSetList", signature="x",
    function(basetype, x, use.names=TRUE, ...)
        standardGeneric("XStringSetList")
)

.newCompressedList <- function(basetype, x, use.names)
{
    if (use.names)
        names <- names(x)
    if (basetype == "DNA")
        IRanges:::newCompressedList("DNAStringSetList", x, seq_len(length(x)), 
                                    names)
    else
        warning("XStringSetList currently supports basetype = 'DNA' only")
} 

.makeListOfXStringSets <- function(basetype, x) 
{
    if (basetype == "DNA")
        lapply(x, as, "DNAStringSet")
    else
        warning("XStringSetList currently supports basetype = 'DNA' only")
}

setMethod("XStringSetList", "list",
    function(basetype, x, use.names=TRUE, ...)
    {
        ok <- sapply(x, is, "XStringSet")
        if (!all(ok))
            x <- .makeListOfXStringSets(basetype, x)
        .newCompressedList(basetype, x, use.names)
    }
)

#setMethod("XStringSetList", "character",
#    function(basetype, x, use.names=TRUE)
#    {
#        x <- .charToXStringSet(basetype, x, start=NA,
#                                      width=NA, use.names=use.names)
#        callGeneric(basetype, x, use.names)
#    }
#)
#
#setMethod("XStringSetList", "XStringSet",
#    function(basetype, x, use.names=TRUE)
#    {
#        .newCompressedList(basetype, x, use.names)
#    }
#)

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

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The user interfaces to the XStringSetList() constructor.
###

DNAStringSetList <- function(..., use.names=TRUE)
                    {
                        listData <- list(...)
                        XStringSetList("DNA", listData, use.names=use.names)
                    }

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method.
###

setMethod("show", "XStringSetList",
    function(object)
    {
        cat(class(object), " of length ", length(object), "\n", sep = "")
        IRanges:::.showAtomicList(object, minLines=10)
    }
)

