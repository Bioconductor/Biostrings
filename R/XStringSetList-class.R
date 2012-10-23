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
    function(seqtype, x, use.names=TRUE, ...)
        standardGeneric("XStringSetList")
)

.newCompressedList <- function(seqtype, x, use.names)
{
    if (seqtype != "DNA")
        stop("XStringSetList() currently supports 'seqtype=\"DNA\"' only")
    if (!use.names)
        names(x) <- NULL
    IRanges:::newList("DNAStringSetList", x)
} 

.makeListOfXStringSets <- function(seqtype, x) 
{
    if (seqtype != "DNA")
        stop("XStringSetList() currently supports 'seqtype=\"DNA\"' only")
    lapply(x, as, "DNAStringSet")
}

setMethod("XStringSetList", "list",
    function(seqtype, x, use.names=TRUE, ...)
    {
        ok <- sapply(x, is, "XStringSet")
        if (!all(ok))
            x <- .makeListOfXStringSets(seqtype, x)
        .newCompressedList(seqtype, x, use.names)
    }
)

#setMethod("XStringSetList", "character",
#    function(seqtype, x, use.names=TRUE)
#    {
#        x <- .charToXStringSet(seqtype, x, start=NA,
#                                      width=NA, use.names=use.names)
#        callGeneric(seqtype, x, use.names)
#    }
#)
#
#setMethod("XStringSetList", "XStringSet",
#    function(seqtype, x, use.names=TRUE)
#    {
#        .newCompressedList(seqtype, x, use.names)
#    }
#)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "seqtype" and "seqtype<-" methods.
###

setMethod("seqtype", "XStringSetList", function(x) seqtype(unlist(x)))

### Downgrades 'x' to a B/DNA/RNA/AAStringSetList instance!
setReplaceMethod("seqtype", "XStringSetList",
    function(x, value)
    {
        ## could be done with 'seqtype(unlisted(x)) <- value'
        ## if `unlisted<-` was available
        ans_class <- paste(value, "StringSetList", sep="")
        ans_unlistData <- unlist(x)
        seqtype(ans_unlistData) <- value
        unsafe.newXStringSetList(ans_class, ans_unlistData, x@partitioning)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###
### No "[" method for now.
###

### Returns an XStringSet object of the same seqtype as 'x'.
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

### Display in a DataTable.
showAsCell <- IRanges:::showAsCell
setMethod("showAsCell", "XStringSetList",
     function(object) showAsCell(CharacterList(object))
)

