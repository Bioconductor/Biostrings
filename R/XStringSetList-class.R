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

partitioning <- function(x) .Defunct("PartitioningByEnd")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The XStringSetList() constructor. NOT exported.
###

.new_XStringSetList_from_List <- function(seqtype, x)
{
    unlisted_x <- unlist(x, use.names=FALSE)
    unlisted_ans <- XStringSet(seqtype, unlisted_x)
    relist(unlisted_ans, x)
}

.new_XStringSetList_from_list <- function(seqtype, x)
{
    x_eltlens <- elementLengths(x)
    empty_idx <- which(x_eltlens == 0L)
    if (length(empty_idx) != 0L) {
        y <- x[-empty_idx]
    } else {
        y <- x
    }
    unlisted_y <- unlist(y, use.names=FALSE, recursive=FALSE)
    if (length(unlisted_y) < sum(x_eltlens)) {
        ## unlist() was not able to fully unlist 'y'. So let's try to turn
        ## each list element into an XStringSet object and then combine
        ## them together.
        y <- lapply(unname(y), XStringSet, seqtype=seqtype)
        unlisted_ans <- do.call(c, y)
    } else {
        unlisted_ans <- XStringSet(seqtype, unlisted_y)
    }
    relist(unlisted_ans, x)
}

XStringSetList <- function(seqtype, ..., use.names=TRUE)
{
    if (!isTRUEorFALSE(use.names))
        stop("'use.names' must be TRUE or FALSE")
    x <- list(...)
    if (length(x) == 1L) {
        x1 <- x[[1L]]
        if (is.list(x1) || (is(x1, "List") && !is(x1, "XStringSet"))) {
            x <- x1
            if (is(x, "List")) {
                if (!use.names)
                    names(x) <- NULL
                return(.new_XStringSetList_from_List(seqtype, x))
            }
        }
    }
    if (!use.names)
        names(x) <- NULL
    .new_XStringSetList_from_list(seqtype, x)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "seqtype" and "seqtype<-" methods.
###

setMethod("seqtype", "XStringSetList",
    function(x) seqtype(unlist(x, use.names=FALSE))
)

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
        unlist(x, use.names=FALSE)[ii]
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### User interface to the XStringSetList() constructor
###

DNAStringSetList <- function(..., use.names=TRUE)
{
    XStringSetList("DNA", ..., use.names=use.names)
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

