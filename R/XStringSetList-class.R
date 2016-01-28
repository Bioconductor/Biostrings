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
### Going from XStringSet to XStringSetList with extractList() and family.
###

setMethod("relistToClass", "XStringSet",
    function(x) paste0(seqtype(x), "StringSetList")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor-like methods.
###

setMethod("nchar", "XStringSetList", IRanges:::nchar_CompressedList)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The XStringSetList() constructor. NOT exported.
###

.new_XStringSetList_from_list <- function(seqtype, x)
{
    x_eltNROWS <- elementNROWS(x)
    empty_idx <- which(x_eltNROWS == 0L)
    if (length(empty_idx) != 0L) {
        y <- x[-empty_idx]
    } else {
        y <- x
    }
    unlisted_y <- unlist(y, use.names=FALSE, recursive=FALSE)
    if (!is.list(unlisted_y) && length(unlisted_y) == sum(x_eltNROWS)) {
        unlisted_ans <- XStringSet(seqtype, unlisted_y)
    } else {
        ## In that case 'length(unlisted_y)' should be < 'sum(x_eltNROWS)'
        ## which means unlist() was not able to fully unlist 'y'. So let's
        ## try to turn each list element into an XStringSet object and then
        ## combine them together. This is of course much slower than if
        ## unlist() had succeeded.
        y <- lapply(unname(y), XStringSet, seqtype=seqtype)
        unlisted_ans <- do.call(c, y)
    }
    relist(unlisted_ans, x)
}

.new_XStringSetList_from_List <- function(seqtype, x)
{
    unlisted_x <- unlist(x, use.names=FALSE)
    unlisted_ans <- XStringSet(seqtype, unlisted_x)
    ans <- relist(unlisted_ans, x)
    ## relist() puts the names back but not the metadata columns.
    mcols(ans) <- mcols(x)
    ans
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
        ## Could be done elegantly with 'seqtype(unlisted(x)) <- value'
        ## if `unlisted<-` was available.
        unlisted_ans <- unlist(x, use.names=FALSE)
        seqtype(unlisted_ans) <- value
        ans <- relist(unlisted_ans, x)
        ## relist() puts the names back but not the metadata columns.
        mcols(ans) <- mcols(x)
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### User interface to the XStringSetList() constructor
###

BStringSetList <- function(..., use.names=TRUE)
    XStringSetList("B", ..., use.names=use.names)

DNAStringSetList <- function(..., use.names=TRUE)
    XStringSetList("DNA", ..., use.names=use.names)

RNAStringSetList <- function(..., use.names=TRUE)
    XStringSetList("RNA", ..., use.names=use.names)

AAStringSetList <- function(..., use.names=TRUE)
    XStringSetList("AA", ..., use.names=use.names)


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
setMethod("showAsCell", "XStringSetList",
     function(object) showAsCell(CharacterList(object))
)

