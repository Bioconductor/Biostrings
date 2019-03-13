### =========================================================================
### strsplit() / unstrsplit()
### -------------------------------------------------------------------------


.strsplit.useAsDefault <- function(x, ...) base::strsplit(x, ...)
setGeneric("strsplit", signature="x",
    function(x, ...) standardGeneric("strsplit"),
    useAsDefault=.strsplit.useAsDefault
)

setMethod("strsplit", "XStringSet",
    function(x, split, ...)
    {
        mi <- vmatchPattern(split, x, ...)
        at <- gaps(as(mi, "CompressedIRangesList"),
                   start=1L, end=width(x))
        extractAt(x, at)
    }
)

setMethod("unstrsplit", "XStringSetList",
    function(x, sep="")
    {
        x_seqtype <- seqtype(x)
        sep <- XString(x_seqtype, sep)
        .Call("XStringSetList_unstrsplit", x, sep, x_seqtype,
              PACKAGE="Biostrings")
    }
)

setMethod("unstrsplit", "XStringSet",
    function(x, sep="")
    {
        x
    }
)

