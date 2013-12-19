### =========================================================================
### "unstrsplit" methods
### -------------------------------------------------------------------------


setMethod("unstrsplit", "XStringSetList",
    function(x, sep="-")
    {
        x_seqtype <- seqtype(x)
        sep <- XString(x_seqtype, sep)
        .Call("XStringSetList_unstrsplit", x, sep, x_seqtype,
              PACKAGE="Biostrings")
    }
)

### We want this method to use the same 'sep' default as the method for list
### objects defined in IRanges.
setMethod("unstrsplit", "BStringSetList",
    function(x, sep=",") callNextMethod(x, sep)
)

