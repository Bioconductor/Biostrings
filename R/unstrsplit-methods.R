### =========================================================================
### "unstrsplit" methods
### -------------------------------------------------------------------------


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
