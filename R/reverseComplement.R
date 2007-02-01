### Dangerous, there is no guarantee that DNA_STRING_CODEC and RNA_STRING_CODEC
### are complementary. FIX AS SOON AS POSSIBLE!
.DNAComplementLookup <- function()
{
    lkup <- DNA_STRING_CODEC@dec_lkup
    lkup[lkup %in% letterAsByteVal("T")] <- letterAsByteVal("U")
    RNA_STRING_CODEC@enc_lkup[lkup + 1]
}

setGeneric("reverse", function(x, ...) standardGeneric("reverse"))

setMethod("reverse", "BString",
    function(x, ...)
    {
        lx <- length(x)
        data <- CharBuffer(lx)
        CharBuffer.reverseCopy(data, x@offset + 1, x@offset + lx, src=x@data)
        ## class(x) can be "BString", "DNAString", "RNAString" or "AAString"
        new(class(x), data)
    }
)

setMethod("reverse", "BStringViews",
    function(x, ...)
    {
        subject <- reverse(x@subject)
        ls <- subject@length
        start <- ls - x@end + 1
        end <- ls - x@start + 1
        views(subject, start, end)
    }
)

setGeneric(
    "complement",
    function(x) standardGeneric("complement")
)

setMethod("complement", "DNAString",
    function(x)
    {
        lx <- length(x)
        data <- CharBuffer(lx)
        lkup <- .DNAComplementLookup()
        CharBuffer.copy(data, x@offset + 1, x@offset + lx, src=x@data, lkup=lkup)
        DNAString(data)
    }
)

setMethod("complement", "BStringViews",
    function(x)
    {
        subject <- complement(x@subject)
        views(subject, x@start, x@end)
    }
)

setGeneric(
    "reverseComplement",
    function(x) standardGeneric("reverseComplement")
)

setMethod("reverseComplement", "DNAString",
    function(x)
    {
        lx <- length(x)
        data <- CharBuffer(lx)
        lkup <- .DNAComplementLookup()
        CharBuffer.reverseCopy(data, x@offset + 1, x@offset + lx, src=x@data, lkup=lkup)
        DNAString(data)
    }
)

setMethod("reverseComplement", "BStringViews",
    function(x)
    {
        subject <- reverseComplement(x@subject)
        ls <- subject@length
        start <- ls - x@end + 1
        end <- ls - x@start + 1
        views(subject, start, end)
    }
)

