DNAComplementHash <- function()
{
    old <- DNAString(toString(DNA_STRING_CODEC@letters))
    new <- DNAString(chartr("U", "T", toString(RNA_STRING_CODEC@letters)))
    buildCodecHashTable(old@data[], new@data[], 0)
}

setGeneric(
    "reverse",
    function(x) standardGeneric("reverse")
)

setMethod("reverse", "BString",
    function(x)
    {
        lx <- length(x)
        data <- CharBuffer(lx)
        CharBuffer.reverseCopy(data, x@offset + 1, x@offset + lx, src=x@data)
        # class(x) can be "BString", "DNAString" or "RNAString"
        new(class(x), data)
    }
)

setMethod("reverse", "BStringViews",
    function(x)
    {
        subject <- reverse(x@subject)
        ls <- subject@length
        first <- ls - x@last + 1
        last <- ls - x@first + 1
        views(subject, first, last)
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
        hash <- DNAComplementHash()
        CharBuffer.copy(data, x@offset + 1, x@offset + lx, src=x@data, hash=hash)
        DNAString(data)
    }
)

setMethod("complement", "BStringViews",
    function(x)
    {
        subject <- complement(x@subject)
        views(subject, x@first, x@last)
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
        hash <- DNAComplementHash()
        CharBuffer.reverseCopy(data, x@offset + 1, x@offset + lx, src=x@data, hash=hash)
        DNAString(data)
    }
)

setMethod("reverseComplement", "BStringViews",
    function(x)
    {
        subject <- reverseComplement(x@subject)
        ls <- subject@length
        first <- ls - x@last + 1
        last <- ls - x@first + 1
        views(subject, first, last)
    }
)

