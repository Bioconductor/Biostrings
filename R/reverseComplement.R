setGeneric(
    "reverse",
    function(x) standardGeneric("reverse")
)

setMethod("reverse", "BString",
    function(x)
    {
        x[length(x):1]
    }
)

setMethod("reverse", "BStringViews",
    function(x)
    {
        subject <- reverse(x@subject)
        ls <- subject@length
        lx <- length(x)
        first <- ls - x@last[lx:1] + 1
        last <- ls - x@first[lx:1] + 1
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
        old <- toString(DNA_STRING_CODEC@letters)
        new <- chartr("U", "T", toString(RNA_STRING_CODEC@letters))
        DNAString(chartr(old, new, toString(x)))
    }
)

setMethod("complement", "BStringViews",
    function(x)
    {
        subject <- complement(x@subject)
        views(subject, x@first, x@last)
    }
)

reverseComplement <- function(x)
{
    reverse(complement(x))
}

