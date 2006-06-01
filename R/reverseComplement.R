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

setGeneric(
    "reverseComplement",
    function(x) standardGeneric("reverseComplement")
)

setMethod("reverseComplement", "DNAString",
    function(x)
    {
        reverse(complement(x))
    }
)

