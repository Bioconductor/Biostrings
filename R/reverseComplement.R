setGeneric(
    "reverseComplement",
    function(x) standardGeneric("reverseComplement")
)

setMethod("reverseComplement", "DNAString",
    function(x)
    {
        RNAString(x[length(x):1])
    }
)

setMethod("reverseComplement", "RNAString",
    function(x)
    {
        DNAString(x[length(x):1])
    }
)
