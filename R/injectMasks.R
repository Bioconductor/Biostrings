### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "injectMasks" generic function and methods.
###

setGeneric("injectMasks", signature="x",
    function(x, letter="+") standardGeneric("injectMasks")
)

setMethod("injectMasks", "MaskedDNAString",
    function(x, letter="+")
    {
        stop("not ready yet")
    }
)

