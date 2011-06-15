### =========================================================================
### Sequence reversing and complementing
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "reverse" methods.
###

setMethod("reverse", "MaskedXString",
    function(x, ...)
    {
        x@unmasked <- reverse(unmasked(x))
        x@masks <- reverse(masks(x))
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "complement" generic and methods.
###

setGeneric("complement", signature="x",
    function(x, ...) standardGeneric("complement")
)

setMethod("complement", "DNAString",
    function(x, ...) xvcopy(x, lkup=getDNAComplementLookup())
)

setMethod("complement", "RNAString",
    function(x, ...) xvcopy(x, lkup=getRNAComplementLookup())
)

setMethod("complement", "DNAStringSet",
    function(x, ...) xvcopy(x, lkup=getDNAComplementLookup())
)

setMethod("complement", "RNAStringSet",
    function(x, ...) xvcopy(x, lkup=getRNAComplementLookup())
)

setMethod("complement", "XStringViews",
    function(x, ...)
    {
        x@subject <- complement(subject(x))
        x
    }
)

setMethod("complement", "MaskedDNAString",
    function(x, ...)
    {
        x@unmasked <- complement(unmasked(x))
        x
    }
)

setMethod("complement", "MaskedRNAString",
    function(x, ...)
    {
        x@unmasked <- complement(unmasked(x))
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "reverseComplement" generic and methods.
###
### We could just do this:
###   reverseComplement <- function(x) reverse(complement(x))
### But we want to perform only 1 copy of the sequence data!
### With the above implementation, reverseComplement(x) would copy the
### sequence data in 'x' twice: a first (temporary) copy to get the
### complement, followed by a second (final) copy to reverse it. Remember
### that the sequence data can be very big e.g. 250MB for Human chr1!
###

setGeneric("reverseComplement", signature="x",
    function(x, ...) standardGeneric("reverseComplement")
)

setMethod("reverseComplement", "DNAString",
    function(x, ...) xvcopy(x, lkup=getDNAComplementLookup(), reverse=TRUE)
)

setMethod("reverseComplement", "RNAString",
    function(x, ...) xvcopy(x, lkup=getRNAComplementLookup(), reverse=TRUE)
)

setMethod("reverseComplement", "DNAStringSet",
    function(x, ...) xvcopy(x, lkup=getDNAComplementLookup(), reverse=TRUE)
)

setMethod("reverseComplement", "RNAStringSet",
    function(x, ...) xvcopy(x, lkup=getRNAComplementLookup(), reverse=TRUE)
)

setMethod("reverseComplement", "XStringViews",
    function(x, ...)
    {
        x@subject <- reverseComplement(subject(x))
        x@ranges <- reverse(ranges(x), start=1L, end=length(subject(x)))
        x
    }
)

setMethod("reverseComplement", "MaskedDNAString",
    function(x, ...)
    {
        x@unmasked <- reverseComplement(unmasked(x))
        x@masks <- reverse(masks(x))
        x
    }
)

setMethod("reverseComplement", "MaskedRNAString",
    function(x, ...)
    {
        x@unmasked <- reverseComplement(unmasked(x))
        x@masks <- reverse(masks(x))
        x
    }
)

