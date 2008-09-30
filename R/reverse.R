### =========================================================================
### The reverse() & related functions
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "reverse" methods.
###

setMethod("reverse", "XString",
    function(x, ...)
        XString.tr(x, reverse=TRUE)
)

setMethod("reverse", "XStringSet",
    function(x, ...)
    {
        x@super <- reverse(super(x))
        x@ranges <- reverse(x@ranges, start=1L, end=length(super(x)))
        x
    }
)

setMethod("reverse", "XStringViews",
    function(x, ...)
    {
        x@subject <- reverse(subject(x))
        callNextMethod(x, start=1L, end=length(subject(x)))
    }
)

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
    function(x, ...)
        XString.tr(x, lkup=getDNAComplementLookup())
)

setMethod("complement", "RNAString",
    function(x, ...)
        XString.tr(x, lkup=getRNAComplementLookup())
)

setMethod("complement", "DNAStringSet",
    function(x, ...)
    {
        x@super <- complement(super(x))
        x
    }
)

setMethod("complement", "RNAStringSet",
    function(x, ...)
    {
        x@super <- complement(super(x))
        x
    }
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
### Why don't we just do this:
###   reverseComplement <- function(x) reverse(complement(x))
### Because we want to perform only 1 copy of the sequence data!
### With the above implementation, reverseComplement(x) would copy the
### sequence data in 'x' twice: a first (temporary) copy to get the
### complement, followed by a second (final) copy to reverse it. Remember
### that the sequence data can be very big e.g. 250MB for Human chr1!
###

setGeneric("reverseComplement", signature="x",
    function(x, ...) standardGeneric("reverseComplement")
)

setMethod("reverseComplement", "DNAString",
    function(x, ...)
        XString.tr(x, lkup=getDNAComplementLookup(), reverse=TRUE)
)

setMethod("reverseComplement", "RNAString",
    function(x, ...)
        XString.tr(x, lkup=getRNAComplementLookup(), reverse=TRUE)
)

### For the following methods, doing this:
###   x@super <- reverseComplement(super(x))
###   x@ranges <- reverse(x@ranges, start=1L, end=length(super(x)))
### achieves our 1-copy goal and therefore is twice more efficient than doing
### this:
###   reverse(complement(x))
### or this:
###   x@super <- complement(super(x))
###   reverse(x)

.IRanges.reverse <- selectMethod("reverse", "IRanges")

setMethod("reverseComplement", "DNAStringSet",
    function(x, ...)
    {
        x@super <- reverseComplement(super(x))
        x@ranges <- reverse(x@ranges, start=1L, end=length(super(x)))
        x
    }
)

setMethod("reverseComplement", "RNAStringSet",
    function(x, ...)
    {
        x@super <- reverseComplement(super(x))
        x@ranges <- reverse(x@ranges, start=1L, end=length(super(x)))
        x
    }
)

setMethod("reverseComplement", "XStringViews",
    function(x, ...)
    {
        x@subject <- reverseComplement(subject(x))
        .IRanges.reverse(x, start=1L, end=length(subject(x)))
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Some convenience wrappers for different kinds of DNA <-> RNA
### transformations.
###

transcribe <- function(x)
{
    if (!is(x, "DNAString")) stop("transcribe() only works on DNA input")
    RNAString(complement(x))
}

cDNA <- function(x)
{
    if (!is(x, "RNAString")) stop("cDNA() only works on RNA input")
    DNAString(complement(x))
}

dna2rna <- function(x)
{
    if (!is(x, "DNAString")) stop("dna2rna() only works on DNA input")
    RNAString(x)
}
 
rna2dna <- function(x)
{
    if (!is(x, "RNAString")) stop("rna2dna() only works on RNA input")
    DNAString(x)
}

