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
        callNextMethod(x, start=1L, end=length(super(x)))
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

setMethod("reverseComplement", "DNAStringSet",
    function(x, ...)
    {
        x@super <- complement(super(x))
		reverse(x, start=1L, end=length(super(x)))
    }
)

setMethod("reverseComplement", "RNAStringSet",
    function(x, ...)
    {
        x@super <- complement(super(x))
		reverse(x, start=1L, end=length(super(x)))
    }
)

setMethod("reverseComplement", "XStringViews",
    function(x, ...)
    {
        x@subject <- complement(subject(x))
		reverse(x, start=1L, end=length(subject(x)))
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

