transcribe = function(x) {
   if( !(is(x, "DNAString")) ) stop("transcribe only works on DNA input")
   RNAString(x)
}

cDNA = function(x) {
   if( !(is(x, "RNAString")) ) stop("cDNA only works on RNA input")
   DNAString(x)
}

dna2rna = function( x ) {
    lx = length(x)
    data <- XRaw(lx)
    lkup <- getDNAComplementLookup()
    XRaw.copy(data, x@offset + 1, x@offset + lx, src=x@data, lkup=lkup)
    RNAString(data)
}
 
rna2dna = function( x ) {
    lx = length(x)
    data <- XRaw(lx)
    lkup <- getRNAComplementLookup()
    XRaw.copy(data, x@offset + 1, x@offset + lx, src=x@data, lkup=lkup)
    DNAString(data)
}
 
setGeneric("reverse", signature="x",
    function(x, ...) standardGeneric("reverse")
)

setMethod("reverse", "BString",
    function(x, ...)
    {
        lx <- length(x)
        data <- XRaw(lx)
        XRaw.reverseCopy(data, x@offset + 1, x@offset + lx, src=x@data)
        ## class(x) can be "BString", "DNAString", "RNAString" or "AAString"
        new(class(x), data)
    }
)

setMethod("reverse", "BStringViews",
    function(x, ...)
    {
        subject <- reverse(x@subject)
        ls <- subject@length
        start <- ls - x@views$end + 1
        end <- ls - x@views$start + 1
        views(subject, start, end)
    }
)

setGeneric("complement", signature="x",
    function(x, ...) standardGeneric("complement")
)

setMethod("complement", "DNAString",
    function(x, ...)
    {
        lx <- length(x)
        data <- XRaw(lx)
        lkup <- getDNAComplementLookup()
        XRaw.copy(data, x@offset + 1, x@offset + lx, src=x@data, lkup=lkup)
        DNAString(data)
    }
)


setMethod("complement", "RNAString",
    function(x, ...)
    {
        lx <- length(x)
        data <- XRaw(lx)
        lkup <- getRNAComplementLookup()
        XRaw.copy(data, x@offset + 1, x@offset + lx, src=x@data, lkup=lkup)
        RNAString(data)
    }
)

setMethod("complement", "BStringViews",
    function(x, ...)
    {
        subject <- complement(x@subject)
        views(subject, x@views$start, x@views$end)
    }
)

setGeneric("reverseComplement", signature="x",
    function(x, ...) standardGeneric("reverseComplement")
)

setMethod("reverseComplement", "DNAString",
    function(x, ...)
    {
        lx <- length(x)
        data <- XRaw(lx)
        lkup <- getDNAComplementLookup()
        XRaw.reverseCopy(data, x@offset + 1, x@offset + lx, src=x@data, lkup=lkup)
        DNAString(data)
    }
)

setMethod("reverseComplement", "RNAString",
    function(x, ...)
    {
        lx <- length(x)
        data <- XRaw(lx)
        lkup <- getRNAComplementLookup()
        XRaw.reverseCopy(data, x@offset + 1, x@offset + lx, src=x@data, lkup=lkup)
        RNAString(data)
    }
)

setMethod("reverseComplement", "BStringViews",
    function(x, ...)
    {
        subject <- reverseComplement(x@subject)
        ls <- subject@length
        start <- ls - x@views$end + 1
        end <- ls - x@views$start + 1
        views(subject, start, end)
    }
)

