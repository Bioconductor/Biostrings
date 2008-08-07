### =========================================================================
### The reverse() generic & related functions
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "reverse" generic and methods.
###

setGeneric("reverse", signature="x",
    function(x, ...) standardGeneric("reverse")
)

### This method does NOT preserve normality.
### TODO: Move this function to the IRanges package.
.IRanges.reverse <- function(x, ...)
{
    args <- extraArgsAsList(NULL, ...)
    argnames <- names(args)
    n2p <- match(c("start", "end", "use.names"), argnames)
    if (is.na(n2p[1]))
        stop("'start' must be specified for \"reverse\" method for IRanges objects")
    start <- normargSingleStart(args[[n2p[1]]])
    if (is.na(n2p[2]))
        stop("'end' must be specified for \"reverse\" method for IRanges objects")
    end <- normargSingleEnd(args[[n2p[2]]])
    if (!is.na(n2p[3]) && !normargUseNames(args[[n2p[3]]])) {
        ## TEMPORARY WORKAROUND until .IRanges.reverse is moved to the IRanges package
        #unsafe.names(x) <- NULL
        x <- IRanges:::`unsafe.names<-`(x, NULL)
    }
    ## TEMPORARY WORKAROUND until .IRanges.reverse is moved to the IRanges package
    #unsafe.start(x) <- start + end - end(x)
    x <- IRanges:::`unsafe.start<-`(x, start + end - end(x))
    x
}

### TODO: Move this method to the IRanges package.
setMethod("reverse", "IRanges", .IRanges.reverse)

### TODO: Move this method to the IRanges package.
setMethod("reverse", "NormalIRanges",
    function(x, ...)
    {
        ## callNextMethod() temporarily breaks 'x' as a NormalIRanges object
        ## because the returned ranges are ordered from right to left.
        x <- callNextMethod()
        ## TEMPORARY WORKAROUND until this method is moved to the IRanges package
        #unsafe.update(x, start=rev(start(x)), width=rev(width(x)), names=rev(names(x)))
        IRanges:::unsafe.update(x, start=rev(start(x)), width=rev(width(x)), names=rev(names(x)))
    }
)

setMethod("reverse", "MaskCollection",
    function(x, ...)
    {
        start <- 1L
        end <- width(x)
        x@nir_list <- lapply(nir_list(x),
            function(nir) reverse(nir, start=start, end=end)
        )
        x
    }
)

setMethod("reverse", "XString",
    function(x, ...)
        XString.tr(x, reverse=TRUE)
)

setMethod("reverse", "XStringSet",
    function(x, ...)
    {
        x@super <- reverse(super(x))
        .IRanges.reverse(x, start=1L, end=length(super(x)))
    }
)

setMethod("reverse", "XStringViews",
    function(x, ...)
    {
        x@subject <- reverse(subject(x))
        .IRanges.reverse(x, start=1L, end=length(subject(x)))
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
        x@super <- reverseComplement(super(x))
        .IRanges.reverse(x, start=1L, end=length(super(x)))
    }
)

setMethod("reverseComplement", "RNAStringSet",
    function(x, ...)
    {
        x@super <- reverseComplement(super(x))
        .IRanges.reverse(x, start=1L, end=length(super(x)))
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

