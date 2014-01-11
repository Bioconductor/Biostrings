### =========================================================================
### MaskedXString objects
### -------------------------------------------------------------------------
###
### Note that thinking of MaskedXString objects as particular XString objects
### might seem natural and therefore making the MaskedXString class an
### extension of the XString class is tempting... But in the end there would
### be no benefit in doing so because inheritance would not give us
### anything good out-of-the-box. This is because almost all XString methods
### are doing the wrong thing for MaskedXString objects and hence would need
### to be overwritten so it would not save us any work. Furthermore, some of
### the XString methods should be disabled for MaskedXString objects. But if
### MaskedXString objects are considered XString objects, then these methods
### need to be overwritten and call stop() so that they fail on MaskedXString
### objects. By not making the MaskedXString class an extension of the XString
### class, we avoid all this mess.
###

### Not good (MaskedXString extends XString).
#setClass("MaskedXString",
#    contains="XString",
#    representation(
#        "VIRTUAL",
#        masks="MaskCollection"
#    )
#)
### 4 direct "MaskedXString" derivations (no additional slot)
#setClass("MaskedBString", contains=c("MaskedXString", "BString"))
#setClass("MaskedDNAString", contains=c("MaskedXString", "DNAString"))
#setClass("MaskedRNAString", contains=c("MaskedXString", "RNAString"))
#setClass("MaskedAAString", contains=c("MaskedXString", "AAString"))

### Better (MaskedXString does NOT extend XString).
setClass("MaskedXString",
    representation(
        "VIRTUAL",
        unmasked="XString",
        masks="MaskCollection"
    )
)

setClass("MaskedBString",
    contains="MaskedXString",
    representation(
        unmasked="BString"
    )
)
setClass("MaskedDNAString",
    contains="MaskedXString",
    representation(
        unmasked="DNAString"
    )
)
setClass("MaskedRNAString",
    contains="MaskedXString",
    representation(
        unmasked="RNAString"
    )
)
setClass("MaskedAAString",
    contains="MaskedXString",
    representation(
        unmasked="AAString"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor-like methods.
###

setGeneric("unmasked", function(x) standardGeneric("unmasked"))
setMethod("unmasked", "MaskedXString", function(x) x@unmasked)
setMethod("unmasked", "XString", function(x) x)  # no-op

setGeneric("masks", function(x) standardGeneric("masks"))
setMethod("masks", "XString", function(x) NULL)
setMethod("masks", "MaskedXString", function(x) x@masks)

setMethod("length", "MaskedXString", function(x) length(unmasked(x)))

setMethod("maskedwidth", "MaskedXString", function(x) maskedwidth(collapse(masks(x))))

setMethod("maskedratio", "MaskedXString", function(x) maskedratio(collapse(masks(x))))

setMethod("nchar", "MaskedXString",
    function(x, type="chars", allowNA=FALSE)
    {
        length(x) - maskedwidth(x)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.MaskedXString.unmasked <- function(object)
{
    if (!is(unmasked(object), "XString"))
        return("the 'unmasked' slot must contain an XString object")
    if (!is(object, paste("Masked", xsbaseclass(object), sep="")))
        return("bad XString base type for the unmasked sequence")
    if (length(object) != width(masks(object)))
        return("the length of the object and the width of the mask collection differ")
    NULL
}

.valid.MaskedXString.masks <- function(object)
{
    masks <- masks(object)
    if (!is(masks, "MaskCollection"))
        return("the 'masks' slot must contain a MaskCollection object")
    if (width(masks) != length(object))
        return("the length of the object and the width of the mask collection differ")
    NULL
}

.valid.MaskedXString <- function(object)
{
    c(.valid.MaskedXString.unmasked(object),
      .valid.MaskedXString.masks(object))
}

setValidity("MaskedXString",
    function(object)
    {
        problems <- .valid.MaskedXString(object)
        if (is.null(problems)) TRUE else problems
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "seqtype" and "seqtype<-" methods.
###

setMethod("seqtype", "MaskedXString", function(x) seqtype(unmasked(x)))

### Downgrades 'x' to a MaskedB/DNA/RNA/AAString instance!
setReplaceMethod("seqtype", "MaskedXString",
    function(x, value)
    {
        ## could be done with 'seqtype(unmasked(x)) <- value'
        ## if `unmasked<-` was available
        unmasked <- unmasked(x)
        seqtype(unmasked) <- value
        ans_class <- paste("Masked", value, "String", sep="")
        new(ans_class, unmasked=unmasked, masks=masks(x))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

### From MaskedXString objects to MaskedXString objects.
setAs("MaskedXString", "MaskedBString",
    function(from) {seqtype(from) <- "B"; from}
)
setAs("MaskedXString", "MaskedDNAString",
    function(from) {seqtype(from) <- "DNA"; from}
)
setAs("MaskedXString", "MaskedRNAString",
    function(from) {seqtype(from) <- "RNA"; from}
)
setAs("MaskedXString", "MaskedAAString",
    function(from) {seqtype(from) <- "AA"; from}
)

### From XString objects to MaskedXString objects.
setAs("BString", "MaskedBString",
    function(from)
    {
        masks <- new("MaskCollection", width=length(from))
        new("MaskedBString", unmasked=from, masks=masks)
    }
)
setAs("DNAString", "MaskedDNAString",
    function(from)
    {
        masks <- new("MaskCollection", width=length(from))
        new("MaskedDNAString", unmasked=from, masks=masks)
    }
)
setAs("RNAString", "MaskedRNAString",
    function(from)
    {
        masks <- new("MaskCollection", width=length(from))
        new("MaskedRNAString", unmasked=from, masks=masks)
    }
)
setAs("AAString", "MaskedAAString",
    function(from)
    {
        masks <- new("MaskCollection", width=length(from))
        new("MaskedAAString", unmasked=from, masks=masks)
    }
)

### From MaskedXString objects to XString objects.
setAs("MaskedBString", "BString",
    function(from) unmasked(from)
)
setAs("MaskedDNAString", "DNAString",
    function(from) unmasked(from)
)
setAs("MaskedRNAString", "RNAString",
    function(from) unmasked(from)
)
setAs("MaskedAAString", "AAString",
    function(from) unmasked(from)
)

### Dispatch on 'x' (see generic in XString-class.R).
setMethod("XString", "MaskedXString",
    function(seqtype, x, start=NA, end=NA, width=NA)
        XString(seqtype, unmasked(x), start=start, end=end, width=width)
)

### From a MaskedXString object to a MaskCollection object.
setAs("MaskedXString", "MaskCollection",
    function(from) masks(from)
)

### From a MaskedXString object to a NormalIRanges object.
setAs("MaskedXString", "NormalIRanges",
    function(from) as(masks(from), "NormalIRanges")
)

### From a MaskedXString object to an XStringViews object.
setAs("MaskedXString", "XStringViews",
    function(from)
    {
        views <- gaps(collapse(masks(from)))[[1]]
        unsafe.newXStringViews(unmasked(from), start(views), width(views))
    }
)

setAs("MaskedXString", "Views", function(from) as(from, "XStringViews"))

### NOT exported.
toXStringViewsOrXString <- function(x)
{
    x0 <- unmasked(x)
    mask1 <- collapse(masks(x))
    if (isEmpty(mask1))
        return(x0)
    views <- gaps(mask1)[[1]]
    unsafe.newXStringViews(x0, start(views), width(views))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The transformation methods (endomorphisms) "collapse" and "gaps".
###

setMethod("collapse", "MaskedXString",
    function(x)
    {
        x@masks <- collapse(masks(x))
        x
    }
)

### 'start' and 'end' are ignored.
setMethod("gaps", "MaskedXString",
    function(x, start=NA, end=NA)
    {
        x@masks <- gaps(collapse(masks(x)))
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "subseq" method for MaskedXString objects.
###
### This method is used by the toSeqSnippet() function when called
### on a MaskedXString object.
###

setMethod("subseq", "MaskedXString",
    function(x, start=NA, end=NA, width=NA)
    {
        x@unmasked <- subseq(unmasked(x), start=start, end=end, width=width)
        x@masks <- narrow(masks(x), start=start, end=end, width=width)
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "as.character" and "toString" methods.
###

### This method is used by the toSeqSnippet() function when called
### on a MaskedXString object.
setMethod("as.character", "MaskedXString",
    function(x)
    {
        ans <- as.character(unmasked(x))
        nir0 <- as(x, "NormalIRanges")
        for (i in seq_len(length(nir0))) {
            strip <- paste(rep.int("#", width(nir0)[i]), collapse="")
            substr(ans,  start(nir0)[i], end(nir0)[i]) <- strip
        }
        ans
    }
)

setMethod("toString", "MaskedXString", function(x, ...) as.character(x))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method.
###

setMethod("show", "MaskedXString",
    function(object)
    {
        lo <- length(object)
        cat("  ", lo, "-letter \"", class(object), "\" instance (# for masking)\n", sep="")
        cat("seq:", toSeqSnippet(object, getOption("width") - 5))
        cat("\n")
        MaskCollection.show_frame(masks(object))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "masks<-" replacement methods.
###

setGeneric("masks<-", function(x, value) standardGeneric("masks<-"))

### Setting the masks of a MaskedXString object...
setReplaceMethod("masks", signature(x="MaskedXString", value="NULL"),
    function(x, value) unmasked(x)
)
setReplaceMethod("masks", signature(x="MaskedXString", value="MaskCollection"),
    function(x, value)
    {
        if (width(value) != length(x))
            stop("the width of the mask collection must be equal ",
                 "to the length of the sequence")
        x@masks <- value
        x
    }
)

### Setting the masks of an XString object...
setReplaceMethod("masks", signature(x="XString", value="NULL"),
    function(x, value) x
)
setReplaceMethod("masks", signature(x="XString", value="ANY"),
    function(x, value)
    {
        x <- as(x, paste("Masked", xsbaseclass(x), sep=""))
        masks(x) <- value
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "Views" method.
###

setMethod("Views", "MaskedXString",
    function(subject, start=NULL, end=NULL, width=NULL, names=NULL)
    {
        if (any(active(masks(subject))))
            warning("masks were dropped")
        Views(unmasked(subject), start=start, end=end, width=width, names=names)
    }
)

