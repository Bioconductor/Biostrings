### =========================================================================
### XStringSet objects
### -------------------------------------------------------------------------
###
### The XStringSet class is a container for storing a set of XString objects
### of the same subtype (e.g. all elements are BString objects or they are
### all DNAString objects etc).
###
### The current implementation only allows for storage of a set of strings
### that belong to the same XRaw object i.e. all the elements of an XStringSet
### object must be substrings of a common string called the "super string".
### The old XStringList container (Biostrings 2.8) didn't have this limitation:
### it could hold XString objects that pointed to different XRaw objects but
### it was so slow that I decided to replace it by the much more efficient
### XStringSet container.
### Maybe the best of both world, or at least a good enough trade-off, could be
### obtained by defining the XStringSet class like this:
###
###   setClass("XRawViews",
###     contains="LockedIRanges",
###     representation(
###         subject="XRaw"
###     )
###   )
###
###   setClass("XStringSet",
###     representation(
###         "VIRTUAL",
###         xrvlist="list",   # a list of XRawViews objects
###         strong="integer",
###         weak="integer"
###     )
###   )
###
### x@strong and x@weak (both of length the number of elements) are maps used
### to "find" the i-th element in x i.e. x[[i]] can be efficiently extracted
### with:
###     xrv <- x@xrvlist[[x@strong[i]]]
###     start <- start(xrv)[x@weak[i]]
###     width <- width(xrv)[x@weak[i]]
###     new(baseClass, xdata=xrv@subject, offset=start-1L, length=width)
### This new XStringSet container would still be almost as efficient as the old
### one (at least for the restricted set of use cases that the old one was
### supporting i.e. when x@xrvlist holds a list of length 1) and the
### flexibility of the old XStringList container.
###
### Some notable differences between XStringSet and XStringViews objects:
###   - the "show" methods produce different output
###   - the views in an XStringViews object can't have a 0-width
###   - the views in an XStringViews object can be out of limits
###

setClass("XStringSet",
    contains="LockedIRanges",
    representation(
        "VIRTUAL",
        super="XString"
    )
)

setClass("BStringSet",
    contains="XStringSet",
    representation(
        super="BString"
    )
)
setClass("DNAStringSet",
    contains="XStringSet",
    representation(
        super="DNAString"
    )
)
setClass("RNAStringSet",
    contains="XStringSet",
    representation(
        super="RNAString"
    )
)
setClass("AAStringSet",
    contains="XStringSet",
    representation(
        super="AAString"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "super" accessor method (NOT exported).
###

setGeneric("super", function(x) standardGeneric("super"))
setMethod("super", "XStringSet", function(x) x@super)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The core XString API.
###
### The core XString API is the strict minimal set of methods that must work
### for XString, XStringSet, XStringViews and MaskedXString objects.
### It currently consists of the following methods:
###   o NOT exported: baseXStringSubtype, codes, codec, enc_lkup, dec_lkup
###   o exported: alphabet, length, nchar
###

### NOT exported
setMethod("baseXStringSubtype", "XStringSet",
    function(x) baseXStringSubtype(super(x))
)
setMethod("codes", "XStringSet", function(x, ...) codes(super(x), ...))
setMethod("codec", "XStringSet", function(x) codec(super(x)))
setMethod("enc_lkup", "XStringSet", function(x) enc_lkup(super(x)))
setMethod("dec_lkup", "XStringSet", function(x) dec_lkup(super(x)))

### exported
setMethod("alphabet", "XStringSet", function(x) alphabet(super(x)))

setMethod("nchar", "XStringSet",
    function(x, type="chars", allowNA=FALSE) width(x)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Initialization (not intended to be used directly by the user).
###

setMethod("initialize", "XStringSet",
    function(.Object, super, start, width, names)
    {
        .Object <- callNextMethod(.Object, start, width, names)
        .Object@super <- super
        .Object
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helper functions (NOT exported) used by the coercion methods and
### versatile (and user friendly) constructors below.
###

newXStringSet <- function(class, super, ranges, use.names=FALSE, names=NULL)
{
    if (!normargUseNames(use.names))
        names <- NULL
    new(class, super, start(ranges), width(ranges), names)
}

### This is an endomorphism iff 'baseClass' is NULL, otherwise it is NOT!
compactXStringSet <- function(x, baseClass=NULL)
{
    from_baseClass <- baseXStringSubtype(x)
    if (is.null(baseClass))
        to_baseClass <- from_baseClass
    else
        to_baseClass <- baseClass
    lkup <- getXStringSubtypeConversionLookup(from_baseClass, to_baseClass)
    ## The frame is the strict minimal set of regions in 'super(x)' that need
    ## to be copied. Hence compacting 'x' returns an XStringSet object 'y'
    ## where 'length(super(y))' can be significantly smaller than
    ## 'length(super(x))' especially if the elements in 'x' cover a small part
    ## of 'super(x)'.
    frame <- reduce(x, with.inframe.attrib=TRUE)
    xdata <- .Call("new_XRaw_from_XString",
                   super(x), start(frame), width(frame), lkup,
                   PACKAGE="Biostrings")
    ans_super <- new(to_baseClass, xdata=xdata, length=length(xdata))
    ans_ranges <- attr(frame, "inframe")
    if (is.null(baseClass)) {
        ## Endomorphism
        x@super <- ans_super
        IRanges:::unsafe.update(x, start=start(ans_ranges), width=width(ans_ranges))
    } else {
        ## NOT an endomorphism
        ans_class <- paste(to_baseClass, "Set", sep="")
        newXStringSet(ans_class, ans_super, ans_ranges, use.names=TRUE, names=names(x))
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

.as.from_XStringSet_to_XStringSet <- function(x, baseClass)
{
    from_baseClass <- baseXStringSubtype(x)
    lkup <- getXStringSubtypeConversionLookup(from_baseClass, baseClass)
    if (!is.null(lkup))
        return(compactXStringSet(x, baseClass=baseClass))
    ans_class <- paste(baseClass, "Set", sep="")
    if (is(x, ans_class))
        ans_super <- super(x)
    else
        ans_super <- XString(baseClass, super(x))
    new(ans_class, ans_super, start(x), width(x), names(x))
}

setAs("XStringSet", "BStringSet",
    function(from) .as.from_XStringSet_to_XStringSet(from, "BString")
)
setAs("XStringSet", "DNAStringSet",
    function(from) .as.from_XStringSet_to_XStringSet(from, "DNAString")
)
setAs("XStringSet", "RNAStringSet",
    function(from) .as.from_XStringSet_to_XStringSet(from, "RNAString")
)
setAs("XStringSet", "AAStringSet",
    function(from) .as.from_XStringSet_to_XStringSet(from, "AAString")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### More helper functions used by the versatile (and user friendly)
### constructors below.
###

.charToXString <- function(x, safe_locs, class)
{
    proto <- newEmptyXString(class)
    xdata <- .Call("new_XRaw_from_STRSXP",
                   x, start(safe_locs), width(safe_locs), "", enc_lkup(proto),
                   PACKAGE="Biostrings")
    new(class, xdata=xdata, length=length(xdata))
}

.charToXStringSet <- function(x, start, end, width, use.names, baseClass)
{
    class <- paste(baseClass, "Set", sep="")
    safe_locs <- narrow(nchar(x, type="bytes"), start, end, width)
    super <- .charToXString(x, safe_locs, baseClass)
    ranges <- intToAdjacentRanges(width(safe_locs))
    newXStringSet(class, super, ranges, use.names=use.names, names=names(x))
}

.XStringToXStringSet <- function(x, start, end, width, use.names, baseClass)
{
    class <- paste(baseClass, "Set", sep="")
    super <- subseq(x, start=start, end=end, width=width)
    ranges <- new("LockedIRanges", start=1L, width=length(super), check=FALSE)
    newXStringSet(class, super, ranges)
}

.narrowAndCoerceXStringSet <- function(x, start, end, width, use.names, baseClass)
{
    y <- narrow(x, start=start, end=end, width=width, use.names=use.names)
    .as.from_XStringSet_to_XStringSet(y, baseClass)
}

### Canonical conversion from XStringViews to XStringSet
XStringViewsToSet <- function(x, use.names, verbose=TRUE)
{
    ranges <- restrict(as(x, "IRanges"), start=1L, end=nchar(subject(x)),
                       keep.all.ranges=TRUE,
                       use.names=use.names)
    if (verbose && any(width(ranges) < width(x)))
        warning("trimming \"out of limits\" views")
    class <- paste(class(subject(x)), "Set", sep="")
    new(class, subject(x), start(ranges), width(ranges), names(ranges))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The user-friendly versatile constructors.
###
### All these constructors use the uSEW (user-specified Start/End/Width)
### interface.
###

setGeneric("XStringSet", signature="x",
    function(baseClass, x, start=NA, end=NA, width=NA, use.names=TRUE)
        standardGeneric("XStringSet")
)
setMethod("XStringSet", "character",
    function(baseClass, x, start=NA, end=NA, width=NA, use.names=TRUE)
    {
        if (is.null(baseClass))
            baseClass <- "BString"
        .charToXStringSet(x, start, end, width, use.names, baseClass)
    }
)
### Just because of the silly "AsIs" objects found in the probe packages
### (e.g. drosophila2probe$sequence)
setMethod("XStringSet", "AsIs",
    function(baseClass, x, start=NA, end=NA, width=NA, use.names=TRUE)
    {
        if (!is.character(x))
            stop("unsuported input type")
        class(x) <- "character" # keeps the names (unlike as.character())
	.charToXStringSet(x, start, end, width, use.names, baseClass)
    }
)
setMethod("XStringSet", "XString",
    function(baseClass, x, start=NA, end=NA, width=NA, use.names=TRUE)
        .XStringToXStringSet(x, start, end, width, use.names, baseClass)
)
setMethod("XStringSet", "XStringSet",
    function(baseClass, x, start=NA, end=NA, width=NA, use.names=TRUE)
        .narrowAndCoerceXStringSet(x, start, end, width, use.names, baseClass)
)
setMethod("XStringSet", "XStringViews",
    function(baseClass, x, start=NA, end=NA, width=NA, use.names=TRUE)
        .narrowAndCoerceXStringSet(XStringViewsToSet(x, use.names),
                                   start, end, width, TRUE, baseClass)
)

BStringSet <- function(x, start=NA, end=NA, width=NA, use.names=TRUE)
    XStringSet("BString", x, start=start, end=end, width=width,
                             use.names=use.names)
DNAStringSet <- function(x, start=NA, end=NA, width=NA, use.names=TRUE)
    XStringSet("DNAString", x, start=start, end=end, width=width,
                               use.names=use.names)
RNAStringSet <- function(x, start=NA, end=NA, width=NA, use.names=TRUE)
    XStringSet("RNAString", x, start=start, end=end, width=width,
                               use.names=use.names)
AAStringSet <- function(x, start=NA, end=NA, width=NA, use.names=TRUE)
    XStringSet("AAString", x, start=start, end=end, width=width,
                              use.names=use.names)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method.
###

.namesW <- 20

.XStringSet.show_frame_header <- function(iW, widthW, with.names)
{
    cat(format("", width=iW+1),
        format("width", width=widthW, justify="right"),
        sep="")
    if (with.names) {
        cat(format(" seq", width=getOption("width")-iW-widthW-.namesW-1),
            format("names", width=.namesW, justify="left"),
            sep="")
    } else {
        cat(" seq")
    }
    cat("\n")
}

.XStringSet.show_frame_line <- function(x, i, iW, widthW)
{
    width <- nchar(x)[i]
    snippetWidth <- getOption("width") - 2 - iW - widthW
    if (!is.null(names(x)))
        snippetWidth <- snippetWidth - .namesW - 1
    seq_snippet <- toSeqSnippet(x[[i]], snippetWidth)
    if (!is.null(names(x)))
        seq_snippet <- format(seq_snippet, width=snippetWidth)
    cat(format(paste("[", i,"]", sep=""), width=iW, justify="right"), " ",
        format(width, width=widthW, justify="right"), " ",
        seq_snippet,
        sep="")
    if (!is.null(names(x))) {
        snippet_name <- names(x)[i]
        if (is.na(snippet_name))
            snippet_name <- "<NA>"
        else if (nchar(snippet_name) > .namesW)
            snippet_name <- paste(substr(snippet_name, 1, .namesW-3), "...", sep="")
        cat(" ", snippet_name, sep="")
    }
    cat("\n")
}

### 'half_nrow' must be >= 1
.XStringSet.show_frame <- function(x, half_nrow=9L)
{
    lx <- length(x)
    iW <- nchar(as.character(lx)) + 2 # 2 for the brackets
    ncharMax <- max(nchar(x))
    widthW <- max(nchar(ncharMax), nchar("width"))
    .XStringSet.show_frame_header(iW, widthW, !is.null(names(x)))
    if (lx <= 2*half_nrow+1) {
        for (i in seq_len(lx))
            .XStringSet.show_frame_line(x, i, iW, widthW)
    } else {
        for (i in 1:half_nrow)
            .XStringSet.show_frame_line(x, i, iW, widthW)
        cat(format("...", width=iW, justify="right"),
            format("...", width=widthW, justify="right"),
            "...\n")
        for (i in (lx-half_nrow+1L):lx)
            .XStringSet.show_frame_line(x, i, iW, widthW)
    }
}

setMethod("show", "XStringSet",
    function(object)
    {
        cat("  A ", class(object), " instance of length ", length(object), "\n", sep="")
        if (length(object) != 0)
            .XStringSet.show_frame(object)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

### Extract the i-th element of an XStringSet object as an XString object.
### Return an XString object of the same subtype as super(x).
### Example:
###   bs <- BString("ABCD-1234-abcd")
###   bset <- new("BStringSet", bs, 1:8, 2L*(7:0))
###   bset[[3]]
### Supported 'i' types: numeric vector of length 1.
setMethod("[[", "XStringSet",
    function(x, i, j, ...)
    {
        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")
        if (missing(i))
            stop("subscript is missing")
        if (is.character(i))
            stop("cannot subset a ", class(x), " object by names")
        if (!is.numeric(i))
            stop("invalid subscript type")
        if (length(i) < 1L)
            stop("attempt to select less than one element")
        if (length(i) > 1L)
            stop("attempt to select more than one element")
        if (is.na(i))
            stop("subscript cannot be NA")
        if (i < 1L || i > length(x))
            stop("subscript out of bounds")
        start <- start(x)[i]
        end <- end(x)[i]
        XString.substr(super(x), start, end)
    }
)

setReplaceMethod("[[", "XStringSet",
    function(x, i, j,..., value)
    {
        stop("attempt to modify the value of a ", class(x), " instance")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "append" method.
###
### TODO: The current version performs too many copies of the sequence data!
###       This will change with the redesign of the XStringSet class (see
###       long comment at the beginning of this file).
###

setMethod("append", "XStringSet",
    function(x, values, after=length(x))
    {
        if (!is(values, "XStringSet"))
            stop("'values' must be an XStringSet object")
        baseClass <- baseXStringSubtype(x)
        if (baseXStringSubtype(values) != baseClass)
            stop("'x' and 'values' must be XStringSet objects of the same subtype")
        if (!isSingleNumber(after))
            stop("'after' must be a single number")
        if (after != length(x))
            stop("'after != length(x)' is not supported for XStringSet objects, sorry!")
        ans_class <- paste(baseClass, "Set", sep="")
        cx <- compactXStringSet(x)
        cvalues <- compactXStringSet(values, baseClass=baseClass)
        ans_super <- XString.append(super(cx), super(cvalues))
        ans_start <- c(start(cx), start(cvalues) + length(super(cx)))
        ans_width <- c(width(cx), width(cvalues))
        ans_names <- c(names(cx), names(cvalues))
        new(ans_class, ans_super, ans_start, ans_width, ans_names)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Equality.
###

### WON'T START THIS UNLESS SOMEONE HAS A USE CASE...
### Look at XStringViews-class.R for how to do this.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Other methods.
###

setMethod("as.list", "XStringSet",
    function(x, ...)
    {
        ans <- lapply(seq_len(length(x)), function(i) x[[i]])
        names(ans) <- names(x)
        ans
    }
)

setMethod("as.character", "XStringSet",
    function(x, use.names=TRUE)
    {
        use.names <- normargUseNames(use.names)
        ans <- .Call("XStringSet_as_STRSXP",
                     x, dec_lkup(x),
                     PACKAGE="Biostrings")
        if (use.names)
            names(ans) <- names(x)
        ans
    }
)

