### =========================================================================
### XStringSet objects
### -------------------------------------------------------------------------
###
### The XStringSet class is a container for storing a set of XString objects
### of the same subtype (e.g. all elements are BString objects or they are
### all DNAString objects).
###
### The current implementation only allows for storage of a set of strings
### that belong to the same XRaw object i.e. all the elements of an XStringSet
### object must be substrings of a common string called the "super string".
### So for storing XString objects that point to different XRaw objects, the
### user must use an XStringList container.
### There are 3 problems with this:
###   1. It's not user-friendly. The responsability of choosing between
###      XStringSet and XStringList based on such an obscure criteria ("are
###      my XString objects sharing the same XRaw object?") should not be left
###      to the user.
###   2. It's currently not possible to add new elements to an existing
###      XStringSet object. Well, in fact it can be done, but only if the new
###      elements belong to the XRaw object shared by the existing elements.
###      Such restriction would not make sense from a user point of view.
###   3. The XStringList container is not as efficient as the XStringSet
###      container.
### This could be changed (and maybe this will show up in the 2.9 series) by
### using something like this for the XStringSet class:
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
### This new XStringSet container would combine the efficiency of the old one
### and the flexibility of the XStringList container (which can then be
### removed).
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
### Accessor methods.
###
### Like for XStringList and XStringViews objects, the strict minimum of
### methods that must work with XStringSet objects is:
###   length, width, nchar, names
### Note that XStringSet objects inherit the "length", "width" and "names"
### methods from the .IRanges class.
###

### NOT exported
setGeneric("super", function(x) standardGeneric("super"))
setMethod("super", "XStringSet", function(x) x@super)

setMethod("nchar", "XStringSet",
    function(x, type="chars", allowNA=FALSE) width(x)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "alphabet" method.
###

setMethod("alphabet", "XStringSet", function(x) alphabet(super(x)))


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
### Helper functions used by the versatile constructors below.
###

newXStringSet <- function(class, super, ranges, use.names=FALSE, names=NULL)
{
    if (!normalize.use.names(use.names))
        names <- NULL
    new(class, super, start(ranges), width(ranges), names)
}

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
    ranges <- new("IRanges", start=1L, width=length(super), check=FALSE)
    newXStringSet(class, super, ranges)
}

.narrowXStringSet <- function(x, start, end, width, use.names, baseClass)
{
    class <- paste(baseClass, "Set", sep="")
    lkup <- getXStringSubtypeConversionLookup(class(super(x)), baseClass)
    if (!is.null(lkup)) {
        ## The frame is the strict minimal region of the original data that
        ## needs to be copied.
        ## This will be paticularly useful (and will significantly reduce the
        ## memory footprint) when the sequences in 'x' point to regions in
        ## 'super(x)' that have a lot of overlapping.
        safe_locs <- narrow(x, start, end, width)
        frame <- reduce(safe_locs, with.inframe.attrib=TRUE)
        xdata <- .Call("new_XRaw_from_XString",
                       super(x), start(frame), width(frame), lkup,
                       PACKAGE="Biostrings")
        super <- new(baseClass, xdata=xdata, length=length(xdata))
        ranges <- attr(frame, "inframe")
    } else {
        super <- XString(baseClass, super(x))
        ranges <- narrow(x, start, end, width)
    }
    newXStringSet(class, super, ranges, use.names=use.names, names=names(x))
}

### Canonical conversion from XStringViews to XStringSet
XStringViewsToSet <- function(x, use.names, verbose=TRUE)
{
    ranges <- restrict(as(x, "IRanges"), 1L, nchar(subject(x)),
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
        .charToXStringSet(x, start, end, width, use.names, baseClass)
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
        .narrowXStringSet(x, start, end, width, use.names, baseClass)
)
setMethod("XStringSet", "XStringViews",
    function(baseClass, x, start=NA, end=NA, width=NA, use.names=TRUE)
        .narrowXStringSet(XStringViewsToSet(x, use.names),
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
        use.names <- normalize.use.names(use.names)
        ans <- sapply(seq_len(length(x)), function(i) as.character(x[[i]]))
        if (use.names)
            names(ans) <- names(x)
        ans
    }
)

