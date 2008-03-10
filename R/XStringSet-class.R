### =========================================================================
### XStringSet objects
### -------------------------------------------------------------------------
###
### The XStringSet class is a container for storing a set of XString objects
### of the same subtype (e.g. all BString instances or all DNAString instances).
### The current implementation only allows the storage of a set of strings
### that are all substrings of a common string called the super string.
### For storing XString objects that point to different XRaw objects, the user
### should use an XStringList container.
###
### Of course, the responsability of choosing between XStringSet and
### XStringList based on such an obscure criteria ("are my XString objects
### sharing the same XRaw object?") should not be left to the user.
### In the future this unfriendly situation could be worked around by renaming
### the XStringSet class -> SubstrSet and making XStringSet an interface
### only (i.e. a virtual class with no slots) that would be shared by specific
### XStringSet-like containers like XStringList and SubstrSet.
###   (1) XStringList: the current XStringList container where each sequence
###       points to its own XRaw object. Maybe some sequences are in fact
###       sharing the same XRaw object but it doesn't matter, except when they
###       all share the same, then the XStringList object can be easily and
###       very efficiently converted to a SubstrSet object.
###   (2) SubstrSet: this would remain the most efficient (i.e. compact and
###       fast) XStringSet-like container and the one still used in most use
###       cases.
### Other XStringSet-like containers could be added later e.g. the
### OnfileXStringList container: would use some delayed loading mechanism like
### what is currently in use for the BSgenome stuff. Still need to think more
### about the pros and cons of having such container though...
###
### Conversion between BStringViews, XStringSet and XStringList objects:
###   o BStringViews --> XStringSet: should always work, except when views
###     in BStringViews are out of limits.
###   o XStringSet --> BStringViews: doesn't make sense.
###   o XStringSet --> XStringList: should always work (no restriction) but
###     I can't think of any use case where doing this would be a good idea!
###     (the XStringSet container is much more efficient).
###   o XStringList --> XStringSet: would be fast and easy to implement when
###     all the sequences in XStringList share the same XRaw object (then
###     no need to copy any data). When they don't share the same XRaw object,
###     then all the sequence data in XStringList need to be copied,
###     concatenated and stuffed into a new XRaw object.
###
### Some notable differences between XStringSet and BStringViews objects:
###   - the "show" methods produce very different output
###   - the views in a BStringViews object can't have a 0-width
###   - the views in a BStringViews object can be out of limits
###

setClass("XStringSet",
    contains=".IRanges",
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
### Like for XStringList and BStringViews objects, the strict minimum of
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

.newXStringSet <- function(class, super, ranges, x, use.names)
{
    use.names <- normalize.use.names(use.names)
    if (use.names) ans_names <- names(x) else ans_names <- NULL
    new(class, super, start(ranges), width(ranges), ans_names)
}

.charToXString <- function(x, safe_locs, class)
{
    proto <- new(class, XRaw(0), 0L, 0L, check=FALSE)
    data <- .Call("STRSXP_to_XRaw",
                  x, start(safe_locs), width(safe_locs), "", enc_lkup(proto),
                  PACKAGE="Biostrings")
    new(class, data, 0L, length(data), check=FALSE)
}

.charToXStringSet <- function(x, start, end, width, use.names, baseClass)
{
    class <- paste(baseClass, "Set", sep="")
    safe_locs <- narrow(nchar(x, type="bytes"), start, end, width)
    super <- .charToXString(x, safe_locs, baseClass)
    ranges <- intToAdjacentRanges(width(safe_locs))
    .newXStringSet(class, super, ranges, x, use.names)
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
        frame <- reduce(safe_locs, TRUE)
        data <- .Call("XString_to_XRaw",
                      super(x), start(frame), width(frame), lkup,
                      PACKAGE="Biostrings")
        super <- new(baseClass, data, 0L, length(data), check=FALSE)
        ranges <- attr(frame, "inframe")
    } else {
        super <- XString(baseClass, super(x))
        ranges <- narrow(x, start, end, width)
    }
    .newXStringSet(class, super, ranges, x, use.names)
}

### Canonical conversion from BStringViews to XStringSet
BStringViewsToSet <- function(x, use.names, verbose=TRUE)
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
### Just because of those silly "AsIs" objects found in the probe packages
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
setMethod("XStringSet", "XStringSet",
    function(baseClass, x, start=NA, end=NA, width=NA, use.names=TRUE)
        .narrowXStringSet(x, start, end, width, use.names, baseClass)
)
setMethod("XStringSet", "BStringViews",
    function(baseClass, x, start=NA, end=NA, width=NA, use.names=TRUE)
        .narrowXStringSet(BStringViewsToSet(x, use.names),
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
    width <- nchar(x[[i]])
    snippetWidth <- getOption("width") - 2 - iW - widthW
    if (!is.null(names(x)))
        snippetWidth <- snippetWidth - .namesW - 1
    seq_snippet <- XString.get_snippet(x[[i]], snippetWidth)
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
        stop("attempt to modify the value of a ", sQuote(class(x)), " object")
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

