### =========================================================================
### The BStringSet class
### -------------------------------------------------------------------------
###
### The BStringSet class is a container for storing a set of BString objects
### of the same class (e.g. all BString instances or all DNAString instances).
### The current implementation only allows the storage of a set of strings
### that are all substrings of a common string called the super string.
### For storing BString objects that point to different XRaw objects, the user
### should use the BStringList container.
###
### Of course, the responsability of choosing between BStringSet and
### BStringList based on such an obscure criteria ("are my BString objects
### sharing the same XRaw object?") should not be left to the user.
### In the future this unfriendly situation could be worked around by renaming
### the BStringSet class -> SubstrSet and making BStringSet an interface only
### (i.e. a virtual class with no slots) that would be shared by specific
### BStringSet-like containers like BStringList and SubstrSet.
###   (1) BStringList: the current BStringList container where each sequence
###       points to its own XRaw object. Maybe some sequences are in fact
###       sharing the same XRaw object but it doesn't matter, except when they
###       all share the same, then the BStringList object can be easily and
###       very efficiently converted to a SubstrSet object.
###   (2) SubstrSet: this would remain the most efficient (i.e. compact and
###       fast) BStringSet-like container and the one still used in most use
###       cases.
### Other BStringSet-like containers could be added later e.g. the
### OnfileBStringList container: would use some delayed loading mechanism like
### what is currently in use for the BSgenome stuff. Still need to think more
### about the pros and cons of having such container though...
### Also worth considering: could the current BStringViews container inherit
### the BStringSet interface?
###
### Conversion between BStringViews, BStringSet and BStringList objects:
###   o BStringViews --> BStringSet: should always work, except when views in
###     BStringViews are out of limits.
###   o BStringSet --> BStringViews: doesn't make sense.
###   o BStringSet --> BStringList: should always work (no restriction) but I
###     can't think of any use case where doing this would be a good idea!
###     (the BStringSet container is much more efficient).
###   o BStringList --> BStringSet: would be fast and easy to implement when
###     all the sequences in BStringList share the same XRaw object (then
###     no need to copy any data). When they don't share the same XRaw object,
###     then all the sequence data in BStringList need to be copied,
###     concatenated and stuffed into a new XRaw object.
###
### Some notable differences between BStringSet and BStringViews objects:
###   - the "show" methods produce very different output
###   - the views in a BStringViews object can't have a 0-width
###   - the views in a BStringViews object can be out of limits
###

setClass("BStringSet",
    contains="IRanges",
    representation(
        super="BString"
    )
)

### 3 direct "BStringSet" derivations (no additional slot)
setClass("DNAStringSet", contains="BStringSet")
setClass("RNAStringSet", contains="BStringSet")
setClass("AAStringSet", contains="BStringSet")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods.
###
### Like for BStringList and BStringViews objects, the strict minimum of
### methods that must work with BStringSet objects is:
###   length, width, nchar, names
### Note that BStringSet objects inherit the "length", "width" and "names"
### methods from the IRanges class.
###

### NOT exported
setGeneric("super", function(x) standardGeneric("super"))
setMethod("super", "BStringSet", function(x) x@super)

setMethod("nchar", "BStringSet", function(x, type="chars", allowNA=FALSE) width(x))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Initialization (not intended to be used directly by the user).
###

setMethod("initialize", "BStringSet",
    function(.Object, super, start=integer(0), width=integer(0), names=NULL, check=TRUE)
    {
        .Object <- callNextMethod(.Object, start=start, width=width, names=names, check=check)
        if (check) {
            if (!is(super, "BString"))
                stop("'super' must be a BString object")
            if (length(.Object) != 0) {
                if (min(start(.Object)) < 1)
                    stop("bad spanning (some ranges are starting before the start of 'super')")
                if (max(end(.Object)) > nchar(super))
                    stop("bad spanning (some ranges are ending after the end of 'super')")
            }
        }
        slot(.Object, "super", check=FALSE) <- super
        .Object
    }
)
setMethod("initialize", "DNAStringSet",
    function(.Object, super, start=integer(0), width=integer(0), names=NULL, check=TRUE)
    {
        if (check) {
            if (!is(super, "DNAString"))
                stop("'super' must be a DNAString object")
        }
        callNextMethod(.Object, super, start=start, width=width, names=names, check=check)
    }
)
setMethod("initialize", "RNAStringSet",
    function(.Object, super, start=integer(0), width=integer(0), names=NULL, check=TRUE)
    {
        if (check) {
            if (!is(super, "RNAString"))
                stop("'super' must be a RNAString object")
        }
        callNextMethod(.Object, super, start=start, width=width, names=names, check=check)
    }
)
setMethod("initialize", "AAStringSet",
    function(.Object, super, start=integer(0), width=integer(0), names=NULL, check=TRUE)
    {
        if (check) {
            if (!is(super, "AAString"))
                stop("'super' must be a AAString object")
        }
        callNextMethod(.Object, super, start=start, width=width, names=names, check=check)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helper functions used by the versatile constructors below.
###

.newBStringSet <- function(class, super, ranges, x, use.names)
{
    use.names <- normalize.use.names(use.names)
    if (use.names) ans_names <- names(x) else ans_names <- NULL
    new(class, super, start=start(ranges),
                      width=width(ranges),
                      names=ans_names,
                      check=FALSE)
}

.charToBString <- function(x, safe_locs, class)
{
    proto <- new(class, XRaw(0), 0L, 0L, check=FALSE)
    data <- .Call("STRSXP_to_XRaw",
                  x, start(safe_locs), width(safe_locs), "", enc_lkup(proto),
                  PACKAGE="Biostrings")
    new(class, data, 0L, length(data), check=FALSE)
}

.charToBStringSet <- function(x, start, end, width, use.names, baseClass, check)
{
    class <- paste(baseClass, "Set", sep="")
    safe_locs <- narrow(nchar(x, type="bytes"), start, end, width)
    super <- .charToBString(x, safe_locs, baseClass)
    ranges <- intToAdjacentRanges(width(safe_locs))
    .newBStringSet(class, super, ranges, x, use.names)
}

.narrowBStringSet <- function(x, start, end, width, use.names, baseClass, check)
{
    class <- paste(baseClass, "Set", sep="")
    lkup <- getBStringTypeConversionLookup(class(super(x)), baseClass)
    if (!is.null(lkup)) {
        ## The frame is the strict minimal region of the original data that
        ## needs to be copied.
        ## This will be paticularly useful (and will significantly reduce the
        ## memory footprint) when the sequences in 'x' point to regions in
        ## 'super(x)' that have a lot of overlapping.
        safe_locs <- narrow(x, start, end, width)
        frame <- reduce(safe_locs, TRUE)
        data <- .Call("BString_to_XRaw",
                      super(x), start(frame), width(frame), lkup,
                      PACKAGE="Biostrings")
        super <- new(baseClass, data, 0L, length(data), check=FALSE)
        ranges <- attr(frame, "inframe")
    } else {
        super <- mkBString(baseClass, super(x))
        ranges <- narrow(x, start, end, width)
    }
    .newBStringSet(class, super, ranges, x, use.names)
}

### Canonical conversion from BStringViews to a BStringSet (or derived)
BStringViewsToSet <- function(x, use.names, verbose=TRUE)
{
    ranges <- restrict(as(x, "IRanges"), 1L, nchar(subject(x)),
                       keep.all.ranges=TRUE,
                       use.names=use.names)
    if (verbose && any(width(ranges) < width(x)))
        warning("trimming \"out of limits\" views")
    class <- paste(class(subject(x)), "Set", sep="")
    new(class, subject(x), start=start(ranges), width=width(ranges),
               names=names(ranges), check=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The user-friendly versatile constructors.
###
### All these constructors use the uSEW (user-specified Start/End/Width)
### interface.
###

setGeneric("BStringSet", signature="x",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
        standardGeneric("BStringSet")
)
setGeneric("DNAStringSet", signature="x",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
        standardGeneric("DNAStringSet")
)
setGeneric("RNAStringSet", signature="x",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
        standardGeneric("RNAStringSet")
)
setGeneric("AAStringSet", signature="x",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
        standardGeneric("AAStringSet")
)

setMethod("BStringSet", "character",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
        .charToBStringSet(x, start, end, width, use.names, "BString", check)
)
setMethod("DNAStringSet", "character",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
        .charToBStringSet(x, start, end, width, use.names, "DNAString", check)
)
setMethod("RNAStringSet", "character",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
        .charToBStringSet(x, start, end, width, use.names, "RNAString", check)
)
setMethod("AAStringSet", "character",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
        .charToBStringSet(x, start, end, width, use.names, "AAString", check)
)

### Just because of those silly "AsIs" objects found in the probe packages
### (e.g. drosophila2probe$sequence)
setMethod("BStringSet", "AsIs",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
    {
        if (!is.character(x))
            stop("unsuported input type")
        class(x) <- "character" # keeps the names (unlike as.character())
	BStringSet(x, start=start, end=end, width=width,
                   use.names=use.names, check=check)
    }
)
setMethod("DNAStringSet", "AsIs",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
    {
        if (!is.character(x))
            stop("unsuported input type")
        class(x) <- "character" # keeps the names (unlike as.character())
	DNAStringSet(x, start=start, end=end, width=width,
                     use.names=use.names, check=check)
    }
)
setMethod("RNAStringSet", "AsIs",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
    {
        if (!is.character(x))
            stop("unsuported input type")
        class(x) <- "character" # keeps the names (unlike as.character())
	RNAStringSet(x, start=start, end=end, width=width,
                     use.names=use.names, check=check)
    }
)
setMethod("AAStringSet", "AsIs",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
    {
        if (!is.character(x))
            stop("unsupported input type")
        class(x) <- "character" # keeps the names (unlike as.character())
	AAStringSet(x, start=start, end=end, width=width,
                    use.names=use.names, check=check)
    }
)

setMethod("BStringSet", "BStringSet",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
        .narrowBStringSet(x, start, end, width, use.names, "BString", check)
)
setMethod("DNAStringSet", "BStringSet",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
        .narrowBStringSet(x, start, end, width, use.names, "DNAString", check)
)
setMethod("RNAStringSet", "BStringSet",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
        .narrowBStringSet(x, start, end, width, use.names, "RNAString", check)
)
setMethod("AAStringSet", "BStringSet",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
        .narrowBStringSet(x, start, end, width, use.names, "AAString", check)
)

setMethod("BStringSet", "BStringViews",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
        .narrowBStringSet(BStringViewsToSet(x, use.names),
                          start, end, width, TRUE, "BString", check)
)
setMethod("DNAStringSet", "BStringViews",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
        .narrowBStringSet(BStringViewsToSet(x, use.names),
                          start, end, width, TRUE, "DNAString", check)
)
setMethod("RNAStringSet", "BStringViews",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
        .narrowBStringSet(BStringViewsToSet(x, use.names),
                          start, end, width, TRUE, "RNAString", check)
)
setMethod("AAStringSet", "BStringViews",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
        .narrowBStringSet(BStringViewsToSet(x, use.names),
                          start, end, width, TRUE, "AAString", check)
)

### Not exported
mkBStringSet <- function(baseClass, x)
{
    class <- paste(baseClass, "Set", sep="")
    do.call(class, list(x=x))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method.
###

.namesW <- 20

BStringSet.show_frame_header <- function(iW, widthW, with.names)
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

BStringSet.show_frame_line <- function(x, i, iW, widthW)
{
    width <- nchar(x[[i]])
    snippetWidth <- getOption("width") - 2 - iW - widthW
    if (!is.null(names(x)))
        snippetWidth <- snippetWidth - .namesW - 1
    seq_snippet <- BString.get_snippet(x[[i]], snippetWidth)
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
BStringSet.show_frame <- function(x, half_nrow=9L)
{
    lx <- length(x)
    iW <- nchar(as.character(lx)) + 2 # 2 for the brackets
    ncharMax <- max(nchar(x))
    widthW <- max(nchar(ncharMax), nchar("width"))
    BStringSet.show_frame_header(iW, widthW, !is.null(names(x)))
    if (lx <= 2*half_nrow+1) {
        for (i in seq_len(lx))
            BStringSet.show_frame_line(x, i, iW, widthW)
    } else {
        for (i in 1:half_nrow)
            BStringSet.show_frame_line(x, i, iW, widthW)
        cat(format("...", width=iW, justify="right"),
            format("...", width=widthW, justify="right"),
            "...\n")
        for (i in (lx-half_nrow+1L):lx)
            BStringSet.show_frame_line(x, i, iW, widthW)
    }
}

setMethod("show", "BStringSet",
    function(object)
    {
        cat("  A ", class(object), " instance of length ", length(object), "\n", sep="")
        if (length(object) != 0)
            BStringSet.show_frame(object)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

### Extract the i-th element of a BStringSet object as a BString object.
### Return a BString object of the same class as super(x).
### Example:
###   bs <- BString("ABCD-1234-abcd")
###   bset <- new("BStringSet", bs, 1:8, 2L*(7:0))
###   bset[[3]]
### Supported 'i' types: numeric vector of length 1.
setMethod("[[", "BStringSet",
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
        BString.substr(super(x), start, end)
    }
)

setReplaceMethod("[[", "BStringSet",
    function(x, i, j,..., value)
    {
        stop("attempt to modify the value of a ", sQuote(class(x)), " object")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Equality.
###

### WON'T START THIS UNLESS SOMEONE HAS A USE CASE...
### Look at BStringViews-class.R for how to do this.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Other methods.
###

setMethod("as.list", "BStringSet",
    function(x, ...)
    {
        ans <- lapply(seq_len(length(x)), function(i) x[[i]])
        names(ans) <- names(x)
        ans
    }
)

setMethod("as.character", "BStringSet",
    function(x, use.names=TRUE)
    {
        use.names <- normalize.use.names(use.names)
        ans <- sapply(seq_len(length(x)), function(i) as.character(x[[i]]))
        if (use.names)
            names(ans) <- names(x)
        ans
    }
)

