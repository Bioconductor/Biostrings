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
    contains="SeqLocs",
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

setGeneric("super", function(x) standardGeneric("super"))
setMethod("super", "BStringSet", function(x) x@super)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Initialization (not intended to be used directly by the user).
###

setMethod("initialize", "BStringSet",
    function(.Object, super, start=integer(0), nchar=integer(0), names=NULL, check=TRUE)
    {
        .Object <- callNextMethod(.Object, start=start, nchar=nchar, names=names, check=check)
        if (check) {
            if (!is(super, "BString"))
                stop("'super' must be a BString object")
            if (length(.Object) != 0 && max(end(.Object)) > nchar(super))
                stop("some start/nchar locations are ending after the end of 'super'")
        }
        slot(.Object, "super", check=FALSE) <- super
        .Object
    }
)
setMethod("initialize", "DNAStringSet",
    function(.Object, super, start=integer(0), nchar=integer(0), names=NULL, check=TRUE)
    {
        if (check) {
            if (!is(super, "DNAString"))
                stop("'super' must be a DNAString object")
        }
        callNextMethod(.Object, super, start, nchar, names, check)
    }
)
setMethod("initialize", "RNAStringSet",
    function(.Object, super, start=integer(0), nchar=integer(0), names=NULL, check=TRUE)
    {
        if (check) {
            if (!is(super, "RNAString"))
                stop("'super' must be a RNAString object")
        }
        callNextMethod(.Object, super, start, nchar, names, check)
    }
)
setMethod("initialize", "AAStringSet",
    function(.Object, super, start=integer(0), nchar=integer(0), names=NULL, check=TRUE)
    {
        if (check) {
            if (!is(super, "AAString"))
                stop("'super' must be a AAString object")
        }
        callNextMethod(.Object, super, start, nchar, names, check)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helper functions used by the versatile constructors below.
###

.charseqsToBString <- function(x, start, nchar, class)
{
    proto <- new(class, XRaw(0), 0L, 0L, check=FALSE)
    data <- .Call("STRSXP_to_XRaw",
                  x, start, nchar, enc_lkup(proto),
                  PACKAGE="Biostrings")
    new(class, data, 0L, length(data), check=FALSE)
}

.BStringSetToBString <- function(x, start, nchar, class)
{
    proto <- new(class, XRaw(0), 0L, 0L, check=FALSE)
    data <- .Call("BStringSet_to_XRaw",
                  x, start, nchar, enc_lkup(proto),
                  PACKAGE="Biostrings")
    new(class, data, 0L, length(data), check=FALSE)
}

.charseqsToBStringSet <- function(x, start, end, nchar, baseClass, check)
{
    locs <- SEN2locs(start, end, nchar, nchar(x, type="bytes"), check=check)
    class <- paste(baseClass, "Set", sep="")
    super <- .charseqsToBString(x, locs$start, locs$nchar, baseClass)
    new(class, super, start=getStartForAdjacentSeqs(locs$nchar),
                      nchar=locs$nchar,
                      names=names(x),
                      check=FALSE)
}

.narrowBStringSet <- function(x, start, end, nchar, baseClass, check)
{
    locs <- SEN2locs(start, end, nchar, nchar(x), check=check)
    class <- paste(baseClass, "Set", sep="")
    if (copyDataOnBStringTypeConversion(class(super(x)), baseClass)) {
        super <- .BStringSetToBString(x, locs$start, locs$nchar, baseClass)
        ans_start <- getStartForAdjacentSeqs(locs$nchar)
    } else {
        super <- mkBString(baseClass, super(x))
        ans_start <- start(x)+locs$start-1L
    }
    new(class, super, start=ans_start,
                      nchar=locs$nchar,
                      names=names(x),
                      check=FALSE)
}

.BStringViewsToBStringSet <- function(x, start, end, nchar, baseClass, check)
{
    locs <- SEN2locs(start, end, nchar, width(x), check=check)
    class <- paste(baseClass, "Set", sep="")
    if (copyDataOnBStringTypeConversion(class(subject(x)), baseClass)) {
        stop("converting a BStringViews object with a ", class(subject(x)), " subject\n",
             "  into a ", class, " object is not supported yet, sorry!\n",
             "  For now you need to convert it into a ",
             paste(class(subject(x)), "Set", sep=""), " object first")
    } else {
        super <- mkBString(baseClass, subject(x))
        ans_start <- start(x)+locs$start-1L
    }
    new(class, super, start=ans_start,
                      nchar=locs$nchar,
                      names=desc(x),
                      check=TRUE) # TRUE for catching out of limits views
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The user-friendly versatile constructors.
###
### All these constructors use the SEN (Start/End/Nchar) interface.
###

setGeneric("BStringSet", signature="x",
    function(x, start=NA, end=NA, nchar=NA, check=TRUE) standardGeneric("BStringSet"))
setGeneric("DNAStringSet", signature="x",
    function(x, start=NA, end=NA, nchar=NA, check=TRUE) standardGeneric("DNAStringSet"))
setGeneric("RNAStringSet", signature="x",
    function(x, start=NA, end=NA, nchar=NA, check=TRUE) standardGeneric("RNAStringSet"))
setGeneric("AAStringSet", signature="x",
    function(x, start=NA, end=NA, nchar=NA, check=TRUE) standardGeneric("AAStringSet"))

setMethod("BStringSet", "character",
    function(x, start=NA, end=NA, nchar=NA, check=TRUE)
        .charseqsToBStringSet(x, start, end, nchar, "BString", check)
)
setMethod("DNAStringSet", "character",
    function(x, start=NA, end=NA, nchar=NA, check=TRUE)
        .charseqsToBStringSet(x, start, end, nchar, "DNAString", check)
)
setMethod("RNAStringSet", "character",
    function(x, start=NA, end=NA, nchar=NA, check=TRUE)
        .charseqsToBStringSet(x, start, end, nchar, "RNAString", check)
)
setMethod("AAStringSet", "character",
    function(x, start=NA, end=NA, nchar=NA, check=TRUE)
        .charseqsToBStringSet(x, start, end, nchar, "AAString", check)
)

### Just because of those silly "AsIs" objects found in the probe packages
### (e.g. drosophila2probe$sequence)
setMethod("BStringSet", "AsIs",
    function(x, start=NA, end=NA, nchar=NA, check=TRUE)
    {
        if (!is.character(x))
            stop("unsuported input type")
        class(x) <- "character" # keeps the names (unlike as.character())
	BStringSet(x, start=start, end=end, nchar=nchar, check=check)
    }
)
setMethod("DNAStringSet", "AsIs",
    function(x, start=NA, end=NA, nchar=NA, check=TRUE)
    {
        if (!is.character(x))
            stop("unsuported input type")
        class(x) <- "character" # keeps the names (unlike as.character())
	DNAStringSet(x, start=start, end=end, nchar=nchar, check=check)
    }
)
setMethod("RNAStringSet", "AsIs",
    function(x, start=NA, end=NA, nchar=NA, check=TRUE)
    {
        if (!is.character(x))
            stop("unsuported input type")
        class(x) <- "character" # keeps the names (unlike as.character())
	RNAStringSet(x, start=start, end=end, nchar=nchar, check=check)
    }
)
setMethod("AAStringSet", "AsIs",
    function(x, start=NA, end=NA, nchar=NA, check=TRUE)
    {
        if (!is.character(x))
            stop("unsupported input type")
        class(x) <- "character" # keeps the names (unlike as.character())
	AAStringSet(x, start=start, end=end, nchar=nchar, check=check)
    }
)

setMethod("BStringSet", "BStringSet",
    function(x, start=NA, end=NA, nchar=NA, check=TRUE)
        .narrowBStringSet(x, start, end, nchar, "BString", check)
)
setMethod("DNAStringSet", "BStringSet",
    function(x, start=NA, end=NA, nchar=NA, check=TRUE)
        .narrowBStringSet(x, start, end, nchar, "DNAString", check)
)
setMethod("RNAStringSet", "BStringSet",
    function(x, start=NA, end=NA, nchar=NA, check=TRUE)
        .narrowBStringSet(x, start, end, nchar, "RNAString", check)
)
setMethod("AAStringSet", "BStringSet",
    function(x, start=NA, end=NA, nchar=NA, check=TRUE)
        .narrowBStringSet(x, start, end, nchar, "AAString", check)
)

setMethod("BStringSet", "BStringViews",
    function(x, start=NA, end=NA, nchar=NA, check=TRUE)
        .BStringViewsToBStringSet(x, start, end, nchar, "BString", check)
)
setMethod("DNAStringSet", "BStringViews",
    function(x, start=NA, end=NA, nchar=NA, check=TRUE)
        .BStringViewsToBStringSet(x, start, end, nchar, "DNAString", check)
)
setMethod("RNAStringSet", "BStringViews",
    function(x, start=NA, end=NA, nchar=NA, check=TRUE)
        .BStringViewsToBStringSet(x, start, end, nchar, "RNAString", check)
)
setMethod("AAStringSet", "BStringViews",
    function(x, start=NA, end=NA, nchar=NA, check=TRUE)
        .BStringViewsToBStringSet(x, start, end, nchar, "AAString", check)
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

BStringSet.show_frame_header <- function(iW, ncharW, with.names)
{
    cat(format("", width=iW+1),
        format("nchar", width=ncharW, justify="right"),
        sep="")
    if (with.names) {
        cat(format(" seq", width=getOption("width")-iW-ncharW-.namesW-1),
            format("names", width=.namesW, justify="left"),
            sep="")
    } else {
        cat(" seq")
    }
    cat("\n")
}

BStringSet.show_frame_line <- function(x, i, iW, ncharW)
{
    nchar <- nchar(x[[i]])
    snippetWidth <- getOption("width") - 2 - iW - ncharW
    if (!is.null(names(x)))
        snippetWidth <- snippetWidth - .namesW - 1
    seq_snippet <- BString.get_snippet(x[[i]], snippetWidth)
    if (!is.null(names(x)))
        seq_snippet <- format(seq_snippet, width=snippetWidth)
    cat(format(paste("[", i,"]", sep=""), width=iW, justify="right"), " ",
        format(nchar, width=ncharW, justify="right"), " ",
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
    ncharW <- max(nchar(ncharMax), nchar("nchar"))
    BStringSet.show_frame_header(iW, ncharW, !is.null(names(x)))
    if (lx <= 2*half_nrow+1) {
        for (i in seq_len(lx))
            BStringSet.show_frame_line(x, i, iW, ncharW)
    } else {
        for (i in 1:half_nrow)
            BStringSet.show_frame_line(x, i, iW, ncharW)
        cat(format("...", width=iW, justify="right"),
            format("...", width=ncharW, justify="right"),
            "...\n")
        for (i in (lx-half_nrow+1L):lx)
            BStringSet.show_frame_line(x, i, iW, ncharW)
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

