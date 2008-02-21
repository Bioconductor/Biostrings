### =========================================================================
### The BStringSet class
### -------------------------------------------------------------------------
###

setClass("BStringSet",
    contains="SubstrLocs",
    representation(
        super="BString"
    )
)

### 3 direct "BStringSet" derivations (no additional slot)
setClass("DNAStringSet", contains="BStringSet")
setClass("RNAStringSet", contains="BStringSet")
setClass("AAStringSet", contains="BStringSet")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Initialization.
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
### Accessor methods.
###

setGeneric("super", function(x) standardGeneric("super"))
setMethod("super", "BStringSet", function(x) x@super)


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
        cat(format("", width=getOption("width")-iW-ncharW-.namesW-1),
            format("names", width=.namesW, justify="left"),
            sep="")
    }
    cat("\n")
}

BStringSet.show_frame_line <- function(x, i, iW, ncharW)
{
    nchar <- nchar(x[[i]])
    snippetWidth <- getOption("width") - 2 - iW - ncharW
    if (!is.null(names(x)))
        snippetWidth <- snippetWidth - .namesW - 1
    cat(format(paste("[", i,"]", sep=""), width=iW, justify="right"), " ",
        format(nchar, width=ncharW, justify="right"), " ",
        BString.get_snippet(x[[i]], snippetWidth),
        sep="")
    if (!is.null(names(x))) {
        snippetDesc <- names(x)[i]
        if (nchar(snippetDesc) > .namesW)
            snippetDesc <- paste(substr(snippetDesc, 1, .namesW-3), "...", sep="")
        cat(" ", snippetDesc, sep="")
    }
    cat("\n")
}

### 'half_nrow' must be >= 1
BStringSet.show_frame <- function(x, half_nrow=9L)
{
    lx <- length(x)
    if (lx == 0)
        return
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
        if (missing(i) || !is.numeric(i))
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

