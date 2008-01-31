### =========================================================================
### The BStringList class
### -------------------------------------------------------------------------
###

setClass("BStringList",
    representation(
        bstrings="list"
    )
)

### 3 direct "BStringList" derivations (no additional slot)
setClass("DNAStringList", contains="BStringList")
setClass("RNAStringList", contains="BStringList")
setClass("AAStringList", contains="BStringList")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Initialization.
###

setMethod("initialize", "BStringList",
    function(.Object, bstrings)
    {
        if (!is.list(bstrings) || !all(sapply(bstrings, function(x) is(x, "BString"))))
                stop("'bstrings' must be a list of BString objects")
        .Object@bstrings <- bstrings
        .Object
    }
)

setMethod("initialize", "DNAStringList",
    function(.Object, bstrings)
    {
        if (!is.list(bstrings) || !all(sapply(bstrings, function(x) is(x, "DNAString"))))
                stop("'bstrings' must be a list of DNAString objects")
        .Object@bstrings <- bstrings
        .Object
    }
)

setMethod("initialize", "RNAStringList",
    function(.Object, bstrings)
    {
        if (!is.list(bstrings) || !all(sapply(bstrings, function(x) is(x, "RNAString"))))
                stop("'bstrings' must be a list of RNAString objects")
        .Object@bstrings <- bstrings
        .Object
    }
)

setMethod("initialize", "AAStringList",
    function(.Object, bstrings)
    {
        if (!is.list(bstrings) || !all(sapply(bstrings, function(x) is(x, "AAString"))))
                stop("'bstrings' must be a list of AAString objects")
        .Object@bstrings <- bstrings
        .Object
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Some convenience constructors.
###

setGeneric("BStringList", function(src) standardGeneric("BStringList"))
setGeneric("DNAStringList", function(src) standardGeneric("DNAStringList"))
setGeneric("RNAStringList", function(src) standardGeneric("RNAStringList"))
setGeneric("AAStringList", function(src) standardGeneric("AAStringList"))

setMethod("BStringList", "list",
    function(src)
    {
        bstrings <- lapply(src, BString)
        new("BStringList", bstrings)
    }
)
setMethod("DNAStringList", "list",
    function(src)
    {
        bstrings <- lapply(src, DNAString)
        new("DNAStringList", bstrings)
    }
)
setMethod("RNAStringList", "list",
    function(src)
    {
        bstrings <- lapply(src, RNAString)
        new("RNAStringList", bstrings)
    }
)
setMethod("AAStringList", "list",
    function(src)
    {
        bstrings <- lapply(src, AAString)
        new("AAStringList", bstrings)
    }
)

### It is important that this works fine on a character vector and on a
### BStringViews object. For those 2 types, as.list() does actually the right
### thing (note that we rely on our "as.list" method for BStringViews objects).
setMethod("BStringList", "ANY", function(src) BStringList(as.list(src)))
setMethod("DNAStringList", "ANY", function(src) DNAStringList(as.list(src)))
setMethod("RNAStringList", "ANY", function(src) RNAStringList(as.list(src)))
setMethod("AAStringList", "ANY", function(src) AAStringList(as.list(src)))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods.
###

setMethod("as.list", "BStringList", function(x) x@bstrings)

setMethod("length", "BStringList", function(x) length(as.list(x)))

setMethod("nchar", "BStringList",
    function(x, type = "chars", allowNA = FALSE)
    {
        if (length(x) == 0)
            return(integer(0))
        sapply(as.list(x), nchar)
    }
)

setGeneric("desc", function(x) standardGeneric("desc"))
setMethod("desc", "BStringList", function(x) names(as.list(x)))

setGeneric("desc<-", signature="x", function(x, value) standardGeneric("desc<-"))
setReplaceMethod("desc", "BStringList",
    function(x, value)
    {
        names(x@bstrings) <- value
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method.
###

.descW <- 20

BStringList.show_frame_header <- function(iW, ncharW, with.desc)
{
    cat(format("", width=iW+1),
        format("nchar", width=ncharW, justify="right"),
        sep="")
    if (with.desc) {
        cat(format("", width=getOption("width")-iW-ncharW-.descW-1),
            format("desc", width=.descW, justify="left"),
            sep="")
    }
    cat("\n")
}

BStringList.show_frame_line <- function(x, i, iW, ncharW)
{
    nchar <- nchar(x[[i]])
    snippetWidth <- getOption("width") - 2 - iW - ncharW
    if (!is.null(desc(x)))
        snippetWidth <- snippetWidth - .descW - 1
    cat(format(paste("[", i,"]", sep=""), width=iW, justify="right"), " ",
        format(nchar, width=ncharW, justify="right"), " ",
        BString.get_snippet(x[[i]], snippetWidth),
        sep="")
    if (!is.null(desc(x))) {
        snippetDesc <- desc(x)[i]
        if (nchar(snippetDesc) > .descW)
            snippetDesc <- paste(substr(snippetDesc, 1, .descW-3), "...", sep="")
        cat(" ", snippetDesc, sep="")
    }
    cat("\n")
}

### 'half_nrow' must be >= 1
BStringList.show_frame <- function(x, half_nrow=9L)
{
    lx <- length(x)
    if (lx == 0)
        return
    iW <- nchar(as.character(lx)) + 2 # 2 for the brackets
    ncharMax <- max(nchar(x))
    ncharW <- max(nchar(ncharMax), nchar("nchar"))
    BStringList.show_frame_header(iW, ncharW, !is.null(desc(x)))
    if (lx <= 2*half_nrow+1) {
        for (i in seq_len(lx))
            BStringList.show_frame_line(x, i, iW, ncharW)
    } else {
        for (i in 1:half_nrow)
            BStringList.show_frame_line(x, i, iW, ncharW)
        cat(format("...", width=iW, justify="right"),
            format("...", width=ncharW, justify="right"),
            "...\n")
        for (i in (lx-half_nrow+1L):lx)
            BStringList.show_frame_line(x, i, iW, ncharW)
    }
}

setMethod("show", "BStringList",
    function(object)
    {
        cat("  A ", class(object), " instance of length ", length(object), "\n", sep="")
        BStringList.show_frame(object)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

setMethod("[", "BStringList",
    function(x, i, j, ..., drop)
    {
        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")
        if (missing(i))
            return(x)
        x@bstrings <- x@bstrings[i]
        x
    }
)

setReplaceMethod("[", "BStringList",
    function(x, i, j,..., value)
    {
        stop("attempt to modify the value of a ", sQuote(class(x)), " object")
    }
)

### Extract the i-th element of a BStringList object as a BString object.
setMethod("[[", "BStringList",
    function(x, i, j, ...)
    {
        if (!missing(j) || length(list(...)) > 0 || missing(i))
            stop("invalid subsetting")
        as.list(x)[[i]]
    }
)

setReplaceMethod("[[", "BStringList",
    function(x, i, j,..., value)
    {
        stop("attempt to modify the value of a ", sQuote(class(x)), " object")
    }
)

