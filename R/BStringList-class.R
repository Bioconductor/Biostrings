### =========================================================================
### The BStringList class
### -------------------------------------------------------------------------
###
### We use the "has a" rather than the "is a" approach for this class.
### The "is a" approach would be doing:
###
###   setClass("BStringList", contains="list")
###
### At first sight, this approach seems to offer the following big advantage
### over the "has a" approach: BStringList objects inherit the full "list
### semantic" out-of-the-box! In other words, the "list API" (i.e. all the
### functions/methods that work on a list, e.g. [, [[, lapply(), rev(), etc...,
### there are a lot) still work on a BStringList object.
### Unfortunately, most of the time, they don't do the right thing.
### For exampple:
### 1. The user can easily screw up a BStringList object 'x' with
###    x[1] <- something, or x[[1]] <- something
###    This needs to be prevented by redefining the "[<-" and the "[[<-"
###    methods for BStringList objects.
### 2. [, rev(), and any method that you would expect to return an object of
###    the same class as the input object, unfortunately won't. This is because
###    they tend to drop the class attribute.
###    Here again, this must be prevented by redefining all these methods.
### So the party is almost over: we end up in a situation where we need to
### redefine almost every "list" method. And there can be a lot of them, in
### many different places (standard methods + methods defined in contributed
### packages), and new ones can show up in the future... Let's face it: there
### is no chance we could guarantee that the BStringList API does (and will
### always do) the right thing.
###
### Long story short: THE "HAS A" APPROACH IS MUCH SAFER. Nothing works
### out-of-the-box (which is in fact a good thing), and new methods are added
### when the need raises. This way we have full control on the entire
### BStringList API (which should stay reasonably small), I we can make it safe.
###

setClass("BStringList",
    representation(
        seqs="list"
    )
)

### 3 direct "BStringList" derivations (no additional slot)
setClass("DNAStringList", contains="BStringList")
setClass("RNAStringList", contains="BStringList")
setClass("AAStringList", contains="BStringList")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods.
###
### Like for BStringSet and BStringViews objects, the strict minimum of
### methods that must work with BStringList objects is:
###   length, width, nchar, names
###

setMethod("length", "BStringList", function(x) length(x@seqs))

setMethod("nchar", "BStringList",
    function(x, type="chars", allowNA=FALSE)
        .Call("BStrings_to_nchars", x@seqs, PACKAGE="Biostrings")
)

setMethod("width", "BStringList", function(x) nchar(x))

setMethod("names", "BStringList", function(x) names(x@seqs))

setReplaceMethod("names", "BStringList",
    function(x, value)
    {
        names(x@seqs) <- value
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Initialization (not intended to be used directly by the user).
###

setMethod("initialize", "BStringList",
    function(.Object, seqs, use.names=TRUE, check=TRUE)
    {
        if (check) {
            if (!is.list(seqs) || !all(sapply(seqs, function(x) is(x, "BString"))))
                stop("'seqs' must be a list of BString objects")
        }
        use.names <- normalize.use.names(use.names)
        if (!use.names) names(seqs) <- NULL
        slot(.Object, "seqs", check=FALSE) <- seqs
        .Object
    }
)
setMethod("initialize", "DNAStringList",
    function(.Object, seqs, use.names=TRUE, check=TRUE)
    {
        if (check) {
            if (!is.list(seqs) || !all(sapply(seqs, function(x) is(x, "DNAString"))))
                stop("'seqs' must be a list of DNAString objects")
        }
        use.names <- normalize.use.names(use.names)
        if (!use.names) names(seqs) <- NULL
        callNextMethod(.Object, seqs, check=FALSE)
    }
)
setMethod("initialize", "RNAStringList",
    function(.Object, seqs, use.names=TRUE, check=TRUE)
    {
        if (check) {
            if (!is.list(seqs) || !all(sapply(seqs, function(x) is(x, "RNAString"))))
                stop("'seqs' must be a list of RNAString objects")
        }
        use.names <- normalize.use.names(use.names)
        if (!use.names) names(seqs) <- NULL
        callNextMethod(.Object, seqs, check=FALSE)
    }
)
setMethod("initialize", "AAStringList",
    function(.Object, seqs, use.names=TRUE, check=TRUE)
    {
        if (check) {
            if (!is.list(seqs) || !all(sapply(seqs, function(x) is(x, "AAString"))))
                stop("'seqs' must be a list of AAString objects")
        }
        use.names <- normalize.use.names(use.names)
        if (!use.names) names(seqs) <- NULL
        callNextMethod(.Object, seqs, check=FALSE)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helper functions used by the versatile constructors below.
###

.charToBStringList <- function(x, start, end, width, use.names, baseClass, check)
{
    safe_locs <- narrow(nchar(x, type="bytes"), start, end, width)
    use.names <- normalize.use.names(use.names)
    proto <- new(baseClass, XRaw(0), 0L, 0L, check=FALSE)
    data <- .Call("STRSXP_to_XRaw",
                  x, start(safe_locs), width(safe_locs), "", enc_lkup(proto),
                  PACKAGE="Biostrings")
    locs2 <- intToAdjacentIntervals(width(safe_locs))
    ans <- .Call("XRaw_to_BStringList",
                 data, start(locs2), width(locs2), proto,
                 PACKAGE="Biostrings")
    if (use.names) names(ans) <- names(x)
    ans
}

.narrowBStringList <- function(x, start, end, width, use.names, baseClass, check)
{
    safe_locs <- narrow(nchar(x), start, end, width)
    use.names <- normalize.use.names(use.names)
    class <- paste(baseClass, "List", sep="")
    if (class(x) == class)
        proto <- NULL
    else
        proto <- new(baseClass, XRaw(0), 0L, 0L, check=FALSE)
    ans <- .Call("narrow_BStringList",
                 x, start(safe_locs), width(safe_locs), proto,
                 PACKAGE="Biostrings")
    if (use.names) names(ans) <- names(x)
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The user-friendly versatile constructors.
###
### All these constructors use the uSEW (user-specified Start/End/Width)
### interface.
###

setGeneric("BStringList", signature="x",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
        standardGeneric("BStringList")
)
setGeneric("DNAStringList", signature="x",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
        standardGeneric("DNAStringList")
)
setGeneric("RNAStringList", signature="x",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
        standardGeneric("RNAStringList")
)
setGeneric("AAStringList", signature="x",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
        standardGeneric("AAStringList")
)

setMethod("BStringList", "character",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
        .charToBStringList(x, start, end, width, use.names, "BString", check)
)
setMethod("DNAStringList", "character",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
        .charToBStringList(x, start, end, width, use.names, "DNAString", check)
)
setMethod("RNAStringList", "character",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
        .charToBStringList(x, start, end, width, use.names, "RNAString", check)
)
setMethod("AAStringList", "character",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
        .charToBStringList(x, start, end, width, use.names, "AAString", check)
)

setMethod("BStringList", "BStringList",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
    {
        if (is(x, "DNAStringList") || is(x, "RNAStringList"))
            stop("not supported yet, sorry!")
        .narrowBStringList(x, start, end, width, use.names, "BString", check)
    }
)
setMethod("DNAStringList", "BStringList",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
    {
        if (is(x, "AAStringList"))
            stop("incompatible input type")
        if (!is(x, "DNAStringList") && !is(x, "RNAStringList"))
            stop("not supported yet, sorry!")
        .narrowBStringList(x, start, end, width, use.names, "DNAString", check)
    }
)
setMethod("RNAStringList", "BStringList",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
    {
        if (is(x, "AAStringList"))
            stop("incompatible input type")
        if (!is(x, "DNAStringList") && !is(x, "RNAStringList"))
            stop("not supported yet, sorry!")
        .narrowBStringList(x, start, end, width, use.names, "RNAString", check)
    }
)
setMethod("AAStringList", "BStringList",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
    {
        if (is(x, "DNAStringList") || is(x, "RNAStringList"))
            stop("incompatible input type")
        .narrowBStringList(x, start, end, width, use.names, "AAString", check)
    }
)

### By "vector", we mean at least "list".
setMethod("BStringList", "vector",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
    {
        if (!is.na(end))
            stop("non-NA 'end' not yet supported by this method, sorry!")
        x <- lapply(x, function(seq) BString(seq, start=start, width=width, check=check))
        new("BStringList", x, use.names=use.names, check=FALSE)
    }
)
setMethod("DNAStringList", "vector",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
    {
        if (!is.na(end))
            stop("non-NA 'end' not yet supported by this method, sorry!")
        x <- lapply(x, function(seq) DNAString(seq, start=start, width=width, check=check))
        new("DNAStringList", x, use.names=use.names, check=FALSE)
    }
)
setMethod("RNAStringList", "vector",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
    {
        if (!is.na(end))
            stop("non-NA 'end' not yet supported by this method, sorry!")
        x <- lapply(x, function(seq) RNAString(seq, start=start, width=width, check=check))
        new("RNAStringList", x, use.names=use.names, check=FALSE)
    }
)
setMethod("AAStringList", "vector",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
    {
        if (!is.na(end))
            stop("non-NA 'end' not yet supported by this method, sorry!")
        x <- lapply(x, function(seq) AAString(seq, start=start, width=width, check=check))
        new("AAStringList", x, use.names=use.names, check=FALSE)
    }
)

### By "ANY", we mean at least "BStringViews" and those silly "AsIs" objects
### found in the probe packages (e.g. drosophila2probe$sequence).
### Note that we rely on the "as.list" methods defined for those objects.
### Performance (on mustafa):
###   > library(drosophila2probe)
###   > dict0 <- drosophila2probe$sequence
###   > system.time(DNAStringList(dict0[1:100000], start=11L, width=as.integer(NA)))
###      user  system elapsed 
###     4.513   0.056   4.761 
###

### Workaround for handling the "AsIs" objects found in the probe packages
.normalize.x <- function(x)
    if (is.character(x)) as.character(x) else as.list(x)

setMethod("BStringList", "ANY",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
        BStringList(.normalize.x(x), start=start, end=end, width=width,
                    use.names=use.names, check=check)
)
setMethod("DNAStringList", "ANY",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
        DNAStringList(.normalize.x(x), start=start, end=end, width=width,
                      use.names=use.names, check=check)
)
setMethod("RNAStringList", "ANY",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
        RNAStringList(.normalize.x(x), start=start, end=end, width=width,
                      use.names=use.names, check=check)
)
setMethod("AAStringList", "ANY",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE, check=TRUE)
        AAStringList(.normalize.x(x), start=start, end=end, width=width,
                     use.names=use.names, check=check)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method.
###

.namesW <- 20

BStringList.show_frame_header <- function(iW, widthW, with.names)
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

BStringList.show_frame_line <- function(x, i, iW, widthW)
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
        snippet_names <- names(x)[i]
        if (is.na(snippet_names))
            snippet_names <- "<NA>"
        else if (nchar(snippet_names) > .namesW)
            snippet_names <- paste(substr(snippet_names, 1, .namesW-3), "...", sep="")
        cat(" ", snippet_names, sep="")
    }
    cat("\n")
}

### 'half_nrow' must be >= 1
BStringList.show_frame <- function(x, half_nrow=9L)
{
    lx <- length(x)
    iW <- nchar(as.character(lx)) + 2 # 2 for the brackets
    ncharMax <- max(nchar(x))
    widthW <- max(nchar(ncharMax), nchar("width"))
    BStringList.show_frame_header(iW, widthW, !is.null(names(x)))
    if (lx <= 2*half_nrow+1) {
        for (i in seq_len(lx))
            BStringList.show_frame_line(x, i, iW, widthW)
    } else {
        for (i in 1:half_nrow)
            BStringList.show_frame_line(x, i, iW, widthW)
        cat(format("...", width=iW, justify="right"),
            format("...", width=widthW, justify="right"),
            "...\n")
        for (i in (lx-half_nrow+1L):lx)
            BStringList.show_frame_line(x, i, iW, widthW)
    }
}

setMethod("show", "BStringList",
    function(object)
    {
        cat("  A ", class(object), " instance of length ", length(object), "\n", sep="")
        if (length(object) != 0)
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
        x@seqs <- x@seqs[i]
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
        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")
        if (missing(i))
            stop("subscript is missing")
        if (is.character(i))
            stop("cannot subset a ", class(x), " object by names")
        as.list(x)[[i]]
    }
)

setReplaceMethod("[[", "BStringList",
    function(x, i, j,..., value)
    {
        stop("attempt to modify the value of a ", sQuote(class(x)), " object")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Other methods.
###

setMethod("as.list", "BStringList", function(x) x@seqs)

