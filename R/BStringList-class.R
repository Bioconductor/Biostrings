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
### Initialization.
###

setMethod("initialize", "BStringList",
    function(.Object, seqs, check=TRUE)
    {
        if (check) {
            if (!is.list(seqs) || !all(sapply(seqs, function(x) is(x, "BString"))))
                stop("'seqs' must be a list of BString objects")
        }
        slot(.Object, "seqs", check=FALSE) <- seqs
        .Object
    }
)
setMethod("initialize", "DNAStringList",
    function(.Object, seqs, check=TRUE)
    {
        if (check) {
            if (!is.list(seqs) || !all(sapply(seqs, function(x) is(x, "DNAString"))))
                stop("'seqs' must be a list of DNAString objects")
        }
        callNextMethod(.Object, seqs, check=FALSE)
    }
)
setMethod("initialize", "RNAStringList",
    function(.Object, seqs, check=TRUE)
    {
        if (check) {
            if (!is.list(seqs) || !all(sapply(seqs, function(x) is(x, "RNAString"))))
                stop("'seqs' must be a list of RNAString objects")
        }
        callNextMethod(.Object, seqs, check=FALSE)
    }
)
setMethod("initialize", "AAStringList",
    function(.Object, seqs, check=TRUE)
    {
        if (check) {
            if (!is.list(seqs) || !all(sapply(seqs, function(x) is(x, "AAString"))))
                stop("'seqs' must be a list of AAString objects")
        }
        callNextMethod(.Object, seqs, check=FALSE)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Core constructors (not exported, used by the exported constructors below).
###

charseqsToBStringList <- function(seqs, start=1L, nchar=NA, class="BString", check=TRUE)
{
    if (check) {
        ## Only limited checking here, more is done at the C level 
        if (!isSingleNumber(start))
            stop("'start' must be a single integer")
        if (!is.integer(start))
            start <- as.integer(start)
        if (!isSingleNumberOrNA(nchar))
            stop("'nchar' must be a single integer or NA")
        if (!is.integer(nchar))
            nchar <- as.integer(nchar)
    }
    proto <- new(class, data=XRaw(0), 0L, 0L, check=FALSE)
    seqs <- .Call("charseqs_to_BStringList",
                  seqs, start, nchar, enc_lkup(proto), proto,
                  PACKAGE="Biostrings")
    new(paste(class, "List", sep=""), seqs, check=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The exported constructors.
###

setGeneric("BStringList", signature="seqs",
    function(seqs, start=1, nchar=NA, check=TRUE) standardGeneric("BStringList"))
setGeneric("DNAStringList", signature="seqs",
    function(seqs, start=1, nchar=NA, check=TRUE) standardGeneric("DNAStringList"))
setGeneric("RNAStringList", signature="seqs",
    function(seqs, start=1, nchar=NA, check=TRUE) standardGeneric("RNAStringList"))
setGeneric("AAStringList", signature="seqs",
    function(seqs, start=1, nchar=NA, check=TRUE) standardGeneric("AAStringList"))

setMethod("BStringList", "character",
    function(seqs, start=1, nchar=NA, check=TRUE)
    {
        charseqsToBStringList(seqs, start=start, nchar=nchar, class="BString", check=check)
    }
)
setMethod("DNAStringList", "character",
    function(seqs, start=1, nchar=NA, check=TRUE)
    {
        charseqsToBStringList(seqs, start=start, nchar=nchar, class="DNAString", check=check)
    }
)
setMethod("RNAStringList", "character",
    function(seqs, start=1, nchar=NA, check=TRUE)
    {
        charseqsToBStringList(seqs, start=start, nchar=nchar, class="RNAString", check=check)
    }
)
setMethod("AAStringList", "character",
    function(seqs, start=1, nchar=NA, check=TRUE)
    {
        charseqsToBStringList(seqs, start=start, nchar=nchar, class="AAString", check=check)
    }
)

### By "vector", we mean at least "list".
setMethod("BStringList", "vector",
    function(seqs, start=1, nchar=NA, check=TRUE)
    {
        seqs <- lapply(seqs, function(seq) BString(seq, start=start, nchar=nchar, check=check))
        new("BStringList", seqs, check=FALSE)
    }
)
setMethod("DNAStringList", "vector",
    function(seqs, start=1, nchar=NA, check=TRUE)
    {
        seqs <- lapply(seqs, function(seq) DNAString(seq, start=start, nchar=nchar, check=check))
        new("DNAStringList", seqs, check=FALSE)
    }
)
setMethod("RNAStringList", "vector",
    function(seqs, start=1, nchar=NA, check=TRUE)
    {
        seqs <- lapply(seqs, function(seq) RNAString(seq, start=start, nchar=nchar, check=check))
        new("RNAStringList", seqs, check=FALSE)
    }
)
setMethod("AAStringList", "vector",
    function(seqs, start=1, nchar=NA, check=TRUE)
    {
        seqs <- lapply(seqs, function(seq) AAString(seq, start=start, nchar=nchar, check=check))
        new("AAStringList", seqs, check=FALSE)
    }
)

### By "ANY", we mean at least "BStringList", "BStringViews" and those silly
### "AsIs" objects found in the probe packages (e.g. drosophila2probe$sequence).
### Note that we rely on the "as.list" methods defined for those objects.
### Performance (on mustafa):
###   > library(drosophila2probe)
###   > dict0 <- drosophila2probe$sequence
###   > system.time(DNAStringList(dict0[1:100000], start=11L, nchar=as.integer(NA)))
###      user  system elapsed 
###     4.513   0.056   4.761 
###

### Workaround for handling the "AsIs" objects found in the probe packages
.normalize.seqs <- function(seqs)
    if (is.character(seqs)) as.character(seqs) else as.list(seqs)

setMethod("BStringList", "ANY",
    function(seqs, start=1, nchar=NA, check=TRUE)
        BStringList(.normalize.seqs(seqs), start=start, nchar=nchar, check=check)
)
setMethod("DNAStringList", "ANY",
    function(seqs, start=1, nchar=NA, check=TRUE)
        DNAStringList(.normalize.seqs(seqs), start=start, nchar=nchar, check=check)
)
setMethod("RNAStringList", "ANY",
    function(seqs, start=1, nchar=NA, check=TRUE)
        RNAStringList(.normalize.seqs(seqs), start=start, nchar=nchar, check=check)
)
setMethod("AAStringList", "ANY",
    function(seqs, start=1, nchar=NA, check=TRUE)
        AAStringList(.normalize.seqs(seqs), start=start, nchar=nchar, check=check)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods.
###

setMethod("as.list", "BStringList", function(x) x@seqs)

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
        names(x@seqs) <- value
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

