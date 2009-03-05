### =========================================================================
### XString objects
### -------------------------------------------------------------------------
###
### The XSequence class is a general container for storing an "external
### string".
###

setClass("XString", contains="XRaw", representation("VIRTUAL"))

### XString subclasses (no additional slots)
setClass("BString", contains="XString")
setClass("DNAString", contains="XString")
setClass("RNAString", contains="XString")
setClass("AAString", contains="XString")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "newEmptyXString" constructor.
### For internal use only. No need to export.
###
### Note that this cannot be made the prototype part of the XString class
### definition (and trying to do so will cause an error at installation time)
### because the DLL of the package needs to be loaded before RawPtr() can be
### called.
###

newEmptyXString <- function(class) new(class, xdata=RawPtr(0))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The core XString API.
###
### The core XString API is the strict minimal set of methods that must work
### for XString, XStringSet, XStringViews and MaskedXString objects.
### It currently consists of the following methods:
###   o NOT exported: xsbasetype
###   o exported: length, nchar
###

### NOT exported
setMethod("xsbasetype", "BString", function(x) "B")
setMethod("xsbasetype", "DNAString", function(x) "DNA")
setMethod("xsbasetype", "RNAString", function(x) "RNA")
setMethod("xsbasetype", "AAString", function(x) "AA")

### exported
setMethod("length", "XString", function(x) x@length)

setMethod("nchar", "XString", function(x, type="chars", allowNA=FALSE) length(x))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "XString.read", "XString.readCodes" and "XString.write" functions.
### NOT exported!
###

XString.read <- function(x, i, imax=integer(0))
{
    RawPtr.read(x@xdata, x@offset + i, x@offset + imax,
                      dec_lkup=xs_dec_lkup(x))
}

XString.readCodes <- function(x, i, imax=integer(0))
{
    RawPtr.readInts(x@xdata, x@offset + i, x@offset + imax)
}

### Only used at initialization time! (XString objects are immutable.)
### 'value' must be a character string (this is not checked).
XString.write <- function(x, i, imax=integer(0), value)
{
    if (missing(i) && missing(imax)) {
        nbytes <- nchar(value, type="bytes")
        if (nbytes == 0)
            return(x)
        ## Write data starting immediately after the last byte in RawPtr object
        ## 'x@xdata' that belongs to the sequence XString object 'x' is
        ## pointing at.
        ## This is safe because RawPtr.write() is protected against subscripts
        ## 'i' and 'imax' being "out of bounds".
        i <- x@length + 1L
        imax <- x@length <- x@length + nbytes
    }
    #cat(x@offset + i, " -- ", x@offset + imax, "\n", sep="")
    RawPtr.write(x@xdata, x@offset + i, x@offset + imax, value=value,
                       enc_lkup=xs_enc_lkup(x))
    x
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "XString.append" function.
### NOT exported!
###

XString.append <- function(x, y)
{
    ans_basetype <- xsbasetype(x)
    if (xsbasetype(y) != ans_basetype)
        stop("'x' and 'y' must be XString objects of the same base type")
    ans_class <- paste(ans_basetype, "String", sep="")
    ans_xdata <- RawPtr.append(x@xdata, x@offset + 1L, x@length,
                               y@xdata, y@offset + 1L, y@length)
    new(ans_class, xdata=ans_xdata, length=length(ans_xdata))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

setMethod("as.character", "XString", function(x) XString.read(x, 1, x@length))

setMethod("toString", "XString", function(x, ...) as.character(x))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helper functions used by the versatile constructors below.
###

.charToRawPtr <- function(x, start=NA, end=NA, width=NA, collapse=NULL, lkup=NULL, check=TRUE)
{
    solved_SEW <- solveUserSEW(nchar(x, type="bytes"),
                               start=start, end=end, width=width)
    .Call("new_RawPtr_from_STRSXP",
          x, start(solved_SEW), width(solved_SEW), collapse, lkup,
          PACKAGE="Biostrings")
}

.copySubRawPtr <- function(x, start=1, nchar=NA, lkup=NULL, check=TRUE)
{
    if (check) {
        start <- normargStart(start)
        nchar <- normargNchar(start, nchar, length(x))
    }
    ans <- RawPtr(nchar)
    RawPtr.copy(ans, start, start + nchar - 1L, src=x, lkup=lkup)
}

charToXString <- function(x, start=NA, end=NA, width=NA, class="BString", check=TRUE)
{
    if (check) {
        if (length(x) == 0)
            stop("no input sequence")
        if (length(x) > 1)
            stop("more than one input sequence")
    }
    lkup <- xs_enc_lkup(newEmptyXString(class))
    xdata <- .charToRawPtr(x, start=start, end=end, width=width, lkup=lkup, check=check)
    new(class, xdata=xdata, length=length(xdata))
}

.XStringToXString <- function(x, start, nchar, class, check)
{
    if (check) {
        start <- normargStart(start)
        nchar <- normargNchar(start, nchar, nchar(x))
    }
    start <- x@offset + start
    to_basetype <- substr(class, 1, nchar(class)-6)  # remove "String" suffix
    lkup <- get_xsbasetypes_conversion_lookup(xsbasetype(x), to_basetype)
    if (is.null(lkup))
        return(new(class, xdata=x@xdata, offset=start-1L, length=nchar))
    xdata <- .copySubRawPtr(x@xdata, start=start, nchar=nchar, lkup=lkup, check=FALSE)
    new(class, xdata=xdata, length=length(xdata))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The user-friendly versatile constructors.
###

setGeneric("XString", signature="x",
    function(class, x, start=1, nchar=NA, check=TRUE) standardGeneric("XString")
)

setMethod("XString", "character",
    function(class, x, start=1, nchar=NA, check=TRUE)
    {
        if (is.null(class))
            class <- "BString"
        charToXString(x, start=start, width=nchar, class=class, check=check)
    }
)

### Just because of the silly "AsIs" objects found in the probe packages
### (e.g. drosophila2probe$sequence)
setMethod("XString", "AsIs",
    function(class, x, start=1, nchar=NA, check=TRUE)
    {
        if (!is.character(x))
            stop("unsuported input type")
        class(x) <- "character" # keeps the names (unlike as.character())
        XString(class, x, start=start, nchar=nchar, check=check)
    }
)

setMethod("XString", "XString",
    function(class, x, start=1, nchar=NA, check=TRUE)
    {
        if (is.null(class))
            class <- class(x)
        .XStringToXString(x, start, nchar, class, check)
    }
)

BString <- function(x, start=1, nchar=NA, check=TRUE)
    XString("BString", x, start=start, nchar=nchar, check=check)

DNAString <- function(x, start=1, nchar=NA, check=TRUE)
    XString("DNAString", x, start=start, nchar=nchar, check=check)

RNAString <- function(x, start=1, nchar=NA, check=TRUE)
    XString("RNAString", x, start=start, nchar=nchar, check=check)

AAString <- function(x, start=1, nchar=NA, check=TRUE)
    XString("AAString", x, start=start, nchar=nchar, check=check)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method.
###

### 'x' must be an XString or MaskedXString object.
toSeqSnippet <- function(x, width)
{
    if (width < 7L)
        width <- 7L
    seqlen <- length(x)
    if (seqlen <= width) {
        as.character(x)
    } else {
        w1 <- (width - 2) %/% 2
        w2 <- (width - 3) %/% 2
        paste(as.character(subseq(x, start=1, width=w1)),
              "...",
              as.character(subseq(x, end=seqlen, width=w2)),
              sep="")
    }
}

setMethod("show", "XString",
    function(object)
    {
        lo <- object@length
        cat("  ", lo, "-letter \"", class(object), "\" instance\n", sep="")
        cat("seq:", toSeqSnippet(object, getOption("width") - 5))
        cat("\n")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

setMethod("[", "XString",
    function(x, i, j, ..., drop)
    {
        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")
        if (missing(i))
            return(x)
        if (!is.numeric(i) || any(is.na(i)))
            stop("invalid subsetting")
        if (any(i < 1) || any(i > length(x)))
            stop("subscript out of bounds")
        xdata <- RawPtr(length(i))
        RawPtr.copy(xdata, x@offset + i, src=x@xdata)
        x@xdata <- xdata
        x@offset <- 0L
        x@length <- length(xdata)
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Equality.
###
### We want:
###   BString("ab") == "ab" # TRUE
###   DNAString("TG") == RNAString("UG") # TRUE!!!
###   library(BSgenome.Hsapiens.UCSC.hg18)
###   dna <- Hsapiens$chr1
###   dna != Hsapiens$chr1 # FALSE
###   dnav <- Views(dna, start=1:7, end=101:107)
###   dnav[[1]] == dnav[[7]] # TRUE
###   dnav <- Views(dna, start=1:7, end=(length(dna)-6):length(dna))
### This is fast:
###   dnav[[1]] == dnav[[7]] # FALSE
### But this would have killed your machine:
###   s1 <- toString(dnav[[1]])
###   s7 <- toString(dnav[[7]])
###   s1 == s7

### 'x' and 'y' must be XString objects
.XString.equal <- function(x, y)
{
    if (x@length != y@length)
        return(FALSE)
    ans <- !RawPtr.compare(x@xdata, x@offset + 1L, y@xdata, y@offset + 1L, x@length)
    as.logical(ans)
}

setMethod("==", signature(e1="XString", e2="XString"),
    function(e1, e2)
    {
        if (!comparable_xsbasetypes(xsbasetype(e1), xsbasetype(e2))) {
            class1 <- class(e1)
            class2 <- class(e2)
            stop("comparison between a \"", class1, "\" instance ",
                 "and a \"", class2, "\" instance ",
                 "is not supported")
        }
        .XString.equal(e1, e2)
    }
)
setMethod("==", signature(e1="BString", e2="character"),
    function(e1, e2)
    {
        if (length(e2) != 1 || e2 %in% c("", NA))
            stop("comparison between a \"BString\" object and a character vector ",
                 "of length != 1 or an empty string or an NA ",
                 "is not supported")
        .XString.equal(e1, BString(e2))
    }
)
setMethod("==", signature(e1="character", e2="BString"),
    function(e1, e2) e2 == e1
)

setMethod("!=", signature(e1="XString", e2="XString"),
    function(e1, e2) !(e1 == e2)
)
setMethod("!=", signature(e1="BString", e2="character"),
    function(e1, e2) !(e1 == e2)
)
setMethod("!=", signature(e1="character", e2="BString"),
    function(e1, e2) !(e1 == e2)
)

