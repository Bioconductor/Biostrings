### =========================================================================
### XString objects
### -------------------------------------------------------------------------

setClass("XString",
    representation(
        "VIRTUAL",
        xdata="XRaw",       # contains the sequence data (external)
        offset="integer",   # a single integer
        length="integer"    # a single integer
    ),
    prototype(
        #xdata=XRaw(0),     # see newEmptyXString() below for why this doesn't
                            # work
        offset=0L,
        length=0L
    )
)

### XString subtypes (direct "XString" derivations with no additional slots)
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
### because the DLL of the package needs to be loaded before XRaw() can be
### called.
###

newEmptyXString <- function(class) new(class, xdata=XRaw(0))


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
setGeneric("baseXStringSubtype", function(x) standardGeneric("baseXStringSubtype"))
setMethod("baseXStringSubtype", "BString",
    function(x) class(newEmptyXString("BString"))
)
setMethod("baseXStringSubtype", "DNAString",
    function(x) class(newEmptyXString("DNAString"))
)
setMethod("baseXStringSubtype", "RNAString",
    function(x) class(newEmptyXString("RNAString"))
)
setMethod("baseXStringSubtype", "AAString",
    function(x) class(newEmptyXString("AAString"))
)

setGeneric("codes", signature="x", function(x, ...) standardGeneric("codes"))
setMethod("codes", "XString", function(x, ...) 0:255)
setMethod("codes", "DNAString", function(x, baseOnly=FALSE) DNAcodes(baseOnly))
setMethod("codes", "RNAString", function(x, baseOnly=FALSE) RNAcodes(baseOnly))

setGeneric("codec", function(x) standardGeneric("codec"))
setMethod("codec", "XString", function(x) NULL)
setMethod("codec", "DNAString", function(x) DNA_STRING_CODEC)
setMethod("codec", "RNAString", function(x) RNA_STRING_CODEC)

setGeneric("enc_lkup", function(x) standardGeneric("enc_lkup"))
setMethod("enc_lkup", "XString",
    function(x)
    {
        codec <- codec(x)
        if (is.null(codec)) NULL else codec@enc_lkup
    }
)

setGeneric("dec_lkup", function(x) standardGeneric("dec_lkup"))
setMethod("dec_lkup", "XString",
    function(x)
    {
        codec <- codec(x)
        if (is.null(codec)) NULL else codec@dec_lkup
    }
)

### exported
setGeneric("alphabet", function(x) standardGeneric("alphabet"))
setMethod("alphabet", "BString", function(x) NULL)
setMethod("alphabet", "DNAString", function(x) DNA_ALPHABET)
setMethod("alphabet", "RNAString", function(x) RNA_ALPHABET)
setMethod("alphabet", "AAString", function(x) AA_ALPHABET)

setMethod("length", "XString", function(x) x@length)

setMethod("nchar", "XString", function(x, type="chars", allowNA=FALSE) length(x))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Some restrictions apply for converting from an XString subtype to
### another or for comparing XString objects of different subtypes. This is
### due to the fact that different XString subtypes can use different
### encodings for their data (or no encoding at all) or simply to the fact
### that the conversion or comparison doesn't make sense from a biological
### perspective.
### The helper functions below are used internally (they are NOT exported) to
### determine those restrictions.
###

compatibleXStringSubtypes <- function(class1, class2)
{
    if (extends(class1, "DNAString") || extends(class1, "RNAString"))
        return(!extends(class2, "AAString"))
    if (extends(class1, "AAString"))
        return(!(extends(class2, "DNAString") || extends(class2, "RNAString")))
    TRUE
}

getXStringSubtypeConversionLookup <- function(from_class, to_class)
{
    if (!compatibleXStringSubtypes(from_class, to_class))
        stop("incompatible XString/XStringSet subtypes")
    from_nucleo <- extends(from_class, "DNAString") || extends(from_class, "RNAString")
    to_nucleo <- extends(to_class, "DNAString") || extends(to_class, "RNAString")
    if (from_nucleo == to_nucleo)
        return(NULL)
    if (extends(to_class, "DNAString"))
        return(DNA_STRING_CODEC@enc_lkup)
    if (extends(to_class, "RNAString"))
        return(RNA_STRING_CODEC@enc_lkup)
    if (extends(from_class, "DNAString"))
        return(DNA_STRING_CODEC@dec_lkup)
    if (extends(from_class, "RNAString"))
        return(RNA_STRING_CODEC@dec_lkup)
    stop("Biostrings internal error, please report") # should never happen
}

comparableXStrings <- function(x1, x2)
{
    is_nucleo1 <- is(x1, "DNAString") || is(x1, "RNAString")
    is_nucleo2 <- is(x2, "DNAString") || is(x2, "RNAString")
    is_nucleo1 == is_nucleo2
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "XString.read", "XString.readCodes" and "XString.write" functions.
### NOT exported!
###

XString.read <- function(x, i, imax=integer(0))
{
    XRaw.read(x@xdata, x@offset + i, x@offset + imax,
                      dec_lkup=dec_lkup(x))
}

XString.readCodes <- function(x, i, imax=integer(0))
{
    XRaw.readInts(x@xdata, x@offset + i, x@offset + imax)
}

### Only used at initialization time! (XString objects are immutable.)
### 'value' must be a character string (this is not checked).
XString.write <- function(x, i, imax=integer(0), value)
{
    if (missing(i) && missing(imax)) {
        nbytes <- nchar(value, type="bytes")
        if (nbytes == 0)
            return(x)
        ## Write data starting immediately after the last byte in XRaw object
        ## 'x@xdata' that belongs to the sequence XString object 'x' is
        ## pointing at.
        ## This is safe because XRaw.write() is protected against subscripts
        ## 'i' and 'imax' being "out of bounds".
        i <- x@length + 1L
        imax <- x@length <- x@length + nbytes
    }
    #cat(x@offset + i, " -- ", x@offset + imax, "\n", sep="")
    XRaw.write(x@xdata, x@offset + i, x@offset + imax, value=value,
                       enc_lkup=enc_lkup(x))
    x
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "XString.append" function.
### NOT exported!
###

XString.append <- function(x, y)
{
    ans_class <- baseXStringSubtype(x)
    if (baseXStringSubtype(y) != ans_class)
        stop("'x' and 'y' must be XString objects of the same subtype")
    ans_xdata <- XRaw.append(x@xdata, x@offset + 1L, x@length,
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

charToXString <- function(x, start=NA, end=NA, width=NA, class="BString", check=TRUE)
{
    if (check) {
        if (length(x) == 0)
            stop("no input sequence")
        if (length(x) > 1)
            stop("more than one input sequence")
    }
    lkup <- enc_lkup(newEmptyXString(class))
    xdata <- charToXRaw(x, start=start, end=end, width=width, lkup=lkup, check=check)
    new(class, xdata=xdata, length=length(xdata))
}

.XStringToXString <- function(x, start, nchar, class, check)
{
    if (check) {
        start <- normargStart(start)
        nchar <- normargNchar(start, nchar, nchar(x))
    }
    start <- x@offset + start
    lkup <- getXStringSubtypeConversionLookup(class(x), class)
    if (is.null(lkup))
        return(new(class, xdata=x@xdata, offset=start-1L, length=nchar))
    xdata <- copySubXRaw(x@xdata, start=start, nchar=nchar, lkup=lkup, check=FALSE)
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
### The "XString.substr" function (NOT exported).
###
### The "XString.substr" function is very fast because it does not copy
### the sequence data. Return an XString object (not vectorized).
### 'start' and 'end' must be single integers verifying:
###   1 <= start AND end <= length(x) AND start <= end + 1
### WARNING: This function is voluntarly unsafe (it doesn't check its
### arguments) because we want it to be the fastest possible!
### The safe (and exported) version of "XString.substr" is the "subseq"
### method for XString objects below.
###

XString.substr <- function(x, start, end)
{
    shift <- start - 1L
    new(class(x), xdata=x@xdata, offset=x@offset+shift, length=end-shift)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "subseq" generic and method for XString objects.
###

setGeneric("subseq", signature="x",
    function(x, start=NA, end=NA, width=NA) standardGeneric("subseq")
)

setMethod("subseq", "XString",
    function(x, start=NA, end=NA, width=NA)
    {
        limits <- new("IRanges", start=1L, width=length(x), check=FALSE)
        limits <- narrow(limits, start=start, end=end, width=width)
        XString.substr(x, start(limits), end(limits))
    }
)


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
        xdata <- XRaw(length(i))
        XRaw.copy(xdata, x@offset + i, src=x@xdata)
        new(class(x), xdata=xdata, length=length(xdata))
    }
)

### The only reason for defining the replacement version of the "[" operator
### is to let the user know that he can't use it:
###   bs <- BString("AbnbIU")
###   bs[2] <- "X" # provokes an error
### If we don't define it, then the user can type the above and believe that
### it actually did something but it didn't.
setReplaceMethod("[", "XString",
    function(x, i, j,..., value)
        stop("attempt to modify the value of a ", class(x), " instance")
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
###   dnav <- views(dna, 1:7, 101:107)
###   dnav[[1]] == dnav[[7]] # TRUE
###   dnav <- views(dna, 1:7, (length(dna)-6):length(dna))
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
    ans <- !XRaw.compare(x@xdata, x@offset + 1L, y@xdata, y@offset + 1L, x@length)
    as.logical(ans)
}

setMethod("==", signature(e1="XString", e2="XString"),
    function(e1, e2)
    {
        if (!comparableXStrings(e1, e2)) {
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

