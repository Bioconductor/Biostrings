### =========================================================================
### XString objects
### -------------------------------------------------------------------------

setClass("XString",
    representation(
        "VIRTUAL",
        data="XRaw",        # contains the string data
        offset="integer",   # a single integer
        length="integer"    # a single integer
    )
)

### 4 direct "XString" derivations (no additional slot)
setClass("BString", contains="XString")
setClass("DNAString", contains="XString")
setClass("RNAString", contains="XString")
setClass("AAString", contains="XString")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "codes", "codec", "enc_lkup" and "dec_lkup" new generics.
### For internal use only. No need to export.
###

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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "alphabet" new generic.
###

setGeneric("alphabet", function(x) standardGeneric("alphabet"))

setMethod("alphabet", "BString", function(x) NULL)
setMethod("alphabet", "DNAString", function(x) DNA_ALPHABET)
setMethod("alphabet", "RNAString", function(x) RNA_ALPHABET)
setMethod("alphabet", "AAString", function(x) AA_ALPHABET)


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
        stop("incompatible XString subtypes")
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
    XRaw.read(x@data, x@offset + i, x@offset + imax,
                      dec_lkup=dec_lkup(x))
}

XString.readCodes <- function(x, i, imax=integer(0))
{
    XRaw.readInts(x@data, x@offset + i, x@offset + imax)
}

### Only used at initialization time! (XString objects are immutable)
### 'value' must be a character string (this is not checked)
XString.write <- function(x, i, imax=integer(0), value)
{
    if (missing(i) && missing(imax)) {
        nbytes <- nchar(value, type="bytes")
        if (nbytes == 0)
            return(x)
        ## Write data starting immediately after the last byte in XRaw object
        ## 'x@data' that belongs to the sequence XString object 'x' is
        ## pointing at.
        ## This is safe because XRaw.write() is protected against subscripts
        ## 'i' and 'imax' being "out of bounds".
        i <- x@length + 1L
        imax <- x@length <- x@length + nbytes
    }
    #cat(x@offset + i, " -- ", x@offset + imax, "\n", sep="")
    XRaw.write(x@data, x@offset + i, x@offset + imax, value=value,
                       enc_lkup=enc_lkup(x))
    x
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Initialization.
###

.normalize.offset <- function(offset)
{
    if (!isSingleNumber(offset))
        stop("'offset' must be a single integer")
    if (!is.integer(offset))
        offset <- as.integer(offset)
    if (offset < 0L)
        stop("'offset' must be a non-negative integer")
    offset
}

.normalize.length <- function(offset, length, data_length)
{
    if (!isSingleNumber(length))
        stop("'length' must be a single integer")
    if (!is.integer(length))
        length <- as.integer(length)
    if (length < 0L)
        stop("'length' must be a non-negative integer")
    if (offset + length > data_length)
        stop("invalid 'length'")
    length
}

setMethod("initialize", "XString",
    function(.Object, data, offset, length, check=TRUE)
    {
        if (check) {
            if (!is(data, "XRaw"))
                stop("'data' must be an XRaw object")
            offset <- .normalize.offset(offset)
            length <- .normalize.length(offset, length, length(data))
        }
        slot(.Object, "data", check=FALSE) <- data
        slot(.Object, "offset", check=FALSE) <- offset
        slot(.Object, "length", check=FALSE) <- length
        .Object
    }
)


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
    lkup <- enc_lkup(new(class, XRaw(0), 0L, 0L, check=FALSE))
    data <- charToXRaw(x, start=start, end=end, width=width, lkup=lkup, check=check)
    new(class, data, 0L, length(data), check=FALSE)
}

.XStringToXString <- function(x, start, nchar, class, check)
{
    if (check) {
        start <- normalize.start(start)
        nchar <- normalize.nchar(start, nchar, nchar(x))
    }
    start <- x@offset + start
    lkup <- getXStringSubtypeConversionLookup(class(x), class)
    if (is.null(lkup))
        return(new(class, x@data, start-1L, nchar, check=FALSE))
    data <- copySubXRaw(x@data, start=start, nchar=nchar, lkup=lkup, check=FALSE)
    new(class, data, 0L, length(data), check=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The user-friendly versatile constructors.
###

setGeneric("XString", signature="x",
    function(class, x, start=1, nchar=NA, check=TRUE) standardGeneric("XString")
)
setMethod("XString", "character",
    function(class, x, start=1, nchar=NA, check=TRUE)
        charToXString(x, start=start, width=nchar, class=class, check=check)
)
setMethod("XString", "XString",
    function(class, x, start=1, nchar=NA, check=TRUE)
        .XStringToXString(x, start, nchar, class, check)
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
### Standard generic methods
###

### Helper function used by the show() method
XString.get_snippet <- function(x, snippetWidth)
{
    if (snippetWidth < 7)
        snippetWidth <- 7
    lx <- x@length
    if (lx <= snippetWidth) {
        toString(x)
    } else {
        w1 <- (snippetWidth - 2) %/% 2
        w2 <- (snippetWidth - 3) %/% 2
        paste(XString.read(x, 1, w1),
              "...",
              XString.read(x, lx - w2 + 1, lx),
              sep="")
    }
}

setMethod("show", "XString",
    function(object)
    {
        lo <- object@length
        cat("  ", lo, "-letter \"", class(object), "\" instance", sep="")
        #if (!is.null(object@codec))
        #    cat(" with alphabet:", toString(object@codec@letters))
        cat("\nseq:", XString.get_snippet(object, getOption("width") - 5))
        cat("\n")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

setMethod("length", "XString", function(x) x@length)

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
        data <- XRaw(length(i))
        XRaw.copy(data, x@offset + i, src=x@data)
        new(class(x), data, 0L, length(data), check=FALSE)
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
    {
        stop(paste("attempt to modify the value of a \"",
                   class(x), "\" instance", sep=""))
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
    one <- as.integer(1)
    ans <- !XRaw.compare(x@data, x@offset + one, y@data, y@offset + one, x@length)
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("as.character", "XString", function(x) XString.read(x, 1, x@length))
setMethod("toString", "XString", function(x, ...) as.character(x))
setMethod("nchar", "XString", function(x, type="chars", allowNA=FALSE) x@length)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "XString.substr" function.
###
### The "XString.substr" function is very fast because it does not copy
### the string data. Return a XString object (not vectorized).
### 'start' and 'end' must be single integers verifying:
###   1 <= start <= end <= length(x)
### WARNING: This function is voluntarly unsafe (it doesn't check its
### arguments) because we want it to be the fastest possible!
### The safe (and exported) version of "XString.substr" is the "subXString"
### function.
XString.substr <- function(x, start, end)
{
    shift <- start - 1L
    new(class(x), x@data, x@offset + shift, end - shift, check=FALSE)
}

