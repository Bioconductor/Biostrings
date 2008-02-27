### =========================================================================
### The BString class
### -------------------------------------------------------------------------

setClass("BString",
    representation(
        data="XRaw",        # contains the string data
        offset="integer",   # a single integer
        length="integer"    # a single integer
    )
)

### 3 direct "BString" derivations (no additional slot)
setClass("DNAString", contains="BString")
setClass("RNAString", contains="BString")
setClass("AAString", contains="BString")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "codes", "codec", "enc_lkup" and "dec_lkup" new generics.
### For internal use only. No need to export.
###

setGeneric("codes", signature="x", function(x, ...) standardGeneric("codes"))
setMethod("codes", "BString", function(x, ...) 0:255)
setMethod("codes", "DNAString", function(x, baseOnly=FALSE) DNAcodes(baseOnly))
setMethod("codes", "RNAString", function(x, baseOnly=FALSE) RNAcodes(baseOnly))

setGeneric("codec", function(x) standardGeneric("codec"))
setMethod("codec", "BString", function(x) NULL)
setMethod("codec", "DNAString", function(x) DNA_STRING_CODEC)
setMethod("codec", "RNAString", function(x) RNA_STRING_CODEC)

setGeneric("enc_lkup", function(x) standardGeneric("enc_lkup"))
setMethod("enc_lkup", "BString",
    function(x)
    {
        codec <- codec(x)
        if (is.null(codec)) NULL else codec@enc_lkup
    }
)

setGeneric("dec_lkup", function(x) standardGeneric("dec_lkup"))
setMethod("dec_lkup", "BString",
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
### The "BString.read", "BString.write" and "BString.readInts" functions.
### NOT exported!
###

BString.read <- function(x, i, imax=integer(0))
{
    XRaw.read(x@data, x@offset + i, x@offset + imax,
                      dec_lkup=dec_lkup(x))
}

BString.readInts <- function(x, i, imax=integer(0))
{
    XRaw.readInts(x@data, x@offset + i, x@offset + imax)
}

### Only used at initialization time! (BString objects are immutable)
### 'value' must be a character string (this is not checked)
BString.write <- function(x, i, imax=integer(0), value)
{
    if (missing(i) && missing(imax)) {
        nbytes <- nchar(value, type="bytes")
        if (nbytes == 0)
            return(x)
        ## Write data starting immediately after the last byte in XRaw object
        ## 'x@data' that belongs to the sequence BString object 'x' is
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

setMethod("initialize", "BString",
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
setMethod("initialize", "DNAString",
    function(.Object, data, offset, length, check=TRUE)
        callNextMethod(.Object, data, offset, length, check=check)
)
setMethod("initialize", "RNAString",
    function(.Object, data, offset, length, check=TRUE)
        callNextMethod(.Object, data, offset, length, check=check)
)
setMethod("initialize", "AAString",
    function(.Object, data, offset, length, check=TRUE)
        callNextMethod(.Object, data, offset, length, check=check)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helper functions used by the versatile constructors below.
###

.charseqToBString <- function(seq, start, nchar, class, check)
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
    proto <- new(class, XRaw(0), 0L, 0L, check=FALSE)
    .Call("charseq_to_BString",
          seq, start, nchar, enc_lkup(proto), proto,
          PACKAGE="Biostrings")
}

copyDataOnBStringConversion <- function(class1, class2)
{
    is_nucl1 <- extends(class1, "DNAString") || extends(class1, "RNAString")
    is_nucl2 <- extends(class2, "DNAString") || extends(class2, "RNAString")
    is_nucl1 != is_nucl2
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The user-friendly versatile constructors.
###

setGeneric("BString", signature="seq",
    function(seq, start=1, nchar=NA, check=TRUE) standardGeneric("BString"))
setGeneric("DNAString", signature="seq",
    function(seq, start=1, nchar=NA, check=TRUE) standardGeneric("DNAString"))
setGeneric("RNAString", signature="seq",
    function(seq, start=1, nchar=NA, check=TRUE) standardGeneric("RNAString"))
setGeneric("AAString", signature="seq",
    function(seq, start=1, nchar=NA, check=TRUE) standardGeneric("AAString"))

setMethod("BString", "character",
    function(seq, start=1, nchar=NA, check=TRUE)
    {
        .charseqToBString(seq, start, nchar, "BString", check)
    }
)
setMethod("DNAString", "character",
    function(seq, start=1, nchar=NA, check=TRUE)
    {
        .charseqToBString(seq, start, nchar, "DNAString", check)
    }
)
setMethod("RNAString", "character",
    function(seq, start=1, nchar=NA, check=TRUE)
    {
        .charseqToBString(seq, start, nchar, "RNAString", check)
    }
)
setMethod("AAString", "character",
    function(seq, start=1, nchar=NA, check=TRUE)
    {
        .charseqToBString(seq, start, nchar, "AAString", check)
    }
)

setMethod("BString", "BString",
    function(seq, start=1, nchar=NA, check=TRUE)
    {
        if (check) {
            start <- normalize.start(start)
            nchar <- normalize.nchar(start, nchar, nchar(seq))
        }
        start <- seq@offset + start
        if (!copyDataOnBStringConversion(class(seq), "BString"))
            return(new("BString", seq@data, start-1L, nchar, check=FALSE))
        lkup <- dec_lkup(seq)
        data <- copySubXRaw(seq@data, start=start, nchar=nchar, lkup=lkup, check=FALSE)
        new("BString", data, 0L, length(data), check=FALSE)
    }
)
setMethod("DNAString", "BString",
    function(seq, start=1, nchar=NA, check=TRUE)
    {
        if (is(seq, "AAString"))
            stop("incompatible input type")
        if (check) {
            start <- normalize.start(start)
            nchar <- normalize.nchar(start, nchar, nchar(seq))
        }
        start <- seq@offset + start
        if (!copyDataOnBStringConversion(class(seq), "DNAString"))
            return(new("DNAString", seq@data, start-1L, nchar, check=FALSE))
        lkup <- DNA_STRING_CODEC@enc_lkup
        data <- copySubXRaw(seq@data, start=start, nchar=nchar, lkup=lkup, check=FALSE)
        new("DNAString", data, 0L, length(data), check=FALSE)
    }
)
setMethod("RNAString", "BString",
    function(seq, start=1, nchar=NA, check=TRUE)
    {
        if (is(seq, "AAString"))
            stop("incompatible input type")
        if (check) {
            start <- normalize.start(start)
            nchar <- normalize.nchar(start, nchar, nchar(seq))
        }
        start <- seq@offset + start
        if (!copyDataOnBStringConversion(class(seq), "RNAString"))
            return(new("RNAString", seq@data, start-1L, nchar, check=FALSE))
        lkup <- RNA_STRING_CODEC@enc_lkup
        data <- copySubXRaw(seq@data, start=start, nchar=nchar, lkup=lkup, check=FALSE)
        new("RNAString", data, 0L, length(data), check=FALSE)
    }
)
setMethod("AAString", "BString",
    function(seq, start=1, nchar=NA, check=TRUE)
    {
        if (is(seq, "DNAString") || is(seq, "RNAString"))
            stop("incompatible input type")
        if (check) {
            start <- normalize.start(start)
            nchar <- normalize.nchar(start, nchar, nchar(seq))
        }
        start <- seq@offset + start
        new("AAString", seq@data, start-1L, nchar, check=FALSE)
    }
)

### Not exported
mkBString <- function(class, seq)
{
    do.call(class, list(seq=seq))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Standard generic methods
###

### Helper function used by the show() method
BString.get_snippet <- function(x, snippetWidth)
{
    if (snippetWidth < 7)
        snippetWidth <- 7
    lx <- x@length
    if (lx <= snippetWidth) {
        toString(x)
    } else {
        w1 <- (snippetWidth - 2) %/% 2
        w2 <- (snippetWidth - 3) %/% 2
        paste(BString.read(x, 1, w1),
              "...",
              BString.read(x, lx - w2 + 1, lx),
              sep="")
    }
}

setMethod("show", "BString",
    function(object)
    {
        lo <- object@length
        cat("  ", lo, "-letter \"", class(object), "\" instance", sep="")
        #if (!is.null(object@codec))
        #    cat(" with alphabet:", toString(object@codec@letters))
        cat("\nseq:", BString.get_snippet(object, getOption("width") - 5))
        cat("\n")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

setMethod("length", "BString", function(x) x@length)

setMethod("[", "BString",
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
        ## class(x) can be "BString" or one of its derivations: "DNAString",
        ## "RNAString" or "AAString".
        new(class(x), data, 0L, length(data), check=FALSE)
    }
)

### The only reason for defining the replacement version of the "[" operator
### is to let the user know that he can't use it:
###   bs <- BString("AbnbIU")
###   bs[2] <- "X" # provokes an error
### If we don't define it, then the user can type the above and believe that
### it actually did something but it didn't.
setReplaceMethod("[", "BString",
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

comparableBStrings <- function(e1, e2)
{
    ## 'e1' and 'e2' are comparable iff they are both encoded or not
    (is(e1, "DNAString") || is(e1, "RNAString")) == (is(e2, "DNAString") || is(e2, "RNAString"))
}

### 'x' and 'y' must be BString objects
.BString.equal <- function(x, y)
{
    if (x@length != y@length)
        return(FALSE)
    one <- as.integer(1)
    ans <- !XRaw.compare(x@data, x@offset + one, y@data, y@offset + one, x@length)
    as.logical(ans)
}

setMethod("==", signature(e1="BString", e2="BString"),
    function(e1, e2)
    {
        if (!comparableBStrings(e1, e2)) {
            class1 <- class(e1)
            class2 <- class(e2)
            stop("comparison between a \"", class1, "\" instance ",
                 "and a \"", class2, "\" instance ",
                 "is not supported")
        }
        .BString.equal(e1, e2)
    }
)
setMethod("==", signature(e1="BString", e2="character"),
    function(e1, e2)
    {
        class1 <- class(e1)
        if (class1 != "BString")
            stop("comparison between a \"", class1, "\" instance ",
                 "and a character vector ",
                 "is not supported")
        if (length(e2) != 1 || e2 %in% c("", NA))
            stop("comparison between a \"BString\" instance and a character vector ",
                 "of length != 1 or an empty string or an NA ",
                 "is not supported")
        .BString.equal(e1, BString(e2))
    }
)
setMethod("==", signature(e1="character", e2="BString"),
    function(e1, e2) e2 == e1
)

setMethod("!=", signature(e1="BString", e2="BString"),
    function(e1, e2) !(e1 == e2)
)
setMethod("!=", signature(e1="BString", e2="character"),
    function(e1, e2) !(e1 == e2)
)
setMethod("!=", signature(e1="character", e2="BString"),
    function(e1, e2) !(e1 == e2)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("as.character", "BString", function(x) BString.read(x, 1, x@length))
setMethod("toString", "BString", function(x, ...) as.character(x))
setMethod("nchar", "BString", function(x, type="chars", allowNA=FALSE) x@length)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "BString.substr" function.
###
### The "BString.substr" function is very fast because it does not copy
### the string data. Return a BString object (not vectorized).
### 'start' and 'end' must be single integers verifying:
###   1 <= start <= end <= length(x)
### WARNING: This function is voluntarly unsafe (it doesn't check its
### arguments) because we want it to be the fastest possible!
### The safe (and exported) version of "BString.substr" is the "subBString"
### function.
BString.substr <- function(x, start, end)
{
    shift <- start - 1L
    new(class(x), x@data, x@offset + shift, end - shift, check=FALSE)
}

