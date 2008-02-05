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
setMethod("codes", "DNAString", function(x, baseOnly=FALSE) DNAcodes(baseOnly))
setMethod("codes", "RNAString", function(x, baseOnly=FALSE) RNAcodes(baseOnly))

setGeneric("codec", function(x) standardGeneric("codec"))
setMethod("codec", "BString", function(x) NULL)
setMethod("codec", "DNAString", function(x) DNA_STRING_CODEC)
setMethod("codec", "RNAString", function(x) RNA_STRING_CODEC)

setGeneric("dec_lkup", function(x) standardGeneric("dec_lkup"))
setMethod("dec_lkup", "BString",
    function(x)
    {
        codec <- codec(x)
        if (is.null(codec)) NULL else codec@dec_lkup
    }
)

setGeneric("enc_lkup", function(x) standardGeneric("enc_lkup"))
setMethod("enc_lkup", "BString",
    function(x)
    {
        codec <- codec(x)
        if (is.null(codec)) NULL else codec@enc_lkup
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
### The "BString.read" and "BString.write" functions (NOT exported).
###

BString.read <- function(x, i, imax=integer(0))
{
    XRaw.read(x@data, x@offset + i, x@offset + imax,
                      dec_lkup=dec_lkup(x))
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
### Constructor-like functions and generics
###

.normalize.start <- function(start)
{
    if (!isSingleNumber(start))
        stop("'start' must be a single integer")
    if (!is.integer(start))
        start <- as.integer(start)
    if (start < 1L)
        stop("'start' must be >= 1")
    start
}

.normalize.nchar <- function(start, nchar, src_nchar)
{
    if (!isSingleNumberOrNA(nchar))
        stop("'nchar' must be a single integer or NA")
    if (is.na(nchar)) {
        nchar <- src_nchar - start + 1L
        if (nchar < 0L)
            stop("cannot read a negative number of letters")
        return(nchar)
    }
    if (!is.integer(nchar))
        nchar <- as.integer(nchar)
    if (nchar < 0L)
        stop("cannot read a negative number of letters")
    end <- start + nchar - 1L
    if (end > src_nchar)
        stop("cannot read beyond the end of 'src'")
    nchar
}

.BString.init_with_XRaw <- function(.Object, src, start, nchar, check)
{
    if (check)
        nchar <- .normalize.nchar(start, nchar, length(src))
    slot(.Object, "data", check=check) <- src
    slot(.Object, "offset", check=check) <- start - 1L
    slot(.Object, "length", check=check) <- nchar
    .Object
}

.BString.init_with_character <- function(.Object, src, start, nchar, check, lkup, verbose)
{
    if (length(src) == 0)
        stop("sorry, don't know what to do when 'src' is a character vector of length 0")
    if (length(src) >= 2)
        stop("see ?BStringList when 'src' is a character vector of length >= 2")
    if (check)
        nchar <- .normalize.nchar(start, nchar, nchar(src, type="bytes"))
    src <- substr(src, start, start + nchar - 1L)
    data <- XRaw(nchar, verbose=verbose)
    XRaw.write(data, 1L, nchar, value=src, enc=lkup)
    .BString.init_with_XRaw(.Object, data, 1L, nchar, FALSE)
}

.BString.init_with_BString_copy <- function(.Object, src, start, nchar, check, lkup, verbose)
{
    if (check)
        nchar <- .normalize.nchar(start, nchar, src@length)
    data <- XRaw(nchar, verbose=verbose)
    XRaw.copy(data, src@offset + start, src@offset + start + nchar - 1L, src@data, lkup=lkup)
    .BString.init_with_XRaw(.Object, data, 1L, nchar, FALSE)
}

.BString.init_with_BString <- function(.Object, src, start, nchar, check, copy.data, verbose)
{
    if (copy.data)
        return(.BString.init_with_BString_copy(.Object, src, start, nchar, check, NULL, verbose))
    if (check)
        nchar <- .normalize.nchar(start, nchar, src@length)
    .BString.init_with_XRaw(.Object, src@data, src@offset + start, nchar, FALSE)
}

.BString.get_init_error_msg <- function(.Object, src)
{
    if (class(src) == "BStringViews") {
        if (class(src@subject) == class(.Object))
            return("please use subject() if you are trying to extract the subject of 'src'")
        else
            return("please use BStringViews() when 'src' is a \"BStringViews\" object")
    }
    paste("sorry, don't know what to do when 'src' ",
          "is of class \"", class(src), "\"", sep="")
}

### Because the 'initialize' method for AAString instances is using 'callNextMethod'
### then '.Object' here can be of class "BString" or "AAString".
setMethod("initialize", "BString",
    function(.Object, src, start=1, nchar=NA, check=TRUE, copy.data=FALSE, verbose=FALSE)
    {
        if (check)
            start <- .normalize.start(start)
        if (is.character(src))
            return(.BString.init_with_character(.Object, src, start, nchar, check, NULL, verbose))
        if (class(src) == "XRaw")
            return(.BString.init_with_XRaw(.Object, src, start, nchar, check))
        if (class(src) %in% c("BString", "AAString"))
            return(.BString.init_with_BString(.Object, src, start, nchar, check, copy.data, verbose))
        if (class(.Object) == "BString" && class(src) %in% c("DNAString", "RNAString"))
            return(.BString.init_with_BString_copy(.Object, src, start, nchar, check, dec_lkup(src), verbose))
        stop(.BString.get_init_error_msg(.Object, src))
    }
)

.BString.init_DNAorRNA <- function(.Object, src, start, nchar, check, copy.data, verbose)
{
    lkup <- enc_lkup(.Object) # for source data encoding
    if (is.character(src))
        return(.BString.init_with_character(.Object, src, start, nchar, check, lkup, verbose))
    if (class(src) == "XRaw")
        return(.BString.init_with_XRaw(.Object, src, start, nchar, check))
    if (class(src) == class(.Object))
        return(.BString.init_with_BString(.Object, src, start, nchar, check, copy.data, verbose))
    if (class(src) == "BString")
        return(.BString.init_with_BString_copy(.Object, src, start, nchar, check, lkup, verbose))
    .BString.get_init_error_msg(.Object, src)
}

setMethod("initialize", "DNAString",
    function(.Object, src, start=1, nchar=NA, check=TRUE, copy.data=FALSE, verbose=FALSE)
    {
        if (check)
            start <- .normalize.start(start)
        if (class(src) == "RNAString")
            return(.BString.init_with_BString(.Object, src, start, nchar, check, copy.data, verbose))
        .Object <- .BString.init_DNAorRNA(.Object, src, start, nchar, check, copy.data, verbose)
        if (is.character(.Object))
            stop(.Object)
        .Object
    }
)

setMethod("initialize", "RNAString",
    function(.Object, src, start=1, nchar=NA, check=TRUE, copy.data=FALSE, verbose=FALSE)
    {
        if (check)
            start <- .normalize.start(start)
        if (class(src) == "DNAString")
            return(.BString.init_with_BString(.Object, src, start, nchar, check, copy.data, verbose))
        .Object <- .BString.init_DNAorRNA(.Object, src, start, nchar, check, copy.data, verbose)
        if (is.character(.Object))
            stop(.Object)
        .Object
    }
)

setMethod("initialize", "AAString",
    function(.Object, src, start=1, nchar=NA, check=TRUE, copy.data=FALSE, verbose=FALSE)
    {
        callNextMethod(.Object, src, start=start, nchar=nchar,
                       check=check, copy.data=copy.data, verbose=verbose)
    }
)

### Some wrappers for compatibility with Biostrings 1.
### To test the speed:
###   big <- paste(sample(c('A','C','G','T'), 10^6, replace=TRUE), collapse="")
###   system.time(d <- DNAString(big))

BString <- function(...)
{
    ans <- try(new("BString", ...), silent=TRUE)
    if (is(ans, "try-error")) stop(ans)
    ans
}

DNAString <- function(...)
{
    ans <- try(new("DNAString", ...), silent=TRUE)
    if (is(ans, "try-error")) stop(ans)
    ans
}

RNAString <- function(...)
{
    ans <- try(new("RNAString", ...), silent=TRUE)
    if (is(ans, "try-error")) stop(ans)
    ans
}

AAString <- function(...)
{
    ans <- try(new("AAString", ...), silent=TRUE)
    if (is(ans, "try-error")) stop(ans)
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Standard generic methods

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
### Subsetting (with decoding)

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
        new(class(x), data)
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

BString.comparable <- function(class1, class2)
{
    if (class1 == class2)
        return(TRUE)
    if (class1 %in% c("DNAString", "RNAString")
     && class2 %in% c("DNAString", "RNAString"))
        return(TRUE)
    FALSE
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
        class1 <- class(e1)
        class2 <- class(e2)
        if (!BString.comparable(class1, class2))
            stop("comparison between a \"", class1, "\" instance ",
                 "and a \"", class2, "\" instance ",
                 "is not supported")
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
    slot(x, "offset", check=FALSE) <- x@offset + start - 1L
    slot(x, "length", check=FALSE) <- end - shift
    x
}

