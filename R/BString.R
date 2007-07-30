### =========================================================================
### The BString class
### -------------------------------------------------------------------------

setClass(
    "BString",
    representation(
        data="XRaw",        # contains the string data
        offset="integer",   # a single integer
        length="integer"    # a single integer
    )
)

### 3 direct "BString" derivations (no additional slot)
setClass("DNAString", representation("BString"))
setClass("RNAString", representation("BString"))
setClass("AAString", representation("BString"))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "codec", "enc_lkup" and "dec_lkup" new generics (NOT exported).
###

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
BString.write <- function(x, i, imax=integer(0), value)
{
    XRaw.write(x@data, x@offset + i, x@offset + imax, value=value,
                       enc_lkup=enc_lkup(x))
    x
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor-like functions and generics

BString.init_with_XRaw <- function(.Object, src)
{
    .Object@data <- src
    .Object@offset <- as.integer(0)
    .Object@length <- length(src)
    .Object
}

BString.init_with_character <- function(.Object, src, lkup=NULL, verbose=FALSE)
{
    if (length(src) == 0)
        stop("sorry, don't know what to do when 'src' is a character vector of length 0")
    if (length(src) >= 2)
        stop("please use BStringViews() when 'src' is a character vector of length >= 2")
    length <- nchar(src)
    data <- XRaw(length, verbose)
    XRaw.write(data, 1, length, value=src, enc=lkup)
    BString.init_with_XRaw(.Object, data)
}

BString.init_with_BString_copy <- function(.Object, src, lkup=NULL, verbose=FALSE)
{
    length <- src@length
    data <- XRaw(length, verbose)
    XRaw.copy(data, src@offset + 1, src@offset + length, src@data, lkup=lkup)
    BString.init_with_XRaw(.Object, data)
}

BString.init_with_BString <- function(.Object, src, copy.data=FALSE, verbose=FALSE)
{
    if (copy.data)
        return(BString.init_with_BString_copy(.Object, src, , verbose))
    .Object@data <- src@data
    .Object@offset <- src@offset
    .Object@length <- src@length
    .Object
}

BString.get_init_error_msg <- function(.Object, src)
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

### Because the 'initialize' method for "AAString" objects is using 'callNextMethod'
### then '.Object' here can be of class "BString" or "AAString".
setMethod("initialize", "BString",
    function(.Object, src, copy.data=FALSE, verbose=FALSE)
    {
        if (is.character(src))
            return(BString.init_with_character(.Object, src, , verbose))
        if (class(src) == "XRaw")
            return(BString.init_with_XRaw(.Object, src))
        if (class(src) %in% c("BString", "AAString"))
            return(BString.init_with_BString(.Object, src, copy.data, verbose))
        if (class(.Object) == "BString" && class(src) %in% c("DNAString", "RNAString"))
            return(BString.init_with_BString_copy(.Object, src, dec_lkup(src), verbose))
        stop(BString.get_init_error_msg(.Object, src))
    }
)

BString.init_DNAorRNA <- function(.Object, src, copy.data, verbose)
{
    lkup <- enc_lkup(.Object) # for source data encoding
    if (is.character(src))
        return(BString.init_with_character(.Object, src, lkup, verbose))
    if (class(src) == "XRaw")
        return(BString.init_with_XRaw(.Object, src))
    if (class(src) == class(.Object))
        return(BString.init_with_BString(.Object, src, copy.data, verbose))
    if (class(src) == "BString")
        return(BString.init_with_BString_copy(.Object, src, lkup, verbose))
    BString.get_init_error_msg(.Object, src)
}

setMethod("initialize", "DNAString",
    function(.Object, src, copy.data=FALSE, verbose=FALSE)
    {
        if (class(src) == "RNAString")
            return(BString.init_with_BString(.Object, src, copy.data, verbose))
        .Object <- BString.init_DNAorRNA(.Object, src, copy.data, verbose)
        if (is.character(.Object))
            stop(.Object)
        .Object
    }
)

setMethod("initialize", "RNAString",
    function(.Object, src, copy.data=FALSE, verbose=FALSE)
    {
        if (class(src) == "DNAString")
            return(BString.init_with_BString(.Object, src, copy.data, verbose))
        .Object <- BString.init_DNAorRNA(.Object, src, copy.data, verbose)
        if (is.character(.Object))
            stop(.Object)
        .Object
    }
)

setMethod("initialize", "AAString",
    function(.Object, src, copy.data=FALSE, verbose=FALSE)
    {
        callNextMethod(.Object, src, copy.data, verbose)
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
        cat("  ", lo, "-letter \"", class(object), "\" object", sep="")
        #if (!is.null(object@codec))
        #    cat(" with alphabet:", toString(object@codec@letters))
        cat("\nValue:", BString.get_snippet(object, 72))
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
        if (any(i < 1) || any(i > length(x)))
            stop("subscript out of bounds")
        data <- XRaw(length(i))
        XRaw.copy(data, x@offset + i, src=x@data)
        ## class(x) can be "BString" or one of its derivations ("DNAString",
        ## "RNAString" or "AAString").
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
                   class(x), "\" object", sep=""))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Equality.
###

### We want:
###   BString("ab") == "ab" # TRUE
###   BString("ab") == "" # Error ("" can't be converted to a BString)
###   DNAString("TG") == "TG" # TRUE
###   "TT" == DNAString("TG") # FALSE
###   "TGG" == DNAString("TG") # FALSE
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

BString.equal <- function(x, y)
{
    if (class(y) != class(x))
        y <- new(class(x), y)
    if (x@length != y@length)
        return(FALSE)
    one <- as.integer(1)
    ans <- !XRaw.compare(x@data, x@offset + one, y@data, y@offset + one, x@length)
    as.logical(ans)
}

setMethod("==", signature(e1="BString", e2="BString"),
    function(e1, e2)
    {
        if (class(e1) == class(e2))
            return(BString.equal(e1, e2))
        if (class(e1) == "BString") # then e2 is a DNAString, RNAString or AAString
            return(BString.equal(e2, e1))
        if (class(e2) == "BString") # then e1 is a DNAString, RNAString or AAString
            return(BString.equal(e1, e2))
        if (class(e1) != "AAString" && class(e2) != "AAString")
            return(BString.equal(e1, e2))
        stop(paste("sorry, don't know how to compare ",
                   "a \"", class(e1), "\" object with ",
                   "a \"", class(e2), "\" object", sep=""))
    }
)
setMethod("!=", signature(e1="BString", e2="BString"),
    function(e1, e2) !(e1 == e2)
)

### These methods are called if at least one side of the "==" (or "!=")
### operator is a "BString" object AND the other side is NOT a "BStringViews"
### object.
setMethod("==", signature(e1="BString"),
    function(e1, e2) BString.equal(e1, e2)
)
setMethod("!=", signature(e1="BString"),
    function(e1, e2) !(e1 == e2)
)

setMethod("==", signature(e2="BString"),
    function(e1, e2) BString.equal(e2, e1)
)
setMethod("!=", signature(e2="BString"),
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
    one <- as.integer(1)
    x@offset <- x@offset + start - one
    x@length <- end - start + one
    x
}

