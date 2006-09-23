# ===========================================================================
# The BString class
# ---------------------------------------------------------------------------

setClass(
    "BString",
    representation(
        data="CharBuffer",      # contains the string data
        offset="integer",       # a single integer
        length="integer"        # a single integer
    )
)

# 3 direct "BString" derivations (no additional slot)
setClass("DNAString", representation("BString"))
setClass("RNAString", representation("BString"))
setClass("AAString", representation("BString"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# The "readChars" and "writeChars" new generics

setGeneric(
    "readChars",
    function(x, i, imax=integer(0)) standardGeneric("readChars")
)
setMethod("readChars", "BString",
    function(x, i, imax)
    {
        CharBuffer.read(x@data, x@offset + i, x@offset + imax)
    }
)
setMethod("readChars", "DNAString",
    function(x, i, imax)
    {
        dec_hash = DNA_STRING_CODEC@dec_hash
        CharBuffer.read(x@data, x@offset + i, x@offset + imax, dec=dec_hash)
    }
)
setMethod("readChars", "RNAString",
    function(x, i, imax)
    {
        dec_hash = RNA_STRING_CODEC@dec_hash
        CharBuffer.read(x@data, x@offset + i, x@offset + imax, dec=dec_hash)
    }
)

# Only used at initialization!
setGeneric(
    "writeChars",
    function(x, i, imax=integer(0), value) standardGeneric("writeChars")
)
setMethod("writeChars", "BString",
    function(x, i, imax, value)
    {
        CharBuffer.write(x@data, x@offset + i, x@offset + imax, value=value)
        x
    }
)
setMethod("writeChars", "DNAString",
    function(x, i, imax, value)
    {
        enc_hash = DNA_STRING_CODEC@enc_hash
        CharBuffer.write(x@data, x@offset + i, x@offset + imax, value=value,
                         enc=enc_hash)
        x
    }
)
setMethod("writeChars", "RNAString",
    function(x, i, imax, value)
    {
        enc_hash = RNA_STRING_CODEC@enc_hash
        CharBuffer.write(x@data, x@offset + i, x@offset + imax, value=value,
                         enc=enc_hash)
        x
    }
)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Accessor methods

# Returns a character-string
setGeneric("letter", function(x, i) standardGeneric("letter"))
setMethod("letter", "BString",
    function(x, i)
    {
        if (!isTRUE(all(i >= 1)) || !isTRUE(all(i <= x@length))) # NA-proof
            stop("subscript out of bounds")
        readChars(x, i)
    }
)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Constructor-like functions and generics

BString.init_with_CharBuffer <- function(.Object, src)
{
    .Object@data <- src
    .Object@offset <- as.integer(0)
    .Object@length <- length(src)
    .Object
}

BString.init_with_character <- function(.Object, src, hash=NULL, verbose=FALSE)
{
    if (length(src) == 0)
        stop("sorry, don't know what to do when 'src' is a character vector of length 0")
    if (length(src) >= 2)
        stop("please use BStringViews() when 'src' is a character vector of length >= 2")
    length <- nchar(src)
    data <- CharBuffer(length, verbose)
    CharBuffer.write(data, 1, length, value=src, enc=hash)
    BString.init_with_CharBuffer(.Object, data)
}

BString.init_with_BString_copy <- function(.Object, src, hash=NULL, verbose=FALSE)
{
    length <- src@length
    data <- CharBuffer(length, verbose)
    CharBuffer.copy(data, src@offset + 1, src@offset + length, src@data, hash=hash)
    BString.init_with_CharBuffer(.Object, data)
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

BString.getInitErrorMsg <- function(.Object, src)
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

# Because the 'initialize' method for "AAString" objects is using 'callNextMethod'
# then '.Object' here can be of class "BString" or "AAString".
setMethod("initialize", "BString",
    function(.Object, src, copy.data=FALSE, verbose=FALSE)
    {
        if (is.character(src))
            return(BString.init_with_character(.Object, src, , verbose))
        if (class(src) == "CharBuffer")
            return(BString.init_with_CharBuffer(.Object, src))
        if (class(src) %in% c("BString", "AAString"))
            return(BString.init_with_BString(.Object, src, copy.data, verbose))
        if (class(.Object) == "BString") {
            if (class(src) == "DNAString") {
                hash <- DNA_STRING_CODEC@dec_hash # for source data decoding
                return(BString.init_with_BString_copy(.Object, src, hash, verbose))
            }
            if (class(src) == "RNAString") {
                hash <- RNA_STRING_CODEC@dec_hash # for source data decoding
                return(BString.init_with_BString_copy(.Object, src, hash, verbose))
            }
        }
        stop(BString.getInitErrorMsg(.Object, src))
    }
)

.initEncodedBString <- function(.Object, src, hash, copy.data, verbose)
{
    if (is.character(src))
        return(BString.init_with_character(.Object, src, hash, verbose))
    if (class(src) == "CharBuffer")
        return(BString.init_with_CharBuffer(.Object, src))
    if (class(src) == class(.Object))
        return(BString.init_with_BString(.Object, src, copy.data, verbose))
    if (class(src) == "BString")
        return(BString.init_with_BString_copy(.Object, src, hash, verbose))
    BString.getInitErrorMsg(.Object, src)
}

setMethod("initialize", "DNAString",
    function(.Object, src, copy.data=FALSE, verbose=FALSE)
    {
        if (class(src) == "RNAString")
            return(BString.init_with_BString(.Object, src, copy.data, verbose))
        hash <- DNA_STRING_CODEC@enc_hash # for source data encoding
        .Object <- .initEncodedBString(.Object, src, hash, copy.data, verbose)
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
        hash <- RNA_STRING_CODEC@enc_hash # for source data encoding
        .Object <- .initEncodedBString(.Object, src, hash, copy.data, verbose)
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

# Some wrappers for compatibility with Biostrings 1.
# To test the speed:
#   big <- paste(sample(c('A','C','G','T'), 10^6, replace=TRUE), collapse="")
#   system.time(d <- DNAString(big))

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


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Standard generic methods

# Helper function used by the show() method
bsSnippet <- function(x, snippetW)
{
    if (snippetW < 7)
        snippetW <- 7
    lx <- x@length
    if (lx <= snippetW) {
        toString(x)
    } else {
        w1 <- (snippetW - 2) %/% 2
        w2 <- (snippetW - 3) %/% 2
        paste(readChars(x, 1, w1),
              "...",
              readChars(x, lx - w2 + 1, lx),
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
        cat("\nValue:", bsSnippet(object, 72))
        cat("\n")
    }
)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Subsetting (with decoding)

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
        data <- CharBuffer(length(i))
        CharBuffer.copy(data, x@offset + i, src=x@data)
        # class(x) can be "BString" or one of its derivations ("DNAString",
        # "RNAString" or "AAString").
        new(class(x), data)
    }
)

# The only reason for defining the replacement version of the "[" operator
# is to let the user know that he can't use it:
#   bs <- BString("AbnbIU")
#   bs[2] <- "X" # provokes an error
# If we don't define it, then the user can type the above and believe that
# it actually did something but it didn't.
setReplaceMethod("[", "BString",
    function(x, i, j,..., value)
    {
        stop(paste("attempt to modify the value of a \"",
                   class(x), "\" object", sep=""))
    }
)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Equality

# We want:
#   BString("ab") == "ab" # TRUE
#   BString("ab") == "" # Error ("" can't be converted to a BString)
#   DNAString("TG") == "TG" # TRUE
#   "TT" == DNAString("TG") # FALSE
#   "TGG" == DNAString("TG") # FALSE
#   DNAString("TG") == RNAString("UG") # TRUE!!!
#   library(BSgenome.Hsapiens.UCSC.hg18)
#   dna <- Hsapiens$chr1[[1]]
#   dna != Hsapiens$chr1[[1]] # FALSE
#   dnav <- views(dna, 1:7, 101:107)
#   dnav[[1]] == dnav[[7]] # TRUE
#   dnav <- views(dna, 1:7, (length(dna)-6):length(dna))
# This is fast:
#   dnav[[1]] == dnav[[7]] # FALSE
# But this would have killed your machine:
#   s1 <- toString(dnav[[1]])
#   s7 <- toString(dnav[[7]])
#   s1 == s7

BString.equal <- function(x, y)
{
    if (class(y) != class(x))
        y <- new(class(x), y)
    if (x@length != y@length)
        return(FALSE)
    one <- as.integer(1)
    ans <- !CharBuffer.compare(x@data, x@offset + one, y@data, y@offset + one, x@length)
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

# These methods are called if at least one side of the "==" (or "!=")
# operator is a "BString" object AND the other side is NOT a "BStringViews"
# object.
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


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("as.character", "BString", function(x) readChars(x, 1, x@length))
setMethod("toString", "BString", function(x) as.character(x))
setMethod("nchar", "BString", function(x, type) x@length)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Other functions and generics

# The "bsSubstr" function is very fast because it does not copy the string
# data. Return a BString object (not vectorized).
# 'first' and 'last' must be single integers verifying:
#   1 <= first <= last <= length(x)
# WARNING: This function is voluntarly unsafe (it doesn't check its
# arguments) because we want it to be the fastest possible!
bsSubstr <- function(x, first, last)
{
    one <- as.integer(1)
    x@offset <- x@offset + first - one
    x@length <- last - first + one
    x
}

# The public (and safe) version of bsSubstr(). Not vectorized.
# We deliberately choose the "NA trick" over defaulting 'first' and 'last'
# to '1' and 'length(x)' because we want to be consistent with what the
# views() function does.
setGeneric(
    "subBString", function(x, first=NA, last=NA) standardGeneric("subBString")
)
setMethod("subBString", "BString",
    function(x, first, last)
    {
        if (!isLooseNumeric(first) || !isLooseNumeric(last))
            stop("'first' and 'last' must be numerics")
        if (length(first) != 1 || length(last) != 1)
            stop("'first' and 'last' must be single numerics")
        if (is.na(first))
            first <- 1
        if (is.na(last))
            last <- x@length
        # This is NA-proof (well, 'first' and 'last' can't be NAs anymore...)
        if (!isTRUE(1 <= first && first <= last && last <= length(x)))
            stop("'first' and 'last' must verify '1 <= first <= last <= length(x)'")
        bsSubstr(x, as.integer(first), as.integer(last))
    }
)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Helper functions for view manipulation

# The 2 functions below convert a given view on a BString object into a
# a character-string.
# They are used as helper functions to display a BStringViews object.
# Both assume that 'first' <= 'last' (so they don't check it) and
# padd the result with spaces to produce the "margin effect"
# if 'first' or 'last' are out of limits.

# nchar(bsView(x, first, last)) is always last-first+1
bsView <- function(x, first, last)
{
    lx <- length(x)
    if (last < 1 || first > lx)
            return(format("", width=last-first+1))
    Lmargin <- ""
    if (first < 1) {
        Lmargin <- format("", width=1-first)
        first <- 1
    }
    Rmargin <- ""
    if (last > lx) {
        Rmargin <- format("", width=last-lx)
        last <- lx
    }
    paste(Lmargin, readChars(x, first, last), Rmargin, sep="")
}

# nchar(bsViewSnippet(x, first, last, snippetW)) is <= snippetW
bsViewSnippet <- function(x, first, last, snippetW)
{
    if (snippetW < 7)
        snippetW <- 7
    width <- last - first + 1
    if (width <= snippetW) {
        bsView(x, first, last)
    } else {
        w1 <- (snippetW - 2) %/% 2
        w2 <- (snippetW - 3) %/% 2
        paste(bsView(x, first, first+w1-1),
              "...",
              bsView(x, last-w2+1, last), sep="")
    }
}

# Assume that 'first1', 'last1', 'first2', 'last2' are single integers
# and that first1 <= last1 and first2 <= last2.
bsIdenticalViews <- function(x1, first1, last1, x2, first2, last2)
{
    one <- as.integer(1)
    w1 <- last1 - first1 + one
    w2 <- last2 - first2 + one
    if (w1 != w2)
        return(FALSE)

    lx1 <- length(x1)
    isBlank1 <- last1 < one || first1 > lx1
    lx2 <- length(x2)
    isBlank2 <- last2 < one || first2 > lx2
    if (isBlank1 && isBlank2)
        return(TRUE)
    if (isBlank1 || isBlank2)
        return(FALSE)

    # Left margin
    LmarginSize1 <- first1 < one
    LmarginSize2 <- first2 < one
    if (LmarginSize1 != LmarginSize2)
        return(FALSE)
    if (LmarginSize1) {
        # Both views have a left margin
        if (first1 != first2)
            return(FALSE)
        first1 <- one
        first2 <- one
    }

    # Right margin
    RmarginSize1 <- last1 > lx1
    RmarginSize2 <- last2 > lx2
    if (RmarginSize1 != RmarginSize2)
        return(FALSE)
    if (RmarginSize1) {
        # Both views have a right margin
        if (last1 - lx1 != last2 - lx2)
            return(FALSE)
        last1 <- lx1
        last2 <- lx2
    }

    # At this point, we can trust that 1 <= first1 <= last1 <= lx1
    # and that 1 <= first2 <= last2 <= lx2 so we can call unsafe
    # function bsSubstr() with no fear...
    bsSubstr(x1, first1, last1) == bsSubstr(x2, first2, last2)
}
