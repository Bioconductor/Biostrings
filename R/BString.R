# ===========================================================================
# The BString class
# ---------------------------------------------------------------------------

setClass(
    "BString",
    representation(
        data="bbuf",            # byte buffer (contains the string data)
        offset="integer",       # a single integer
        length="integer"        # a single integer
    )
)

# 2 direct extensions of the BString class
# No additional slot!
setClass("DNAString", representation("BString"))
setClass("RNAString", representation("BString"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("length", "BString", function(x) x@length)
setMethod("nchar", "BString", function(x, type) x@length)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# The "readChars" and "writeChars" new generics

setGeneric(
    "readChars",
    function(x, i, imax=integer(0)) standardGeneric("readChars")
)
setMethod("readChars", "BString",
    function(x, i, imax)
    {
        bbReadChars(x@data, i + x@offset, imax + x@offset)
    }
)
setMethod("readChars", "DNAString",
    function(x, i, imax)
    {
        dec_hash = DNA_STRING_CODEC@dec_hash
        bbReadChars(x@data, i + x@offset, imax + x@offset, dec=dec_hash)
    }
)
setMethod("readChars", "RNAString",
    function(x, i, imax)
    {
        dec_hash = RNA_STRING_CODEC@dec_hash
        bbReadChars(x@data, i + x@offset, imax + x@offset, dec=dec_hash)
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
        bbWriteChars(x@data, i + x@offset, imax + x@offset, value=value)
        x
    }
)
setMethod("writeChars", "DNAString",
    function(x, i, imax, value)
    {
        enc_hash = DNA_STRING_CODEC@enc_hash
        bbWriteChars(x@data, i + x@offset, imax + x@offset, value=value,
                     enc=enc_hash)
        x
    }
)
setMethod("writeChars", "RNAString",
    function(x, i, imax, value)
    {
        enc_hash = RNA_STRING_CODEC@enc_hash
        bbWriteChars(x@data, i + x@offset, imax + x@offset, value=value,
                     enc=enc_hash)
        x
    }
)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("as.character", "BString", function(x) readChars(x, 1, x@length))
setMethod("toString", "BString", function(x) as.character(x))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Initialization

# Must work at least with 'src' being one of the following:
#   - a single non-empty string (character vector of length 1)
#   - a "bbuf" object
#   - a "BString", "DNAString" or "RNAString" object
setMethod("initialize", "BString",
    function(.Object, src)
    {
        # class(.Object) can be "BString", "DNAString" or "RNAString"!
        if (class(src) == class(.Object))
            return(src)
        .Object@offset <- as.integer(0) # Set me BEFORE calling writeChars()
        if (class(src) == "bbuf") {
            length <- length(src)
            .Object@data <- src
        } else {
            if (is.character(src)) {
                if (length(src) != 1)
                    stop("use BStringViews() when 'src' is a character vector of length != 1")
            } else {
                if (class(src) == "BStringViews") {
                    if (class(src@subject) == class(.Object))
                        stop("use subject() if you want to extract the subject of 'src'")
                    else
                        stop("use BStringViews() when 'src' is a \"BStringViews\" object")
                }
                src <- toString(src)
            }
            length <- nchar(src)
            .Object@data <- bbuf(length)
            writeChars(.Object, 1, length, value=src)
        }
        .Object@length <- length
        .Object
    }
)
BString <- function(...)
{
    new("BString", ...)
}

# Uses global variable DNA_STRING_CODEC to encode source string.
setMethod("initialize", "DNAString",
    function(.Object, src)
    {
        if (class(src) == "RNAString") {
            .Object@data <- src@data
            .Object@offset <- src@offset
            .Object@length <- src@length
            return(.Object)
        }
        callNextMethod(.Object, src)
    }
)

# Uses global variable DNA_STRING_CODEC to encode source string.
setMethod("initialize", "RNAString",
    function(.Object, src)
    {
        if (class(src) == "DNAString") {
            .Object@data <- src@data
            .Object@offset <- src@offset
            .Object@length <- src@length
            return(.Object)
        }
        callNextMethod(.Object, src)
    }
)

# Some wrappers for compatibility with BStrings 1.4.x (BioC 1.7).
# To test the speed:
#   big <- paste(sample(c('A','C','G','T'), 10^6, replace=TRUE), collapse="")
#   system.time(d <- DNAString(big))

DNAString <- function(...)
{
    new("DNAString", ...)
}
RNAString <- function(...)
{
    new("RNAString", ...)
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Subsetting (with decoding)

setMethod("[", "BString",
    function(x, i, j, ..., drop)
    {
        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")
        if (missing(i))
            return(x)
        if (any(i < 1) || any(i > length(x)))
            stop("subscript out of bounds")
        data <- bbuf(length(i))
        bbCopy(data, i + x@offset, src=x@data)
        # class(x) can be "BString", "DNAString" or "RNAString"
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

isLooseNumeric <- function(x)
{
    return(is.numeric(x) || (!is.null(x) && all(is.na(x))))
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
# Display

bsSnippet <- function(x, snippetW)
{
    if (snippetW < 7)
        snippetW <- 7
    l <- length(x)
    if (l <= snippetW) {
        toString(x)
    } else {
        w1 <- (snippetW - 2) %/% 2
        w2 <- (snippetW - 3) %/% 2
        paste(readChars(x, 1, w1),
              "...",
              readChars(x, l - w2 + 1, l),
              sep="")
    }
}

setMethod("show", "BString",
    function(object)
    {
        l <- length(object)
        cat("  ", l, "-letter \"", class(object), "\" object", sep="")
        #if (!is.null(object@codec))
        #    cat(" with alphabet:", toString(object@codec@letters))
        cat("\nValue:", bsSnippet(object, 72))
        cat("\n")
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
#   library(HsapiensGenome)
#   dna <- hsa$chr1[[1]]
#   dna != hsa$chr1[[1]] # FALSE
#   dnav <- views(dna, 1:7, 101:107)
#   dnav[[1]] == dnav[[7]] # TRUE
#   dnav <- views(dna, 1:7, (length(dna)-6):length(dna))
# This is fast:
#   dnav[[1]] == dnav[[7]] # FALSE
# But this would have killed your machine:
#   s1 <- toString(dnav[[1]])
#   s7 <- toString(dnav[[7]])
#   s1 == s7

.different <- function(x, y)
{
    if (class(y) != class(x))
        y <- new(class(x), y)
    if (x@length != y@length)
        return(TRUE)
    one <- as.integer(1)
    ans <- bbCompare(x@data, one + x@offset, y@data, one + y@offset, x@length)
    as.logical(ans)
}

# These methods are called if at least one side of the "==" (or "!=")
# operator is a "BString" object AND the other side is NOT a "BStringViews"
# object.
setMethod("==", signature(e1="BString"),
    function(e1, e2) !.different(e1, e2)
)
setMethod("==", signature(e2="BString"),
    function(e1, e2) !.different(e2, e1)
)

setMethod("!=", signature(e1="BString"),
    function(e1, e2) .different(e1, e2)
)
setMethod("!=", signature(e2="BString"),
    function(e1, e2) .different(e2, e1)
)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# These 2 functions convert a given view on a BString object into a
# a single string (character vector of length 1).
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


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
