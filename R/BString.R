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
# The "codec", "enc_hash" and "dec_hash" new generics (NOT exported)

setGeneric("codec", function(x) standardGeneric("codec"))
setMethod("codec", "BString", function(x) NULL)
setMethod("codec", "DNAString", function(x) DNA_STRING_CODEC)
setMethod("codec", "RNAString", function(x) RNA_STRING_CODEC)

setGeneric("dec_hash", function(x) standardGeneric("dec_hash"))
setMethod("dec_hash", "BString",
    function(x)
    {
        codec <- codec(x)
        if (is.null(codec)) NULL else codec@dec_hash
    }
)

setGeneric("enc_hash", function(x) standardGeneric("enc_hash"))
setMethod("enc_hash", "BString",
    function(x)
    {
        codec <- codec(x)
        if (is.null(codec)) NULL else codec@enc_hash
    }
)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# The "alphabet" generic is defined in BStringCodec.R

setMethod("alphabet", "DNAString", function(x) DNA_ALPHABET)
setMethod("alphabet", "RNAString", function(x) RNA_ALPHABET)
setMethod("alphabet", "AAString", function(x) AA_ALPHABET)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# The "BString.read" and "BString.write" functions (NOT exported)

BString.read <- function(x, i, imax=integer(0))
{
    CharBuffer.read(x@data, x@offset + i, x@offset + imax,
                    dec_hash=dec_hash(x))
}

# Only used at initialization time!
BString.write <- function(x, i, imax=integer(0), value)
{
    CharBuffer.write(x@data, x@offset + i, x@offset + imax, value=value,
                     enc_hash=enc_hash(x))
    x
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Accessor methods

# Returns a character-string
setGeneric("letter", function(x, i) standardGeneric("letter"))
setMethod("letter", "BString",
    function(x, i)
    {
        if (!isTRUE(all(i >= 1)) || !isTRUE(all(i <= x@length))) # NA-proof
            stop("subscript out of bounds")
        BString.read(x, i)
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
        if (class(.Object) == "BString" && class(src) %in% c("DNAString", "RNAString"))
            return(BString.init_with_BString_copy(.Object, src, dec_hash(src), verbose))
        stop(BString.get_init_error_msg(.Object, src))
    }
)

BString.init_DNA_or_RNA <- function(.Object, src, copy.data, verbose)
{
    hash <- enc_hash(.Object) # for source data encoding
    if (is.character(src))
        return(BString.init_with_character(.Object, src, hash, verbose))
    if (class(src) == "CharBuffer")
        return(BString.init_with_CharBuffer(.Object, src))
    if (class(src) == class(.Object))
        return(BString.init_with_BString(.Object, src, copy.data, verbose))
    if (class(src) == "BString")
        return(BString.init_with_BString_copy(.Object, src, hash, verbose))
    BString.get_init_error_msg(.Object, src)
}

setMethod("initialize", "DNAString",
    function(.Object, src, copy.data=FALSE, verbose=FALSE)
    {
        if (class(src) == "RNAString")
            return(BString.init_with_BString(.Object, src, copy.data, verbose))
        .Object <- BString.init_DNA_or_RNA(.Object, src, copy.data, verbose)
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
        .Object <- BString.init_DNA_or_RNA(.Object, src, copy.data, verbose)
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
#   dna <- Hsapiens$chr1
#   dna != Hsapiens$chr1 # FALSE
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
setMethod("as.character", "BString", function(x) BString.read(x, 1, x@length))
setMethod("toString", "BString", function(x) as.character(x))
setMethod("nchar", "BString", function(x, type) x@length)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Other functions and generics

# The "BString.substring" function is very fast because it does not copy
# the string data. Return a BString object (not vectorized).
# 'first' and 'last' must be single integers verifying:
#   1 <= first <= last <= length(x)
# WARNING: This function is voluntarly unsafe (it doesn't check its
# arguments) because we want it to be the fastest possible!
BString.substring <- function(x, first, last)
{
    one <- as.integer(1)
    x@offset <- x@offset + first - one
    x@length <- last - first + one
    x
}

# The public (and safe) version of BString.substring(). Not vectorized.
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
        BString.substring(x, as.integer(first), as.integer(last))
    }
)



