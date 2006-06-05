# ===========================================================================
# The BStringViews class
# ---------------------------------------------------------------------------
# A BStringViews object contains a set of views
# on the same BString object, the subject string.

setClassUnion("BStringUnion", c("BString", "DNAString", "RNAString"))

# See the initialization section below for the integrity checking
# of a BStringViews object.
setClass(
    "BStringViews",
    representation(
        subject="BStringUnion",
        first="integer",
        last="integer",
        desc="character"   # store per-view comment (e.g. from FASTA file)
    )
)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Accessor methods

setGeneric("subject", function(x) standardGeneric("subject"))
setMethod("subject", "BStringViews", function(x) x@subject)

# Names for slots 'first' and 'last' are those of the arguments of the
# substring() function.
# Because the start() and stop() functions are already defined as R standard
# functions, using names 'start' and 'stop', like in the substr() function,
# is a problem if we want to be able to name the following accessor
# functions like the slot they access to.
setGeneric("first", function(x) standardGeneric("first"))
setMethod("first", "BStringViews", function(x) x@first)

setGeneric("last", function(x) standardGeneric("last"))
setMethod("last", "BStringViews", function(x) x@last)

# We choose to call this method 'width' and not 'length' because
# we want to define the length of a BStringViews object as the number
# of views contained in it.
# Another option was to call it 'nchar' but the width of a view is not
# necesarily equal to the number of letters that it contains (this happens
# when the view is out of limits).
setGeneric("width", function(x) standardGeneric("width"))
setMethod("width", "BStringViews", function(x) x@last - x@first + 1)

setGeneric("desc", function(x) standardGeneric("desc"))
setMethod("desc", "BStringViews", function(x) x@desc)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Constructor-like functions and generics

# WARNING: This function is unsafe! (it doesn't check its arguments)
# Only 2 valid ways to use it:
#   new("BStringViews", subject)
#   new("BStringViews", subject, first, last)
# where 'subject' is a BString (or derived) object,
# and 'first' and 'last' are integer vectors of the same length
# such that 'first <= last'.
setMethod("initialize", "BStringViews",
    function(.Object, subject, first, last)
    {
        .Object@subject <- subject
        if (!missing(first)) {
            .Object@first <- first
            .Object@last <- last
        }
        .Object
    }
)

# The 2 functions above share the following properties:
#   - They are exported (and safe).
#   - First argument is 'subject'. It must be a character vector or a BString
#     (or derived) object.
#   - Passing something else to 'subject' provokes an error.
#   - They return a BStringViews object whose 'subject' slot is the object
#     passed in the 'subject' argument.

# Typical use:
#   dna <- DNAString("AA-CC-GG-TT")
# Just one view:
#   dnav1 <- views(dna, 2, 7)
# 9 views, 3 are out of limits:
#   dnav2 <- views(dna, 6:-2, 6:14)
# 5 out of limits views, all have a width of 6:
#   dnav3 <- views(dna, -5:-1, 0:4)
# Same as doing views(dna, 1, length(dna)):
#   dnav4 <- views(dna)
# A BStringViews object with no view:
#   dnav5 <- views(dna, integer(0), integer(0))
views <- function(subject, first=NA, last=NA)
{
    if (class(subject) == "character")
        subject <- BString(subject)
    ans <- new("BStringViews", subject)
    # Integrity checking
    if (!isLooseNumeric(first) || !isLooseNumeric(last))
        stop("'first' and 'last' must be numerics")
    #if (length(first) != length(last))
    #    stop("'first' and 'last' must have the same length")
    if (!is.integer(first))
        first <- as.integer(first)
    first[is.na(first)] <- as.integer(1)
    if (!is.integer(last))
        last <- as.integer(last)
    last[is.na(last)] <- subject@length
    if (length(first) < length(last))
        first <- recycleVector(first, length(last))
    else if (length(last) < length(first))
        last <- recycleVector(last, length(first))
    # The NA-proof version of 'if (any(last < first))'
    if (!isTRUE(all(first <= last)))
        stop("'first' and 'last' must verify 'first <= last'")
    ans@first <- first
    ans@last <- last
    ans
}

# 'width' is the vector of view widths.
# 'gapwidth' is the vector of inter-view widths (recycled).
adjacentViews <- function(subject, width, gapwidth=0)
{
    if (class(subject) == "character")
        subject <- BString(subject)
    ans <- new("BStringViews", subject)
    if (!is.numeric(width) || !isTRUE(all(width >= 1))) # NA-proof
        stop("'width' must be numerics >= 1")
    if (!is.numeric(gapwidth) || !isTRUE(all(gapwidth >= 0))) # NA-proof
        stop("'gapwidth' must be numerics >= 0")
    lw <- length(width)
    if (lw == 0)
        return(ans)
    if (!is.integer(width))
        width <- as.integer(width)
    lg <- length(gapwidth)
    if (!is.integer(gapwidth))
        gapwidth <- as.integer(gapwidth)
    first <- integer(lw)
    last <- integer(lw)
    one <- as.integer(1)
    first[1] <- one
    last[1] <- width[1]
    if (lw >= 2) {
        j <- 1
        for (i in 2:lw) {
            first[i] <- last[i-1] + one + gapwidth[j]
            last[i] <- first[i] + width[i] - one
            if (j < lg) j <- j + 1 else j <- 1
        }
    }
    ans@first <- first
    ans@last <- last
    ans
}


# The BStringViews() generic function
# -----------------------------------
# 'src' should typically be a character vector but the function will also
# work for other kind of input like numeric or even logical vectors.
# 'subjectClass' must be "BString", "DNAString" or "RNAString".
#
# Benchmarks:
#   alphabet <- strsplit("-TGCANBDHKMRSVWY", NULL)[[1]]
#   n <- 40000
#   src <- sapply(1:n, function(i) {paste(sample(alphabet, 250, replace=TRUE), collapse="")})
#   v <- BStringViews(src, "DNAString")
# Comparing BStringViews() speed vs "old" vectorized DNAString() speed:
#       n  BStringViews  "old" DNAString
#   -----  ------------  ---------------
#    5000        0.26 s           4.15 s
#   10000        0.51 s          16.29 s
#   20000        0.99 s          64.85 s
#   40000        1.69 s         488.43 s
# The quadratic behaviour of "old" DNAString() was first reported
# by Wolfgang.
BStringViews <- function(src, subjectClass, sep="")
{
    if (!is.character(sep))
        sep <- toString(sep)
    collapsed <- paste(src, collapse=sep)
    subject <- new(subjectClass, collapsed)
    adjacentViews(subject, nchar(src), nchar(sep))
}

setGeneric(
    "BStringViews",
    function(src, subjectClass, sep="") standardGeneric("BStringViews")
)

# Only FASTA files are supported for now.
# Typical use:
#   srcpath <- system.file("Exfiles", "someORF.fsa", package="Biostrings")
#   f <- file(srcpath)
#   v <- BStringViews(f, "DNAString")
#   close(f)
setMethod("BStringViews", "file",
    function(src, subjectClass, sep)
    {
        fasta <- readFASTA(src)
        src <- sapply(fasta, function(rec) rec$seq)
        desc <- sapply(fasta, function(rec) rec$desc)
        ans <- BStringViews(src, subjectClass, sep) # call the default method
        ans@desc <- desc
        ans
    }
)

# Called when 'src' is a BString (or derived) object.
# When not missing, 'subjectClass' must be "BString", "DNAString"
# or "RNAString".
setMethod("BStringViews", "BString",
    function(src, subjectClass, sep)
    {
        if (!missing(sep)) {
            # The semantic is: views are delimited by the occurences of 'sep'
            # in 'src' (a kind of strsplit() for BString objects).
            # Uncomment when normalize() and ! method are ready (see TODO file):
            #return(!normalize(matchPattern(sep, b, fixed=TRUE)))
            stop("'sep' not yet supported when 'src' is a \"BString\" object")
        }
        if (!missing(subjectClass) && subjectClass != class(src))
            src <- new(subjectClass, src)
        new("BStringViews", src, as.integer(1), src@length)
    }
)

# Called when 'src' is a BStringViews object.
# 'subjectClass' must be "BString", "DNAString" or "RNAString".
# The 'sep' arg is ignored.
setMethod("BStringViews", "BStringViews",
    function(src, subjectClass, sep)
    {
        if (!missing(sep))
            stop("'sep' not supported when 'src' is a \"BStringViews\" object")
        src@subject <- new(subjectClass, src@subject)
        src
    }
)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Standard generic methods

# 2 helper functions used by the show() method
catBsvFrameHeader <- function(iW, firstW, lastW, widthW)
{
    cat(format("", width=iW+1),
        format("first", width=firstW, justify="right"), " ",
        format("last", width=lastW, justify="right"), " ",
        format("width", width=widthW, justify="right"), "\n",
        sep="")
}
catBsvFrameLine <- function(x, i, iW, firstW, lastW, widthW)
{
    first <- x@first[i]
    last <- x@last[i]
    width <- last - first + 1
    snippetW <- 73 - iW - firstW - lastW - widthW
    cat(format(paste("[", i,"]", sep=""), width=iW, justify="right"), " ",
        format(first, width=firstW, justify="right"), " ",
        format(last, width=lastW, justify="right"), " ",
        format(width, width=widthW, justify="right"), " ",
        "|", bsViewSnippet(subject(x), first, last, snippetW), "|\n",
        sep="")
}

setMethod("show", "BStringViews",
    function(object)
    {
        subject <- subject(object)
        lsub <- length(subject)
        cat("  Views on a ", lsub, "-letter ",
            class(subject), " subject", sep="")
        #if (!is.null(subject@codec))
        #    cat(" with alphabet:", toString(subject@codec@letters))
        cat("\nSubject:", bsSnippet(subject, 70))
        cat("\nViews:")
        lo <- length(object)
        if (lo == 0)
            cat(" NONE\n")
        else {
            cat("\n")
            iW <- nchar(as.character(lo)) + 2 # 2 for the brackets
            firstMax <- max(object@first)
            firstW <- max(nchar(firstMax), nchar("first"))
            lastMax <- max(object@last)
            lastW <- max(nchar(lastMax), nchar("last"))
            widthMax <- max(width(object))
            widthW <- max(nchar(widthMax), nchar("width"))
            catBsvFrameHeader(iW, firstW, lastW, widthW)
            if (lo <= 19) {
                for (i in 1:lo)
                    catBsvFrameLine(object, i, iW, firstW, lastW, widthW)
            } else {
                for (i in 1:9)
                    catBsvFrameLine(object, i, iW, firstW, lastW, widthW)
                cat(format("...", width=iW, justify="right"),
                    " ",
                    format("...", width=firstW, justify="right"),
                    " ",
                    format("...", width=lastW, justify="right"),
                    " ",
                    format("...", width=widthW, justify="right"),
                    " ...\n", sep="")
                for (i in (lo-8):lo)
                    catBsvFrameLine(object, i, iW, firstW, lastW, widthW)
            }
        }
    }
)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Subsetting

setMethod("length", "BStringViews",
    function(x)
    {
        length(x@first)
    }
)

setMethod("[", "BStringViews",
    function(x, i, j, ..., drop)
    {
        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")
        if (missing(i))
            return(x)
        lx <- length(x)
        if (is.numeric(i)) {
            if (any(i < -lx) || any(i > lx))
                stop("subscript out of bounds")
        } else if (is.logical(i)) {
            if (length(i) > lx)
                stop("subscript out of bounds")
        } else {
            stop("invalid subscript type")
        }
        x@first <- x@first[i]
        x@last <- x@last[i]
        if (length(x@desc) != 0) {
            x@desc <- x@desc[i]
        }
        x
    }
)

setReplaceMethod("[", "BStringViews",
    function(x, i, j,..., value)
    {
        stop("attempt to modify the value of a \"BStringViews\" object")
    }
)

# Extract the i-th views of a BStringViews object as a BString object.
# Return a BString (or DNAString, or RNAString) object.
# Example:
#   bs <- BString("ABCD-1234-abcd")
#   bsv <- views(bs, 1:7, 13:7)
#   bsv[[3]]
#   bsv[[0]] # Return bs, same as subject(bsv)
#   views(bs)[[1]] # Returns bs too!
setMethod("[[", "BStringViews",
    function(x, i, j, ...)
    {
        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")
        if (missing(i))
            stop("invalid subscript type")
        if (length(i) < 1)
            stop("attempt to select less than one element")
        if (length(i) > 1)
            stop("attempt to select more than one element")
        if (i == 0)
            return(x@subject)
        if (i < 1 || i > length(x))
            stop("subscript out of bounds")
        first <- x@first[i]
        last <- x@last[i]
        if (first < 1 || last > length(x@subject))
            stop("view is out of limits")
        bsSubstr(x@subject, first, last)
    }
)

setReplaceMethod("[[", "BStringViews",
    function(x, i, j,..., value)
    {
        stop("attempt to modify the value of a \"BStringViews\" object")
    }
)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Equality

# Typical use:
#   v <- views(DNAString("TAATAATG"), -2:9, 0:11)
#   v == v[4]
#   v == v[1]
#   v2 <- views(DNAString("G"), 1, 3)
#   v == v2
# Also works if one side is a BString object:
#   v == DNAString("ATG")
#   RNAString("AUG") == v
# Whitespace matters:
#   v == "TG"
# But this doesn't work neither ("TG " can't be converted to a DNAString
# object):
#   v == "TG "

# Returns a logical vector of length max(length(x), length(y)).
# Recycle its arguments.
.equal <- function(x, y)
{
    if (class(y) != "BStringViews")
        y <- BStringViews(y, class(x@subject))
    lx <- length(x)
    ly <- length(y)
    if (lx < ly) {
        tmp <- x
        x <- y
        y <- tmp
        tmp <- lx
        lx <- ly
        ly <- tmp
    }
    if (ly == 0)
        return(logical(0))
    # Now we are sure that lx >= ly >= 1
    ans <- logical(lx)
    j <- 1
    for (i in 1:lx) {
        ans[i] <- bsIdenticalViews(x@subject, x@first[i], x@last[i],
                                   y@subject, y@first[j], y@last[j])
        # Recycle
        if (j < ly) j <- j + 1 else j <- 1
    }
    if (j != 1)
        warning(paste("longer object length",
                      "is not a multiple of shorter object length"))
    ans
}

# These methods are called if at least one side of the "==" (or "!=")
# operator is a "BStringViews" object. They have precedence over the
# corresponding methods defined for "BString" objects, i.e. they will
# be called if one side is a "BStringViews" object and the other side
# is a "BString" object.
setMethod("==", signature(e1="BStringViews", e2="BString"),
    function(e1, e2) .equal(e1, e2)
)
setMethod("==", signature(e1="BString", e2="BStringViews"),
    function(e1, e2) .equal(e2, e1)
)
setMethod("==", signature(e1="BStringViews"),
    function(e1, e2) .equal(e1, e2)
)
setMethod("==", signature(e2="BStringViews"),
    function(e1, e2) .equal(e2, e1)
)

setMethod("!=", signature(e1="BStringViews", e2="BString"),
    function(e1, e2) !.equal(e1, e2)
)
setMethod("!=", signature(e1="BString", e2="BStringViews"),
    function(e1, e2) !.equal(e2, e1)
)
setMethod("!=", signature(e1="BStringViews"),
    function(e1, e2) !.equal(e1, e2)
)
setMethod("!=", signature(e2="BStringViews"),
    function(e1, e2) !.equal(e2, e1)
)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("as.character", "BStringViews",
    function(x)
    {
        lx <- length(x)
        ans <- character(lx)
        if (lx >= 1) {
            for (i in 1:lx) {
                ans[i] <- bsView(x@subject, x@first[i], x@last[i])
            }
        }
        ans
    }
)

setMethod("toString", "BStringViews",
    function(x)
    {
        toString(as.character(x))
    }
)

setMethod("nchar", "BStringViews",
    function(x, type)
    {
        ls <- length(x@subject)
        ifelse(x@last<=ls, x@last, ls) - ifelse(x@first>=1, x@first, 1) + 1
    }
)

setMethod("as.matrix", "BStringViews",
    function(x)
    {
        cbind(x@first, x@last)
    }
)

setMethod("as.list", "BStringViews",
    function(x)
    {
        lx <- length(x)
        ans <- vector("list", lx)
        for (i in 1:lx) {
            ans[[i]] <- x[[i]]
        }
        ans
    }
)

