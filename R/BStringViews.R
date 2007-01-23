### =========================================================================
### The BStringViews class
### -------------------------------------------------------------------------
### A BStringViews object contains a set of views
### on the same BString object, the subject string.

setClassUnion("BStringLike", c("BString", "DNAString", "RNAString", "AAString"))

### See the initialization section below for the integrity checking
### of a BStringViews object.
setClass(
    "BStringViews",
    representation(
        subject="BStringLike",
        start="integer",
        end="integer",
        desc="character"   # store per-view comment (e.g. from FASTA file)
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods

setGeneric("subject", function(x) standardGeneric("subject"))
setMethod("subject", "BStringViews", function(x) x@subject)

### The 'substring' function uses 'first' and 'last'.
### The 'substr' function uses 'start' and 'stop'.
### But we prefer to use 'start' and 'end'. That's all!
setGeneric("start", function(x) standardGeneric("start"))
setMethod("start", "BStringViews", function(x) x@start)

setGeneric("first", function(x) standardGeneric("first"))
setMethod("first", "BStringViews", function(x) {.Deprecated("start"); start(x)})

setGeneric("end", function(x) standardGeneric("end"))
setMethod("end", "BStringViews", function(x) x@end)

setGeneric("last", function(x) standardGeneric("last"))
setMethod("last", "BStringViews", function(x) {.Deprecated("end"); end(x)})

### We choose to call this method 'width' and not 'length' because
### we want to define the length of a BStringViews object as the number
### of views contained in it.
### Another option was to call it 'nchar' but the width of a view is not
### necesarily equal to the number of letters that it contains (this happens
### when the view is out of limits).
setGeneric("width", function(x) standardGeneric("width"))
setMethod("width", "BStringViews", function(x) x@end - x@start + 1)

setGeneric("desc", function(x) standardGeneric("desc"))
setMethod("desc", "BStringViews", function(x) x@desc)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The 'show' method

### The 2 helper functions below convert a given view on a BString object
### into a character-string.
### Both assume that 'start' <= 'end' (so they don't check it) and
### padd the result with spaces to produce the "margin effect"
### if 'start' or 'end' are out of limits.

### nchar(BStringViews.get_view(x, start, end)) is always end-start+1
BStringViews.get_view <- function(x, start, end)
{
    lx <- length(x)
    if (end < 1 || start > lx)
            return(format("", width=end-start+1))
    Lmargin <- ""
    if (start < 1) {
        Lmargin <- format("", width=1-start)
        start <- 1
    }
    Rmargin <- ""
    if (end > lx) {
        Rmargin <- format("", width=end-lx)
        end <- lx
    }
    paste(Lmargin, BString.read(x, start, end), Rmargin, sep="")
}

### nchar(BStringViews.get_snippet(x, start, end, snippetWidth)) is <= snippetWidth
BStringViews.get_snippet <- function(x, start, end, snippetWidth)
{
    if (snippetWidth < 7)
        snippetWidth <- 7
    width <- end - start + 1
    if (width <= snippetWidth) {
        BStringViews.get_view(x, start, end)
    } else {
        w1 <- (snippetWidth - 2) %/% 2
        w2 <- (snippetWidth - 3) %/% 2
        paste(BStringViews.get_view(x, start, start+w1-1),
              "...",
              BStringViews.get_view(x, end-w2+1, end), sep="")
    }
}

BStringViews.show_frame_header <- function(iW, startW, endW, widthW)
{
    cat(format("", width=iW+1),
        format("start", width=startW, justify="right"), " ",
        format("end", width=endW, justify="right"), " ",
        format("width", width=widthW, justify="right"), "\n",
        sep="")
}

BStringViews.show_frame_line <- function(x, i, iW, startW, endW, widthW)
{
    start <- x@start[i]
    end <- x@end[i]
    width <- end - start + 1
    snippetWidth <- 73 - iW - startW - endW - widthW
    cat(format(paste("[", i,"]", sep=""), width=iW, justify="right"), " ",
        format(start, width=startW, justify="right"), " ",
        format(end, width=endW, justify="right"), " ",
        format(width, width=widthW, justify="right"), " ",
        "[", BStringViews.get_snippet(subject(x), start, end, snippetWidth), "]\n",
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
        cat("\nSubject:", BString.get_snippet(subject, 70))
        cat("\nViews:")
        lo <- length(object)
        if (lo == 0)
            cat(" NONE\n")
        else {
            cat("\n")
            iW <- nchar(as.character(lo)) + 2 # 2 for the brackets
            startMax <- max(object@start)
            startW <- max(nchar(startMax), nchar("start"))
            endMax <- max(object@end)
            endW <- max(nchar(endMax), nchar("end"))
            widthMax <- max(width(object))
            widthW <- max(nchar(widthMax), nchar("width"))
            BStringViews.show_frame_header(iW, startW, endW, widthW)
            if (lo <= 19) {
                for (i in 1:lo)
                    BStringViews.show_frame_line(object, i, iW, startW, endW, widthW)
            } else {
                for (i in 1:9)
                    BStringViews.show_frame_line(object, i, iW, startW, endW, widthW)
                cat(format("...", width=iW, justify="right"),
                    " ",
                    format("...", width=startW, justify="right"),
                    " ",
                    format("...", width=endW, justify="right"),
                    " ",
                    format("...", width=widthW, justify="right"),
                    " ...\n", sep="")
                for (i in (lo-8):lo)
                    BStringViews.show_frame_line(object, i, iW, startW, endW, widthW)
            }
        }
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting

setMethod("length", "BStringViews",
    function(x)
    {
        length(x@start)
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
        x@start <- x@start[i]
        x@end <- x@end[i]
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

### Extract the i-th views of a BStringViews object as a BString object.
### Return a "BString" (or one of its derivations) object.
### Example:
###   bs <- BString("ABCD-1234-abcd")
###   bsv <- views(bs, 1:7, 13:7)
###   bsv[[3]]
###   bsv[[0]] # Return bs, same as subject(bsv)
###   views(bs)[[1]] # Returns bs too!
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
        start <- x@start[i]
        end <- x@end[i]
        if (start < 1 || end > length(x@subject))
            stop("view is out of limits")
        BString.substring(x@subject, start, end)
    }
)

setReplaceMethod("[[", "BStringViews",
    function(x, i, j,..., value)
    {
        stop("attempt to modify the value of a \"BStringViews\" object")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Equality

### Assume that 'start1', 'end1', 'start2', 'end2' are single integers
### and that start1 <= end1 and start2 <= end2.
BStringViews.view1_equal_view2 <- function(x1, start1, end1, x2, start2, end2)
{
    one <- as.integer(1)
    w1 <- end1 - start1 + one
    w2 <- end2 - start2 + one
    if (w1 != w2)
        return(FALSE)

    lx1 <- length(x1)
    isBlank1 <- end1 < one || start1 > lx1
    lx2 <- length(x2)
    isBlank2 <- end2 < one || start2 > lx2
    if (isBlank1 && isBlank2)
        return(TRUE)
    if (isBlank1 || isBlank2)
        return(FALSE)

    # Left margin
    LmarginSize1 <- start1 < one
    LmarginSize2 <- start2 < one
    if (LmarginSize1 != LmarginSize2)
        return(FALSE)
    if (LmarginSize1) {
        # Both views have a left margin
        if (start1 != start2)
            return(FALSE)
        start1 <- one
        start2 <- one
    }

    # Right margin
    RmarginSize1 <- end1 > lx1
    RmarginSize2 <- end2 > lx2
    if (RmarginSize1 != RmarginSize2)
        return(FALSE)
    if (RmarginSize1) {
        # Both views have a right margin
        if (end1 - lx1 != end2 - lx2)
            return(FALSE)
        end1 <- lx1
        end2 <- lx2
    }

    # At this point, we can trust that 1 <= start1 <= end1 <= lx1
    # and that 1 <= start2 <= end2 <= lx2 so we can call unsafe
    # function BString.substring() with no fear...
    BString.substring(x1, start1, end1) == BString.substring(x2, start2, end2)
}

### Returns a logical vector of length max(length(x), length(y)).
### Recycle its arguments.
BStringViews.equal <- function(x, y)
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
        ans[i] <- BStringViews.view1_equal_view2(x@subject, x@start[i], x@end[i],
                                                 y@subject, y@start[j], y@end[j])
        # Recycle
        if (j < ly) j <- j + 1 else j <- 1
    }
    if (j != 1)
        warning(paste("longer object length",
                      "is not a multiple of shorter object length"))
    ans
}

### These methods are called if at least one side of the "==" (or "!=")
### operator is a "BStringViews" object. They have precedence over the
### corresponding methods defined for "BString" objects, i.e. they will
### be called if one side is a "BStringViews" object and the other side
### is a "BString" object.
### Typical use:
###   v <- views(DNAString("TAATAATG"), -2:9, 0:11)
###   v == v[4]
###   v == v[1]
###   v2 <- views(DNAString("G"), 1, 3)
###   v == v2
### Also works if one side is a BString object:
###   v == DNAString("ATG")
###   RNAString("AUG") == v
### Whitespace matters:
###   v == "TG"
### But this doesn't work neither ("TG " can't be converted to a DNAString
### object):
###   v == "TG "

setMethod("==", signature(e1="BStringViews", e2="BString"),
    function(e1, e2) BStringViews.equal(e1, e2)
)
setMethod("==", signature(e1="BString", e2="BStringViews"),
    function(e1, e2) BStringViews.equal(e2, e1)
)
setMethod("==", signature(e1="BStringViews"),
    function(e1, e2) BStringViews.equal(e1, e2)
)
setMethod("==", signature(e2="BStringViews"),
    function(e1, e2) BStringViews.equal(e2, e1)
)

setMethod("!=", signature(e1="BStringViews", e2="BString"),
    function(e1, e2) !BStringViews.equal(e1, e2)
)
setMethod("!=", signature(e1="BString", e2="BStringViews"),
    function(e1, e2) !BStringViews.equal(e2, e1)
)
setMethod("!=", signature(e1="BStringViews"),
    function(e1, e2) !BStringViews.equal(e1, e2)
)
setMethod("!=", signature(e2="BStringViews"),
    function(e1, e2) !BStringViews.equal(e2, e1)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Other standard generic methods

setMethod("as.character", "BStringViews",
    function(x)
    {
        lx <- length(x)
        ans <- character(lx)
        if (lx >= 1) {
            for (i in 1:lx) {
                ans[i] <- BStringViews.get_view(x@subject, x@start[i], x@end[i])
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
        ifelse(x@end<=ls, x@end, ls) - ifelse(x@start>=1, x@start, 1) + 1
    }
)

setMethod("as.matrix", "BStringViews",
    function(x, ...)
    {
        cbind(x@start, x@end)
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

