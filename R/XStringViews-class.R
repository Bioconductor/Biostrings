### =========================================================================
### XStringViews objects
### -------------------------------------------------------------------------
###
### The XStringViews class is the basic container for storing a set of views
### (start/end locations) on the same XString object, called the "subject"
### string.
###

setClass("XStringViews",
    contains="Views",
    representation(
        subject="XString"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Unsafe constructor (not exported). Use only when 'start' and 'width' are
### guaranteed to be valid.
###

unsafe.newXStringViews <- function(subject, start, width)
    new2("XStringViews", subject=subject,
                         ranges=IRanges(start=start, width=width),
                         check=FALSE)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### User-friendly constructor.
###

setMethod("Views", "XString",
    function(subject, start=NULL, end=NULL, width=NULL, names=NULL)
        IRanges:::new_Views(subject, start=start,
                                     end=end,
                                     width=width,
                                     names=names,
                                     Class="XStringViews")
)

setMethod("Views", "character",
    function(subject, start=NULL, end=NULL, width=NULL, names=NULL)
    {
        xsubject <- XString(NULL, subject)
        Views(xsubject, start=start, end=end, width=width, names=names)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor-like methods.
###

setMethod("nchar", "XStringViews",
    function(x, type="chars", allowNA=FALSE)
    {
        if (length(x) == 0)
            return(integer(0))
        start0 <- pmax.int(start(x), 1L)
        end0 <- pmin.int(end(x), nchar(subject(x)))
        ans <- end0 - start0 + 1L
        ans[ans < 0L] <- 0L
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "seqtype" and "seqtype<-" methods.
###

setMethod("seqtype", "XStringViews", function(x) seqtype(subject(x)))

### Does NOT downgrade 'x' to an XStringViews instance! (endomorphism)
setReplaceMethod("seqtype", "XStringViews",
    function(x, value)
    {
        ## could be done with 'seqtype(subject(x)) <- value'
        ## if `subject<-` was available
        subject <- subject(x)
        seqtype(subject) <- value
        x@subject <- subject
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

### From XStringViews to XStringSet.
fromXStringViewsToStringSet <- function(x,
                                 out.of.limits=c("ok", "warning", "error"),
                                 use.names=FALSE)
{
    out.of.limits <- match.arg(out.of.limits)
    ans_ranges <- restrict(as(x, "IRanges"), start=1L, end=nchar(subject(x)),
                           keep.all.ranges=TRUE,
                           use.names=use.names)
    if (out.of.limits != "ok" && any(width(ans_ranges) < width(x))) {
        if (out.of.limits == "warning")
            warning("trimming \"out of limits\" views")
        else
            stop("'x' has \"out of limits\" views")
    }
    extractList(subject(x), ans_ranges)
}

### We need this so that B/DNA/RNA/AAStringSet() used below work on an
### XStringViews object.
setMethod("XStringSet", "XStringViews",
    function(seqtype, x, start=NA, end=NA, width=NA, use.names=TRUE)
    {
        y <- fromXStringViewsToStringSet(x, out.of.limits="warning",
                                         use.names=use.names)
        XStringSet(seqtype, y,
                   start=start, end=end, width=width, use.names=TRUE)
    }
)

setAs("XStringViews", "XStringSet",
    function(from) fromXStringViewsToStringSet(from, out.of.limits="warning",
                                               use.names=TRUE))

setAs("XStringViews", "BStringSet", function(from) BStringSet(from))
setAs("XStringViews", "DNAStringSet", function(from) DNAStringSet(from))
setAs("XStringViews", "RNAStringSet", function(from) RNAStringSet(from))
setAs("XStringViews", "AAStringSet", function(from) AAStringSet(from))

### From XStringSet to XStringViews.

.XStringSetAsViews <- function(from) successiveViews(unlist(from), width(from))

setAs("XStringSet", "Views", .XStringSetAsViews)
setAs("XStringSet", "XStringViews", .XStringSetAsViews)

setMethod("as.data.frame", "XStringViews",
          function (x, row.names = NULL, optional = FALSE, ...)
          {
              as.data.frame(as(x, "XStringSet"), row.names, optional, ...)
          })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method.
###

### The 2 helper functions below convert a given view on an XString object
### into a character-string.
### Both assume that 'start' <= 'end' (so they don't check it) and
### padd the result with spaces to produce the "margin effect"
### if 'start' or 'end' are out of limits.

### nchar(XStringViews.get_view(x, start, end)) is always end-start+1
XStringViews.get_view <- function(x, start, end)
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
    paste(Lmargin, XString.read(x, start, end), Rmargin, sep="")
}

### nchar(XStringViews.get_snippet(x, start, end, snippetWidth)) is <= snippetWidth
XStringViews.get_snippet <- function(x, start, end, snippetWidth)
{
    if (snippetWidth < 7)
        snippetWidth <- 7
    width <- end - start + 1
    if (width <= snippetWidth) {
        XStringViews.get_view(x, start, end)
    } else {
        w1 <- (snippetWidth - 2) %/% 2
        w2 <- (snippetWidth - 3) %/% 2
        paste(XStringViews.get_view(x, start, start+w1-1),
              "...",
              XStringViews.get_view(x, end-w2+1, end), sep="")
    }
}

XStringViews.show_vframe_header <- function(iW, startW, endW, widthW)
{
    cat(format("", width=iW+1),
        format("start", width=startW, justify="right"), " ",
        format("end", width=endW, justify="right"), " ",
        format("width", width=widthW, justify="right"), "\n",
        sep="")
}

XStringViews.show_vframe_line <- function(x, i, iW, startW, endW, widthW)
{
    start <- start(x)[i]
    end <- end(x)[i]
    width <- end - start + 1
    snippetWidth <- getOption("width") - 6 - iW - startW - endW - widthW
    cat(format(paste("[", i,"]", sep=""), width=iW, justify="right"), " ",
        format(start, width=startW, justify="right"), " ",
        format(end, width=endW, justify="right"), " ",
        format(width, width=widthW, justify="right"), " ",
        "[", XStringViews.get_snippet(subject(x), start, end, snippetWidth), "]\n",
        sep="")
}

XStringViews.show_vframe <- function(x)
{
    nhead <- get_showHeadLines()
    ntail <- get_showTailLines()
    cat("\nviews:")
    lx <- length(x)
    if (lx == 0)
        cat(" NONE\n")
    else {
        cat("\n")
        iW <- nchar(as.character(lx)) + 2 # 2 for the brackets
        startMax <- max(start(x))
        startW <- max(nchar(startMax), nchar("start"))
        endMax <- max(end(x))
        endW <- max(nchar(endMax), nchar("end"))
        widthMax <- max(width(x))
        widthW <- max(nchar(widthMax), nchar("width"))
        XStringViews.show_vframe_header(iW, startW, endW, widthW)
        if (lx <= nhead + ntail +1) {
            for (i in seq_len(lx))
                XStringViews.show_vframe_line(x, i, iW, startW, endW, widthW)
        } else {
            if (nhead > 0)
                for (i in seq_len(nhead))
                    XStringViews.show_vframe_line(x, i, iW, startW, endW, 
                                                  widthW)
            cat(format("...", width=iW, justify="right"),
                " ",
                format("...", width=startW, justify="right"),
                " ",
                format("...", width=endW, justify="right"),
                " ",
                format("...", width=widthW, justify="right"),
                " ...\n", sep="")
            if (ntail > 0)
                for (i in (lx-ntail+1L):lx)
                    XStringViews.show_vframe_line(x, i, iW, startW, endW, 
                                                  widthW)
        }
    }
}

setMethod("show", "XStringViews",
    function(object)
    {
        subject <- subject(object)
        lsub <- length(subject)
        cat("  Views on a ", lsub, "-letter ", class(subject), " subject", sep="")
        cat("\nsubject:", toSeqSnippet(subject, getOption("width") - 9))
        XStringViews.show_vframe(object)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Equality.
###

### Assume that 'start1', 'end1', 'start2', 'end2' are single integers
### and that start1 <= end1 and start2 <= end2.
XStringViews.view1_equal_view2 <- function(x1, start1, end1, x2, start2, end2)
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
    # and that 1 <= start2 <= end2 <= lx2.
    subseq(x1, start=start1, end=end1) == subseq(x2, start=start2, end=end2)
}

### 'x' and 'y' must be XStringViews objects.
### Returns a logical vector of length max(length(x), length(y)).
### Recycle its arguments.
XStringViews.equal <- function(x, y)
{
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
    for (i in seq_len(lx)) {
        ans[i] <- XStringViews.view1_equal_view2(
                      subject(x), start(x)[i], end(x)[i],
                      subject(y), start(y)[j], end(y)[j])
        # Recycle
        if (j < ly) j <- j + 1 else j <- 1
    }
    if (j != 1)
        warning(paste("longer object length",
                      "is not a multiple of shorter object length"))
    ans
}

### These methods are called if at least one side of the "==" (or "!=")
### operator is an XStringViews object. They have precedence over the
### corresponding methods defined for XString objects, i.e. they will
### be called if one side is an XStringViews object and the other side
### is an XString object.
### Typical use:
###   v <- Views(DNAString("TAATAATG"), start=-2:9, end=0:11)
###   v == v[4]
###   v == v[1]
###   v2 <- Views(DNAString("G"), start=1, end=3)
###   v == v2
### Also works if one side is an XString object:
###   v == DNAString("ATG")
###   RNAString("AUG") == v
### Whitespace matters:
###   v == "TG"
### But this doesn't work neither ("TG " can't be converted to a DNAString
### object):
###   v == "TG "

setMethod("==", signature(e1="XStringViews", e2="XStringViews"),
    function(e1, e2)
    {
        if (!comparable_seqtypes(seqtype(e1), seqtype(e2))) {
            class1 <- class(subject(e1))
            class2 <- class(subject(e2))
            stop("comparison between XStringViews objects with subjects of ",
                 "class \"", class1, "\" and \"", class2, "\" ",
                 "is not supported")
        }
        XStringViews.equal(e1, e2)
    }
)
setMethod("==", signature(e1="XStringViews", e2="XString"),
    function(e1, e2)
    {
        if (!comparable_seqtypes(seqtype(e1), seqtype(e2))) {
            class1 <- class(subject(e1))
            class2 <- class(e2)
            stop("comparison between an XStringViews object with a subject of ",
                 "class \"", class1, "\" and a \"", class2, "\" instance ",
                 "is not supported")
        }
        XStringViews.equal(e1, as(e2, "Views"))
    }
)
setMethod("==", signature(e1="XStringViews", e2="character"),
    function(e1, e2)
    {
        if (!is(subject(e1), "BString"))
            stop("comparison between an XStringViews object with a subject of ",
                 "class \"", class(subject(e1)), "\" and a character vector ",
                 "is not supported")
        if (length(e2) == 0 || any(e2 %in% c("", NA)))
            stop("comparison between an XStringViews object and a character ",
                 "vector of length 0 or with empty strings or NAs ",
                 "is not supported")
        XStringViews.equal(e1, as(BStringSet(e2), "Views"))
    }
)
setMethod("==", signature(e1="XString", e2="XStringViews"),
    function(e1, e2) e2 == e1
)
setMethod("==", signature(e1="character", e2="XStringViews"),
    function(e1, e2) e2 == e1
)

setMethod("!=", signature(e1="XStringViews", e2="XStringViews"),
    function(e1, e2) !(e1 == e2)
)
setMethod("!=", signature(e1="XStringViews", e2="XString"),
    function(e1, e2) !(e1 == e2)
)
setMethod("!=", signature(e1="XStringViews", e2="character"),
    function(e1, e2) !(e1 == e2)
)
setMethod("!=", signature(e1="XString", e2="XStringViews"),
    function(e1, e2) !(e1 == e2)
)
setMethod("!=", signature(e1="character", e2="XStringViews"),
    function(e1, e2) !(e1 == e2)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Other methods.
###

setMethod("as.character", "XStringViews",
    function(x, use.names=TRUE, check.limits=TRUE)
    {
        if (!isTRUEorFALSE(check.limits))
            stop("'check.limits' must be TRUE or FALSE")
        if (check.limits)
            out.of.limits <- "error"
        else
            out.of.limits <- "warning"
        y <- fromXStringViewsToStringSet(x, out.of.limits=out.of.limits,
                                         use.names=use.names)
        as.character(y)
    }
)

setMethod("as.matrix", "XStringViews",
    function(x, use.names=TRUE)
    {
        ## TODO: Supress this warning in BioC 2.12.
        msg <- c("as.matrix() on an XStringViews object 'x' has changed ",
                 "behavior: now the\n  views in 'x' must be of equal width ",
                 "and each view is converted into a row of\n",
                 "  single characters. To achieve the old behavior, ",
                 "do 'as.matrix(ranges(x))'.\n  To supress this warning, do ",
                 "'suppressWarnings(as.matrix(x))'.\n  This warning will be ",
                 "removed in BioC 2.12.")
        warning(msg)
        y <- fromXStringViewsToStringSet(x, out.of.limits="error",
                                         use.names=use.names)
        as.matrix(y)
    }
)

setMethod("toString", "XStringViews",
    function(x, ...)
    {
        toString(as.character(x), ...)
    }
)

