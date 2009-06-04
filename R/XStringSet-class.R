### =========================================================================
### XStringSet objects
### -------------------------------------------------------------------------
###
### The XStringSet class is a container for storing a set of XString objects
### of the same base type (e.g. all elements are BString objects or they are
### all DNAString objects etc).
###
### The current implementation only allows for storage of a set of strings
### that belong to the same RawPtr object i.e. all the elements of an XStringSet
### object must be substrings of a common string called the "super string".
### The old XStringList container (Biostrings 2.8) didn't have this limitation:
### it could hold XString objects that pointed to different RawPtr objects but
### it was so slow that I decided to replace it by the much more efficient
### XStringSet container.
### Maybe the best of both world, or at least a good enough trade-off, could be
### obtained by defining the XStringSet class like this:
###
###   setClass("RawPtrViews",
###     contains="IRanges",
###     representation(
###         subject="RawPtr"
###     )
###   )
###
###   setClass("XStringSet",
###     contains="Sequence",
###     representation(
###         "VIRTUAL",
###         xrvlist="list",   # a list of RawPtrViews objects
###         strong="integer",
###         weak="integer"
###     )
###   )
###
### x@strong and x@weak (both of length the number of elements) are maps used
### to "find" the i-th element in x i.e. x[[i]] can be efficiently extracted
### with:
###     xrv <- x@xrvlist[[x@strong[i]]]
###     start <- start(xrv)[x@weak[i]]
###     width <- width(xrv)[x@weak[i]]
###     new(class, xdata=xrv@subject, offset=start-1L, length=width)
### This new XStringSet container would still be almost as efficient as the old
### one (at least for the restricted set of use cases that the old one was
### supporting i.e. when x@xrvlist holds a list of length 1) and the
### flexibility of the old XStringList container.
###
### Some notable differences between XStringSet and XStringViews objects:
###   - the "show" methods produce different output
###   - the views in an XStringViews object can't have a 0-width
###   - the views in an XStringViews object can be out of limits
###

setClass("XStringSet",
    contains="Sequence",
    representation(
        "VIRTUAL",
        super="XString",
        ranges="IRanges"
    )
)

setClass("BStringSet",
    contains="XStringSet",
    representation(
        super="BString"
    )
)
setClass("DNAStringSet",
    contains="XStringSet",
    representation(
        super="DNAString"
    )
)
setClass("RNAStringSet",
    contains="XStringSet",
    representation(
        super="RNAString"
    )
)
setClass("AAStringSet",
    contains="XStringSet",
    representation(
        super="AAString"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor-like methods.
###

### super() must remain for internal use only. Do not export!
setGeneric("super", function(x) standardGeneric("super"))
setMethod("super", "XStringSet", function(x) x@super)

setMethod("length", "XStringSet", function(x) length(x@ranges))

setMethod("width", "XStringSet", function(x) width(x@ranges))

setMethod("width", "character",
    function(x)
    {
        if (any(is.na(x)))
            stop("NAs in 'x' are not supported")
        nchar(x, type="bytes")
    }
)

setMethod("nchar", "XStringSet",
    function(x, type="chars", allowNA=FALSE) width(x)
)

setMethod("names", "XStringSet", function(x) names(x@ranges))

setReplaceMethod("names", "XStringSet",
    function(x, value)
    {
        names(x@ranges) <- value
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "subseq" endomorphism and related transformations.
###

setMethod("narrow", "character",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE)
    {
        if (!normargUseNames(use.names))
            names(x) <- NULL
        x_width <- width(x)
        solved_SEW <- solveUserSEW(x_width, start=start, end=end, width=width)
        substr(x, start=start(solved_SEW), stop=end(solved_SEW))
    }
)

setMethod("narrow", "XStringSet",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE)
    {
        x@ranges <- narrow(x@ranges, start=start, end=end, width=width,
                           use.names=use.names)
        x
    }
)

setMethod("subseq", "character",
    function(x, start=NA, end=NA, width=NA)
        narrow(x, start=start, end=end, width=width)
)

setMethod("subseq", "XStringSet",
    function(x, start=NA, end=NA, width=NA)
        narrow(x, start=start, end=end, width=width)
)

setMethod("threebands", "character",
    function(x, start=NA, end=NA, width=NA)
    {
        names(x) <- NULL
        x_width <- width(x)
        solved_SEW <- solveUserSEW(x_width, start=start, end=end, width=width)
        left <- substr(x, start=1L, stop=start(solved_SEW)-1L)
        middle <- substr(x, start=start(solved_SEW), stop=end(solved_SEW))
        right <- substr(x, start=end(solved_SEW)+1L, stop=x_width)
        list(left=left, middle=middle, right=right)
    }
)

setMethod("threebands", "XStringSet",
    function(x, start=NA, end=NA, width=NA)
    {
        threeranges <- threebands(x@ranges, start=start, end=end, width=width)
        left <- right <- x
        left@ranges <- threeranges$left
        x@ranges <- threeranges$middle
        right@ranges <- threeranges$right
        list(left=left, middle=x, right=right)
    }
)

setReplaceMethod("subseq", "character",
    function(x, start=NA, end=NA, width=NA, value)
    {
        bands <- threebands(x, start=start, end=end, width=width)
        paste(bands$left, value, bands$right, sep="")
    }
)

setReplaceMethod("subseq", "XStringSet",
    function(x, start=NA, end=NA, width=NA, value)
    {
        bands <- threebands(x, start=start, end=end, width=width)
        if (is.null(value))
            xscat(bands$left, bands$right)
        else
            xscat(bands$left, value, bands$right)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Unsafe constructor (not exported). Use only when 'ranges' is guaranteed
### to contain valid ranges on 'super' i.e. ranges that are within the limits
### of 'super'.
###

unsafe.newXStringSet <- function(super, ranges, use.names=FALSE, names=NULL)
{
    class <- paste(xsbasetype(super), "StringSet", sep="")
    ans <- new2(class, super=super, ranges=ranges, check=FALSE)
    if (normargUseNames(use.names))
        names(ans) <- names
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "xsbasetype" and "xsbasetype<-" methods.
###

setMethod("xsbasetype", "XStringSet", function(x) xsbasetype(super(x)))

### This is an endomorphism iff 'basetype' is NULL, otherwise it is NOT!
.XStringSet.compact <- function(x, basetype=NULL)
{
    from_basetype <- xsbasetype(x)
    from_baseclass <- paste(from_basetype, "String", sep="")
    if (is.null(basetype)) {
        basetype <- from_basetype
        baseclass <- from_baseclass
    } else {
        baseclass <- paste(basetype, "String", sep="")
    }
    lkup <- get_xsbasetypes_conversion_lookup(from_basetype, basetype)
    ## The frame is the strict minimal set of regions in 'super(x)' that need
    ## to be copied. Hence compacting 'x' returns an XStringSet object 'y'
    ## where 'length(super(y))' can be significantly smaller than
    ## 'length(super(x))' especially if the elements in 'x' cover a small part
    ## of 'super(x)'.
    frame <- reduce(x@ranges, with.inframe.attrib=TRUE)
    xdata <- .Call("new_RawPtr_from_XString",
                   super(x), start(frame), width(frame), lkup,
                   PACKAGE="Biostrings")
    ans_super <- new(baseclass, xdata=xdata, length=length(xdata))
    ans_ranges <- attr(frame, "inframe")
    if (is.null(basetype)) {
        ## Endomorphism
        x@super <- ans_super
        x@ranges <- IRanges:::unsafe.update(x@ranges,
                                            start=start(ans_ranges),
                                            width=width(ans_ranges))
        return(x)
    }
    ## NOT an endomorphism! (downgrades 'x' to a B/DNA/RNA/AAStringSet instance)
    unsafe.newXStringSet(ans_super, ans_ranges, use.names=TRUE, names=names(x))
}

### Downgrades 'x' to a B/DNA/RNA/AAStringSet instance!
setReplaceMethod("xsbasetype", "XStringSet",
    function(x, value)
    {
        from_basetype <- xsbasetype(x)
        lkup <- get_xsbasetypes_conversion_lookup(from_basetype, value)
        if (!is.null(lkup))
            return(.XStringSet.compact(x, basetype=value))
        ans_class <- paste(value, "StringSet", sep="")
        if (is(x, ans_class)) {
            ans_super <- super(x)
        } else {
            ans_super <- XString(value, super(x))
        }
        unsafe.newXStringSet(ans_super, x@ranges, use.names=TRUE, names=names(x))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### compact()
###

setGeneric("compact", function(x, ...) standardGeneric("compact"))

setMethod("compact", "XStringSet", .XStringSet.compact)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helper functions (NOT exported) used by the versatile (and user friendly)
### constructors and coercion methods below.
###

.charToXString <- function(x, solved_SEW, basetype)
{
    class <- paste(basetype, "String", sep="")
    proto <- newEmptyXString(class)
    xdata <- .Call("new_RawPtr_from_STRSXP",
                   x, start(solved_SEW), width(solved_SEW), "", xs_enc_lkup(proto),
                   PACKAGE="Biostrings")
    new(class, xdata=xdata, length=length(xdata))
}

.charToXStringSet <- function(x, start, end, width, use.names, basetype)
{
    solved_SEW <- solveUserSEW(width(x), start=start, end=end, width=width)
    ans_super <- .charToXString(x, solved_SEW, basetype)
    ans_ranges <- successiveIRanges(width(solved_SEW))
    unsafe.newXStringSet(ans_super, ans_ranges,
                         use.names=use.names, names=names(x))
}

.XStringToXStringSet <- function(x, start, end, width, use.names, basetype)
{
    ans_super <- XString(basetype, x, start=start, end=end, width=width)
    ans_ranges <- new2("IRanges", start=1L, width=length(ans_super), check=FALSE)
    unsafe.newXStringSet(ans_super, ans_ranges)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The user-friendly versatile constructors.
###
### All these constructors use the uSEW (user-specified Start/End/Width)
### interface.
###

setGeneric("XStringSet", signature="x",
    function(basetype, x, start=NA, end=NA, width=NA, use.names=TRUE)
        standardGeneric("XStringSet")
)
setMethod("XStringSet", "character",
    function(basetype, x, start=NA, end=NA, width=NA, use.names=TRUE)
    {
        if (is.null(basetype))
            basetype <- "B"
        .charToXStringSet(x, start, end, width, use.names, basetype)
    }
)
### Just because of the silly "AsIs" objects found in the probe packages
### (e.g. drosophila2probe$sequence)
setMethod("XStringSet", "AsIs",
    function(basetype, x, start=NA, end=NA, width=NA, use.names=TRUE)
    {
        if (!is.character(x))
            stop("unsuported input type")
        class(x) <- "character" # keeps the names (unlike as.character())
	.charToXStringSet(x, start, end, width, use.names, basetype)
    }
)
setMethod("XStringSet", "XString",
    function(basetype, x, start=NA, end=NA, width=NA, use.names=TRUE)
        .XStringToXStringSet(x, start, end, width, use.names, basetype)
)
setMethod("XStringSet", "XStringSet",
    function(basetype, x, start=NA, end=NA, width=NA, use.names=TRUE)
    {
        ans <- narrow(x, start=start, end=end, width=width, use.names=use.names)
        ## `xsbasetype<-` must be called even when user supplied 'basetype' is
        ## NULL because we want to enforce downgrade to a B/DNA/RNA/AAStringSet
        ## instance
        if (is.null(basetype))
            basetype <- xsbasetype(x)
        xsbasetype(ans) <- basetype
        ans
    }
)

BStringSet <- function(x=character(), start=NA, end=NA, width=NA, use.names=TRUE)
    XStringSet("B", x, start=start, end=end, width=width,
                             use.names=use.names)
DNAStringSet <- function(x=character(), start=NA, end=NA, width=NA, use.names=TRUE)
    XStringSet("DNA", x, start=start, end=end, width=width,
                               use.names=use.names)
RNAStringSet <- function(x=character(), start=NA, end=NA, width=NA, use.names=TRUE)
    XStringSet("RNA", x, start=start, end=end, width=width,
                               use.names=use.names)
AAStringSet <- function(x=character(), start=NA, end=NA, width=NA, use.names=TRUE)
    XStringSet("AA", x, start=start, end=end, width=width,
                              use.names=use.names)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

setAs("XStringSet", "BStringSet",
    function(from) {xsbasetype(from) <- "B"; from}
)
setAs("XStringSet", "DNAStringSet",
    function(from) {xsbasetype(from) <- "DNA"; from}
)
setAs("XStringSet", "RNAStringSet",
    function(from) {xsbasetype(from) <- "RNA"; from}
)
setAs("XStringSet", "AAStringSet",
    function(from) {xsbasetype(from) <- "AA"; from}
)

setAs("character", "BStringSet", function(from) BStringSet(from))
setAs("character", "DNAStringSet", function(from) DNAStringSet(from))
setAs("character", "RNAStringSet", function(from) RNAStringSet(from))
setAs("character", "AAStringSet", function(from) AAStringSet(from))
setAs("character", "XStringSet", function(from) BStringSet(from))

setAs("XString", "BStringSet", function(from) BStringSet(from))
setAs("XString", "DNAStringSet", function(from) DNAStringSet(from))
setAs("XString", "RNAStringSet", function(from) RNAStringSet(from))
setAs("XString", "AAStringSet", function(from) AAStringSet(from))
setAs("XString", "XStringSet", function(from) XStringSet(xsbasetype(from), from))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method.
###

.namesW <- 20

.XStringSet.show_frame_header <- function(iW, widthW, with.names)
{
    cat(format("", width=iW+1),
        format("width", width=widthW, justify="right"),
        sep="")
    if (with.names) {
        cat(format(" seq", width=getOption("width")-iW-widthW-.namesW-1),
            format("names", width=.namesW, justify="left"),
            sep="")
    } else {
        cat(" seq")
    }
    cat("\n")
}

.XStringSet.show_frame_line <- function(x, i, iW, widthW)
{
    width <- nchar(x)[i]
    snippetWidth <- getOption("width") - 2 - iW - widthW
    if (!is.null(names(x)))
        snippetWidth <- snippetWidth - .namesW - 1
    seq_snippet <- toSeqSnippet(x[[i]], snippetWidth)
    if (!is.null(names(x)))
        seq_snippet <- format(seq_snippet, width=snippetWidth)
    cat(format(paste("[", i,"]", sep=""), width=iW, justify="right"), " ",
        format(width, width=widthW, justify="right"), " ",
        seq_snippet,
        sep="")
    if (!is.null(names(x))) {
        snippet_name <- names(x)[i]
        if (is.na(snippet_name))
            snippet_name <- "<NA>"
        else if (nchar(snippet_name) > .namesW)
            snippet_name <- paste(substr(snippet_name, 1, .namesW-3), "...", sep="")
        cat(" ", snippet_name, sep="")
    }
    cat("\n")
}

### 'half_nrow' must be >= 1
.XStringSet.show_frame <- function(x, half_nrow=9L)
{
    lx <- length(x)
    iW <- nchar(as.character(lx)) + 2 # 2 for the brackets
    ncharMax <- max(nchar(x))
    widthW <- max(nchar(ncharMax), nchar("width"))
    .XStringSet.show_frame_header(iW, widthW, !is.null(names(x)))
    if (lx <= 2*half_nrow+1) {
        for (i in seq_len(lx))
            .XStringSet.show_frame_line(x, i, iW, widthW)
    } else {
        for (i in 1:half_nrow)
            .XStringSet.show_frame_line(x, i, iW, widthW)
        cat(format("...", width=iW, justify="right"),
            format("...", width=widthW, justify="right"),
            "...\n")
        for (i in (lx-half_nrow+1L):lx)
            .XStringSet.show_frame_line(x, i, iW, widthW)
    }
}

setMethod("show", "XStringSet",
    function(object)
    {
        cat("  A ", class(object), " instance of length ", length(object), "\n", sep="")
        if (length(object) != 0)
            .XStringSet.show_frame(object)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

### Supported 'i' types: numeric vector, logical vector, NULL and missing.
setMethod("[", "XStringSet",
    function(x, i, j, ..., drop)
    {
        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")
        if (missing(i))
            return(x)
        if (is.character(i))
            stop("cannot subset a ", class(x), " object by names")
        x@ranges <- x@ranges[i]
        x
    }
)

setMethod("rep", "XStringSet",
    function(x, times)
        x[rep.int(seq_len(length(x)), times)]
)

### Return an XString object of the same base type as 'x'.
### Example:
###   bs <- BString("ABCD-1234-abcd")
###   bset <- new("BStringSet", super=bs, start=1:8, width=2L*(7:0))
###   bset[[3]]
setMethod("[[", "XStringSet",
    function(x, i, j, ...)
    {
        i <- IRanges:::checkAndTranslateDbleBracketSubscript(x, i)
        start <- start(x@ranges)[i]
        width <- width(x@ranges)[i]
        subseq(super(x), start=start, width=width)
    }
)

setReplaceMethod("[[", "XStringSet",
    function(x, i, j,..., value)
    {
        stop("attempt to modify the value of a ", class(x), " instance")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "append" method.
###
### TODO: The current version performs too many copies of the sequence data!
###       This will change with the redesign of the XStringSet class (see
###       long comment at the beginning of this file).
###

setMethod("append", c("XStringSet", "XStringSet"),
    function(x, values, after=length(x))
    {
        basetype <- xsbasetype(x)
        if (xsbasetype(values) != basetype)
            stop("'x' and 'values' must be XStringSet objects of the same base type")
        if (!isSingleNumber(after))
            stop("'after' must be a single number")
        if (length(values) == 0)
            return(x)
        if (after != length(x))
            stop("'after != length(x)' is not supported for XStringSet objects, sorry!")
        cx <- compact(x)
        cvalues <- compact(values, basetype=basetype)
        ans_super <- XString.append(super(cx), super(cvalues))
        ans_start <- c(start(cx@ranges), start(cvalues@ranges) + length(super(cx)))
        ans_width <- c(width(cx), width(cvalues))
        ans_ranges <- new2("IRanges", start=ans_start, width=ans_width, check=FALSE)
        ans_names <- c(names(cx), names(cvalues))
        unsafe.newXStringSet(ans_super, ans_ranges, use.names=TRUE, names=ans_names)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Set Operations
###

.XStringSet.SetOperation <- function(x, y, FUN)
{
    basetype <- xsbasetype(x)
    if (xsbasetype(y) != basetype)
        stop("'x' and 'y' must be XStringSet objects of the same base type")
    XStringSet(basetype, FUN(as.character(unique(x)), as.character(unique(y))))
}

setMethod("union", c("XStringSet", "XStringSet"),
    function(x, y) .XStringSet.SetOperation(x, y, FUN = union)
)
setMethod("intersect", c("XStringSet", "XStringSet"),
    function(x, y) .XStringSet.SetOperation(x, y, FUN = intersect)
)
setMethod("setdiff", c("XStringSet", "XStringSet"),
    function(x, y) .XStringSet.SetOperation(x, y, FUN = setdiff)
)
setMethod("setequal", c("XStringSet", "XStringSet"),
    function(x, y) {
        basetype <- xsbasetype(x)
        if (basetype(y) != basetype)
            stop("'x' and 'y' must be XStringSet objects of the same base type")
        setequal(as.character(unique(x)), as.character(unique(y)))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "%in%" and match methods.
###

setMethod("%in%", c("character", "XStringSet"),
    function(x, table)
        XStringSet(xsbasetype(table), x) %in% table
)

setMethod("%in%", c("XString", "XStringSet"),
    function(x, table)
        XStringSet(xsbasetype(table), x) %in% table
)

setMethod("%in%", c("XStringSet", "XStringSet"),
    function(x, table)
        .Call("XStringSet_in_set", x, table, PACKAGE = "Biostrings")
)

setMethod("match", c("character", "XStringSet"),
    function (x, table, nomatch = NA_integer_, incomparables = NULL)
        match(XStringSet(xsbasetype(table), x), table, nomatch = nomatch,
              incomparables = incomparables)
)

setMethod("match", c("XString", "XStringSet"),
    function (x, table, nomatch = NA_integer_, incomparables = NULL)
        match(XStringSet(xsbasetype(table), x), table, nomatch = nomatch,
              incomparables = incomparables)
)

setMethod("match", c("XStringSet", "XStringSet"),
    function (x, table, nomatch = NA_integer_, incomparables = NULL) {
        if (!is.null(incomparables))
            stop("'incomparables' argument is not supported")
        if (!isSingleNumberOrNA(nomatch))
            stop("'nomatch' must be a single integer")
        nomatch <- as.integer(nomatch)
        .Call("XStringSet_match", x, table, nomatch, PACKAGE = "Biostrings")
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Equality.
###

### WON'T START THIS UNLESS SOMEONE HAS A USE CASE...
### Look at XStringViews-class.R for how to do this.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Other coercion methods.
###

### 'unlist(x)' turns XStringSet object 'x' into an XString object.
setMethod("unlist", "XStringSet",
    function(x, recursive=TRUE, use.names=TRUE)
        .Call("XStringSet_unlist", x, PACKAGE="Biostrings")
)

setMethod("as.character", "XStringSet",
    function(x, use.names=TRUE)
    {
        use.names <- normargUseNames(use.names)
        ans <- .Call("XStringSet_as_STRSXP",
                     x, xs_dec_lkup(x),
                     PACKAGE="Biostrings")
        if (use.names)
            names(ans) <- names(x)
        ans
    }
)

setMethod("toString", "XStringSet", function(x, ...) toString(as.character(x), ...))

setMethod("as.matrix", "XStringSet",
    function(x, use.names=TRUE)
    {
        use.names <- normargUseNames(use.names)
        nrow <- length(x)
        if (nrow == 0)
            stop("'x' must contain at least 1 string")
        widths <- width(x)
        ncol <- widths[1]
        if (!all(widths == ncol))
            stop("'x' strings are not equal-width")
        y <- as.character(x, use.names=FALSE)
        y <- unlist(strsplit(y, NULL), recursive=FALSE, use.names=FALSE)
        m <- matrix(y, nrow=nrow, byrow=TRUE)
        if (use.names)
            rownames(m) <- names(x)
        m
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Ordering and related methods.
###

setMethod("order", "XStringSet",
    function(..., na.last=TRUE, decreasing=FALSE)
    {
        if (!missing(na.last) && !isTRUE(na.last))
            warning("argument 'na.last' is ignored when ordering XStringSet objects")
        if (!isTRUEorFALSE(decreasing))
            stop("'decreasing' must be TRUE or FALSE")
        if (decreasing)
            stop("'decreasing=TRUE' is not supported yet, sorry!")
        args <- list(...)
        ## All the arguments are guaranteed to be XStringSet objects
        if (length(args) != 1)
            return(callNextMethod())
        .Call("XStringSet_order", args[[1]], PACKAGE="Biostrings")
    }
)

setMethod("sort", "XStringSet",
    function(x, decreasing=FALSE, ...)
    {
        if (!isTRUEorFALSE(decreasing))
            stop("'decreasing' must be TRUE or FALSE")
        if (decreasing)
            stop("'decreasing=TRUE' is not supported yet, sorry!")
        x[order(x)]
    }
)

setMethod("rank", "XStringSet",
    function(x, na.last = TRUE,
             ties.method = c("average", "first", "random", "max", "min"))
    {
         if (!missing(na.last) && !isTRUE(na.last))
             warning("argument 'na.last' is ignored when ordering XStringSet objects")
         if (!missing(ties.method) && ties.method != "min")
             stop("only the 'min' option to the 'ties.method' argument is supported")
         .Call("XStringSet_rank", x, PACKAGE="Biostrings")
     }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Duplicated and related methods.
###

setMethod("duplicated", "XStringSet",
    function(x, incomparables=FALSE, ...)
        .Call("XStringSet_duplicated", x, PACKAGE="Biostrings")
)

### Should be moved to IRanges and made the default method for Sequence objects
setMethod("unique", "XStringSet",
    function(x, incomparables=FALSE, ...)
        x[!duplicated(x, incomparables=incomparables)]
)

