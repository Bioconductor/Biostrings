### =========================================================================
### XStringSet objects
### -------------------------------------------------------------------------
###

setClass("XStringSet",
    contains="XRawList",
    representation("VIRTUAL"),
    prototype(elementType="XString")
)

### The concrete XStringSet subclasses below have no additional slots.
setClass("BStringSet",
    contains="XStringSet",
    representation(),
    prototype(elementType="BString")
)
setClass("DNAStringSet",
    contains="XStringSet",
    representation(),
    prototype(elementType="DNAString")
)
setClass("RNAStringSet",
    contains="XStringSet",
    representation(),
    prototype(elementType="RNAString")
)
setClass("AAStringSet",
    contains="XStringSet",
    representation(),
    prototype(elementType="AAString")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor-like methods.
###

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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "subseq" endomorphism and related transformations.
###
### Methods for XStringSet objects are inherited from the XVectorList class.
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

setMethod("subseq", "character",
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

setReplaceMethod("subseq", "character",
    function(x, start=NA, end=NA, width=NA, value)
    {
        bands <- threebands(x, start=start, end=end, width=width)
        paste(bands$left, value, bands$right, sep="")
    }
)

### TODO: Make this a method for XVectorList objects and move it to the
### IRanges package (this means the implementation cannot use xscat() anymore).
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
### to contain valid ranges on 'xvector' i.e. ranges that are within the
### limits of 'xvector'.
###

unsafe.newXStringSet <- function(xvector, ranges, use.names=FALSE, names=NULL)
{
    classname <- paste(class(xvector), "Set", sep="")
    ans <- IRanges:::unsafe.newXVectorList1(classname, xvector, ranges)
    if (normargUseNames(use.names))
        names(ans) <- names
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "xsbasetype" and "xsbasetype<-" methods.
###

setMethod("xsbasetype", "XStringSet",
    function(x) xsbasetype(newEmptyXString(elementType(x)))
)

### NOT an endomorphism in general! (Because it downgrades 'x' to a
### B/DNA/RNA/AAStringSet instance.)
setReplaceMethod("xsbasetype", "XStringSet",
    function(x, value)
    {
        lkup <- get_xsbasetypes_conversion_lookup(xsbasetype(x), value)
        if (!is.null(lkup))
            x <- xvcopy(x, lkup=lkup)  # temporarily breaks 'x'!
        ans_class <- paste(value, "StringSet", sep="")
        new2(ans_class, pool=x@pool, ranges=x@ranges, check=FALSE)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The XStringSet() constructor. NOT exported.
###
### This constructor and its helper functions use the uSEW (user-specified
### Start/End/Width) interface.
###

### 'x' must be a character string or an XString object.
.oneSeqToXStringSet <- function(basetype, x, start, end, width, use.names)
{
    ans_xvector <- XString(basetype, x)
    ans_ranges <- solveUserSEW(length(ans_xvector),
                               start=start, end=end, width=width,
                               rep.refwidths=TRUE)
    ## We mimic how substring() replicates the name of a single string (try
    ## 'substring(c(A="abcdefghij"), 2, 6:2)').
    if (!is(x, "XString") && normargUseNames(use.names)) {
        ans_names <- names(x)
        if (!is.null(ans_names))
            ans_names <- rep.int(ans_names, length(ans_ranges))
    } else {
        ans_names <- NULL
    }
    unsafe.newXStringSet(ans_xvector, ans_ranges,
                         use.names=TRUE, names=ans_names)
}

.charToXStringSet <- function(basetype, x, start, end, width, use.names)
{
    if (length(x) == 1L) {
        ans <- .oneSeqToXStringSet(basetype, x, start, end, width, use.names)
        return(ans)
    }
    use.names <- normargUseNames(use.names)
    elementType <- paste(basetype, "String", sep="")
    classname <- paste(elementType, "Set", sep="")
    solved_SEW <- solveUserSEW(width(x), start=start, end=end, width=width)
    ans <- .Call("new_XStringSet_from_CHARACTER",
                 classname, elementType,
                 x, start(solved_SEW), width(solved_SEW),
                 get_xsbasetypes_conversion_lookup("B", basetype),
                 PACKAGE="Biostrings")
    if (use.names)
        names(ans) <- names(x)
    ans
}


setGeneric("XStringSet", signature="x",
    function(basetype, x, start=NA, end=NA, width=NA, use.names=TRUE)
        standardGeneric("XStringSet")
)

setMethod("XStringSet", "factor",
    function(basetype, x, start=NA, end=NA, width=NA, use.names=TRUE)
    {
        if (is.null(basetype))
            basetype <- "B"
        ans <- .charToXStringSet(basetype, levels(x),
                                 start, end, width, use.names)
        ans[as.integer(x)]
    }
)

setMethod("XStringSet", "character",
    function(basetype, x, start=NA, end=NA, width=NA, use.names=TRUE)
    {
        if (is.null(basetype))
            basetype <- "B"
        .charToXStringSet(basetype, x, start, end, width, use.names)
    }
)

setMethod("XStringSet", "XString",
    function(basetype, x, start=NA, end=NA, width=NA, use.names=TRUE)
        .oneSeqToXStringSet(basetype, x, start, end, width, use.names)
)

setMethod("XStringSet", "XStringSet",
    function(basetype, x, start=NA, end=NA, width=NA, use.names=TRUE)
    {
        ans <- narrow(x, start=start, end=end, width=width, use.names=use.names)
        ## `xsbasetype<-` must be called even when 'basetype' is NULL
        ## because we want to enforce downgrade to a B/DNA/RNA/AAStringSet
        ## instance
        if (is.null(basetype))
            basetype <- xsbasetype(x)
        xsbasetype(ans) <- basetype
        ans
    }
)

### 2 extra "XStringSet" methods to deal with the probe sequences stored
### in the *probe annotation packages (e.g. drosophila2probe).

setMethod("XStringSet", "AsIs",
    function(basetype, x, start=NA, end=NA, width=NA, use.names=TRUE)
    {
        if (!is.character(x))
            stop("unsuported input type")
        class(x) <- "character" # keeps the names (unlike as.character())
	.charToXStringSet(basetype, x, start, end, width, use.names)
    }
)

setMethod("XStringSet", "probetable",
    function(basetype, x, start=NA, end=NA, width=NA, use.names=TRUE)
        XStringSet(basetype, x$sequence, start=start, end=end, width=width,
                   use.names=use.names)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The user interfaces to the XStringSet() constructor.
###

BStringSet <- function(x=character(), start=NA, end=NA, width=NA,
                       use.names=TRUE)
    XStringSet("B", x, start=start, end=end, width=width,
                    use.names=use.names)

DNAStringSet <- function(x=character(), start=NA, end=NA, width=NA,
                         use.names=TRUE)
    XStringSet("DNA", x, start=start, end=end, width=width,
                      use.names=use.names)

RNAStringSet <- function(x=character(), start=NA, end=NA, width=NA,
                         use.names=TRUE)
    XStringSet("RNA", x, start=start, end=end, width=width,
                      use.names=use.names)

AAStringSet <- function(x=character(), start=NA, end=NA, width=NA,
                        use.names=TRUE)
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
### unsplit.list.of.XStringSet()
###
### Helper function, not for the end user.
###

unsplit.list.of.XStringSet <- function(class, value, f)
{
    ans <- rep.int(as("", class), length(f))
    unlisted_value <- do.call(c, unname(value))
    idx <- unname(split(seq_len(length(f)), f))
    ans[unlist(idx)] <- unlisted_value
    ans
}


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
    function(x, y, ...) .XStringSet.SetOperation(x, y, FUN = union)
)
setMethod("intersect", c("XStringSet", "XStringSet"),
    function(x, y, ...) .XStringSet.SetOperation(x, y, FUN = intersect)
)
setMethod("setdiff", c("XStringSet", "XStringSet"),
    function(x, y, ...) .XStringSet.SetOperation(x, y, FUN = setdiff)
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
        (match(x, table, nomatch = 0L) > 0L)
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
        if (xsbasetype(x) != xsbasetype(table))
            stop("'x' and 'table' must be XStringSet objects of the same base type")
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
        ans <- .Call("new_CHARACTER_from_XStringSet",
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

setMethod("is.unsorted", "XStringSet",
    function (x, na.rm = FALSE, strictly = FALSE) 
    {
        if (!is.logical(strictly) || length(strictly) != 1 || is.na(strictly))
            stop("'strictly' must be TRUE or FALSE")
        .Call("XStringSet_is_unsorted", x, strictly, PACKAGE="Biostrings")
    }
)

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

### Should be moved to IRanges and made the default method for Vector objects
setMethod("unique", "XStringSet",
    function(x, incomparables=FALSE, ...)
        x[!duplicated(x, incomparables=incomparables)]
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### updateObject()
###

### Update XStringSet objects created before the big change to the XStringSet
### internals ("super" slot replaced by "pool" slot).
### This change happened in Biostrings 2.13.43.
setMethod("updateObject", "XStringSet",
    function(object, ..., verbose=FALSE)
    {
        if (!is(try(object@pool, silent=TRUE), "try-error"))
            return(object)
        ans_xvector <- updateObject(object@super)
        ans_ranges <- updateObject(object@ranges)
        unsafe.newXStringSet(ans_xvector, ans_ranges,
                             use.names=TRUE, names=names(object))
    }
)

