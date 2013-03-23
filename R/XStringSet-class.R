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
    ans_class <- paste(class(xvector), "Set", sep="")
    ans <- IRanges:::unsafe.newXVectorList1(ans_class, xvector, ranges)
    if (normargUseNames(use.names))
        names(ans) <- names
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "seqtype" and "seqtype<-" methods.
###

setMethod("seqtype", "XStringSet",
    function(x) seqtype(new(elementType(x)))
)

### NOT an endomorphism in general! (Because it downgrades 'x' to a
### B/DNA/RNA/AAStringSet instance.)
setReplaceMethod("seqtype", "XStringSet",
    function(x, value)
    {
        lkup <- get_seqtype_conversion_lookup(seqtype(x), value)
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
.oneSeqToXStringSet <- function(seqtype, x, start, end, width, use.names)
{
    ans_xvector <- XString(seqtype, x)
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

.charToXStringSet <- function(seqtype, x, start, end, width, use.names)
{
    if (length(x) == 1L) {
        ans <- .oneSeqToXStringSet(seqtype, x, start, end, width, use.names)
        return(ans)
    }
    use.names <- normargUseNames(use.names)
    ans_elementType <- paste(seqtype, "String", sep="")
    ans_class <- paste(ans_elementType, "Set", sep="")
    solved_SEW <- solveUserSEW(width(x), start=start, end=end, width=width)
    ans <- .Call2("new_XStringSet_from_CHARACTER",
                 ans_class, ans_elementType,
                 x, start(solved_SEW), width(solved_SEW),
                 get_seqtype_conversion_lookup("B", seqtype),
                 PACKAGE="Biostrings")
    if (use.names)
        names(ans) <- names(x)
    ans
}


setGeneric("XStringSet", signature="x",
    function(seqtype, x, start=NA, end=NA, width=NA, use.names=TRUE)
        standardGeneric("XStringSet")
)

setMethod("XStringSet", "character",
    function(seqtype, x, start=NA, end=NA, width=NA, use.names=TRUE)
    {
        if (is.null(seqtype))
            seqtype <- "B"
        .charToXStringSet(seqtype, x, start, end, width, use.names)
    }
)

setMethod("XStringSet", "factor",
    function(seqtype, x, start=NA, end=NA, width=NA, use.names=TRUE)
    {
        if (is.null(seqtype))
            seqtype <- "B"
        if (length(x) < nlevels(x)) {
            ans <- .charToXStringSet(seqtype, as.character(x),
                                     start, end, width, use.names)
            return(ans)
        }
        ## If 'x' has less levels than elements, then it's cheaper to
        ## operate on its levels. In case of equality (i.e. if
        ## length(x) == nlevels(x)), the price is the same but the final
        ## XStringSet object obtained by operating on the levels might use
        ## less memory (if 'x' contains duplicated values).
        ans <- .charToXStringSet(seqtype, levels(x),
                                 start, end, width, use.names)
        ans[as.integer(x)]
    }
)

setMethod("XStringSet", "XString",
    function(seqtype, x, start=NA, end=NA, width=NA, use.names=TRUE)
        .oneSeqToXStringSet(seqtype, x, start, end, width, use.names)
)

setMethod("XStringSet", "XStringSet",
    function(seqtype, x, start=NA, end=NA, width=NA, use.names=TRUE)
    {
        ans <- narrow(x, start=start, end=end, width=width, use.names=use.names)
        ## `seqtype<-` must be called even when 'seqtype' is NULL
        ## because we want to enforce downgrade to a B/DNA/RNA/AAStringSet
        ## instance
        if (is.null(seqtype))
            seqtype <- seqtype(x)
        seqtype(ans) <- seqtype
        ans
    }
)

setMethod("XStringSet", "list",
    function(seqtype, x, start=NA, end=NA, width=NA, use.names=TRUE)
    {
        x_len <- length(x)
        if (x_len == 0L) {
            tmp_elementType <- "BString"
        } else {
            tmp_elementType <- paste(seqtype(x[[1L]]), "String", sep="")
        }
        tmp_class <- paste(tmp_elementType, "Set", sep="")
        tmp <- IRanges:::new_XVectorList_from_list_of_XVector(tmp_class, x)
        XStringSet(seqtype, tmp,
                   start=start, end=end, width=width, use.names=use.names)
    }
)

### 2 extra "XStringSet" methods to deal with the probe sequences stored
### in the *probe annotation packages (e.g. drosophila2probe).

setMethod("XStringSet", "AsIs",
    function(seqtype, x, start=NA, end=NA, width=NA, use.names=TRUE)
    {
        if (!is.character(x))
            stop("unsuported input type")
        class(x) <- "character" # keeps the names (unlike as.character())
	.charToXStringSet(seqtype, x, start, end, width, use.names)
    }
)

setMethod("XStringSet", "probetable",
    function(seqtype, x, start=NA, end=NA, width=NA, use.names=TRUE)
        XStringSet(seqtype, x$sequence, start=start, end=end, width=width,
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
    function(from) {seqtype(from) <- "B"; from}
)
setAs("XStringSet", "DNAStringSet",
    function(from) {seqtype(from) <- "DNA"; from}
)
setAs("XStringSet", "RNAStringSet",
    function(from) {seqtype(from) <- "RNA"; from}
)
setAs("XStringSet", "AAStringSet",
    function(from) {seqtype(from) <- "AA"; from}
)

setAs("character", "BStringSet", function(from) BStringSet(from))
setAs("character", "DNAStringSet", function(from) DNAStringSet(from))
setAs("character", "RNAStringSet", function(from) RNAStringSet(from))
setAs("character", "AAStringSet", function(from) AAStringSet(from))
setAs("character", "XStringSet", function(from) BStringSet(from))

setAs("factor", "BStringSet", function(from) BStringSet(from))
setAs("factor", "DNAStringSet", function(from) DNAStringSet(from))
setAs("factor", "RNAStringSet", function(from) RNAStringSet(from))
setAs("factor", "AAStringSet", function(from) AAStringSet(from))
setAs("factor", "XStringSet", function(from) BStringSet(from))

setAs("XString", "BStringSet", function(from) BStringSet(from))
setAs("XString", "DNAStringSet", function(from) DNAStringSet(from))
setAs("XString", "RNAStringSet", function(from) RNAStringSet(from))
setAs("XString", "AAStringSet", function(from) AAStringSet(from))
setAs("XString", "XStringSet", function(from) XStringSet(seqtype(from), from))

setAs("list", "BStringSet", function(from) BStringSet(from))
setAs("list", "DNAStringSet", function(from) DNAStringSet(from))
setAs("list", "RNAStringSet", function(from) RNAStringSet(from))
setAs("list", "AAStringSet", function(from) AAStringSet(from))


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
### Splitting
###

setMethod("splitAsListReturnedClass", "BStringSet",
    function(x) "BStringSetList"
)
setMethod("splitAsListReturnedClass", "DNAStringSet",
    function(x) "DNAStringSetList"
)
setMethod("splitAsListReturnedClass", "RNAStringSet",
    function(x) "RNAStringSetList"
)
setMethod("splitAsListReturnedClass", "AAStringSet",
    function(x) "AAStringSetList"
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Set Operations
###

.XStringSet.SetOperation <- function(x, y, FUN)
{
    x_seqtype <- seqtype(x)
    if (seqtype(y) != x_seqtype)
        stop("'x' and 'y' must be XStringSet objects containing ",
             "sequences of the same type")
    XStringSet(x_seqtype, FUN(as.character(unique(x)), as.character(unique(y))))
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
    function(x, y)
    {
        x_seqtype <- seqtype(x)
        if (seqtype(y) != x_seqtype)
            stop("'x' and 'y' must be XStringSet objects containing ",
                 "sequences of the same type")
        setequal(as.character(unique(x)), as.character(unique(y)))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Other coercion methods.
###

### 'unlist(x)' turns XStringSet object 'x' into an XString object.
setMethod("unlist", "XStringSet",
    function(x, recursive=TRUE, use.names=TRUE)
        .Call2("XStringSet_unlist", x, PACKAGE="Biostrings")
)

setMethod("as.character", "XStringSet",
    function(x, use.names=TRUE)
    {
        use.names <- normargUseNames(use.names)
        ans <- .Call2("new_CHARACTER_from_XStringSet",
                     x, xs_dec_lkup(x),
                     PACKAGE="Biostrings")
        if (use.names)
            names(ans) <- names(x)
        ans
    }
)

setMethod("as.vector", "XStringSet",
    function(x, mode="any")
    {
        if (!isSingleString(mode))
            stop("'mode' must be a single string")
        if (!(mode %in% c("any", "character")))
            stop("'mode' can only be \"any\" or \"character\" ",
                 "when 'x' is an XStringSet object")
        as.character(x)
    }
)

setMethod("toString", "XStringSet",
    function(x, ...) toString(as.character(x), ...)
)

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

