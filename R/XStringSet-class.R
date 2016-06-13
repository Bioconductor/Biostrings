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
        ans_class <- paste(value, "StringSet", sep="")
        ## Don't try to replace the code below with 'as(x, ans_class)' because
        ## that would introduce a chicken-egg situation ('as(x, ans_class)'
        ## actually calls the seqtype() setter when 'x' is an XStringSet
        ## object).
        lkup <- get_seqtype_conversion_lookup(seqtype(x), value)
        if (!is.null(lkup))
            x <- xvcopy(x, lkup=lkup)  # temporarily breaks 'x'!
        new2(ans_class, pool=x@pool, ranges=x@ranges, check=FALSE)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Going from XString to XStringSet with extractList() and family.
###

setMethod("relistToClass", "XString",
    function(x) paste0(class(x), "Set")
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
        ## Use x[] <- ... to preserve names and any other attribute.
        x[] <- paste(bands$left, value, bands$right, sep="")
        x
    }
)

### TODO: Make this a method for XVectorList objects and move it to the
### IRanges package (this means the implementation cannot use xscat() anymore).
setReplaceMethod("subseq", "XStringSet",
    function(x, start=NA, end=NA, width=NA, value)
    {
        bands <- threebands(x, start=start, end=end, width=width)
        ## Use x[] <- ... to preserve class (endomorphism), names, metadata
        ## columns, and any other attribute.
        if (is.null(value)) {
            x[] <- xscat(bands$left, bands$right)
        } else {
            x[] <- xscat(bands$left, value, bands$right)
        }
        x
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
        x_names <- names(x)
        if (!is.null(x_names)) {
            ans_names <- rep.int(x_names, length(ans_ranges))
            names(ans_ranges) <- ans_names
        }
    }
    extractList(ans_xvector, ans_ranges)
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
        tmp <- XVector:::new_XVectorList_from_list_of_XVector(tmp_class, x)
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
        XStringSet(seqtype, x$sequence,
                   start=start, end=end, width=width, use.names=use.names)
)

### Default method.
setMethod("XStringSet", "ANY",
    function(seqtype, x, start=NA, end=NA, width=NA, use.names=TRUE)
        XStringSet(seqtype, as.character(x),
                   start=start, end=end, width=width, use.names=use.names)
)

setMethod("XStringSet", "missing",
    function(seqtype, x, start=NA, end=NA, width=NA, use.names=TRUE)
        XStringSet(seqtype, NULL)
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

setAs("ANY", "BStringSet", function(from) BStringSet(from))

setAs("ANY", "DNAStringSet", function(from) DNAStringSet(from))

setAs("ANY", "RNAStringSet", function(from) RNAStringSet(from))

setAs("ANY", "AAStringSet", function(from) AAStringSet(from))

setAs("ANY", "XStringSet",
    function(from)
    {
        from_seqtype <- try(seqtype(from), silent=TRUE)
        if (is(from_seqtype, "try-error"))
            from_seqtype <- "B"
        XStringSet(from_seqtype, from)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" and "showAsCell" methods.
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
.XStringSet.show_frame <- function(x, half_nrow=5L)
{
    if (is.null(head_nrow <- getOption("showHeadLines")))
        head_nrow <- half_nrow 
    if (is.null(tail_nrow <- getOption("showTailLines")))
        tail_nrow <- half_nrow
 
    lx <- length(x)
    iW <- nchar(as.character(lx)) + 2 # 2 for the brackets
    ncharMax <- max(nchar(x))
    widthW <- max(nchar(ncharMax), nchar("width"))
    .XStringSet.show_frame_header(iW, widthW, !is.null(names(x)))
    if (lx < (2*half_nrow+1L) | (lx < (head_nrow+tail_nrow+1L))) {
        for (i in seq_len(lx))
            .XStringSet.show_frame_line(x, i, iW, widthW)
    } else {
        if (head_nrow > 0)
            for (i in 1:head_nrow)
                .XStringSet.show_frame_line(x, i, iW, widthW)
        cat(format("...", width=iW, justify="right"),
            format("...", width=widthW, justify="right"),
            "...\n")
        if (tail_nrow > 0)
            for (i in (lx-tail_nrow+1L):lx)
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

setMethod("showAsCell", "XStringSet",
    function(object) 
        vapply(object, toSeqSnippet, character(1), width=23L)
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

setMethod("as.factor", "XStringSet",
    function(x)
    {
        as.factor(as.character(x))
    })

### TODO: Turn this into an S3/S4 combo for as.data.frame.XStringSet
setMethod("as.data.frame", "XStringSet",
    function(x, row.names=NULL, optional=FALSE)
    {
        x <- as.character(x)
        as.data.frame(x, row.names=NULL, optional=optional,
                         stringsAsFactors=FALSE)
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
        names(ans_ranges) <- names(object)
        extractList(ans_xvector, ans_ranges)
    }
)

