### =========================================================================
### Constructor-like functions and generics for XStringViews objects
### -------------------------------------------------------------------------

### The "views" and "adjacentViews" functions below share the following
### properties:
###   - They are exported (and safe).
###   - First argument is 'subject'. It must be a character vector or an
###     XString object.
###   - Passing something else to 'subject' provokes an error.
###   - They return an XStringViews object whose 'subject' slot is the object
###     passed in the 'subject' argument.

.safeMakeViews <- function(subject, start, end)
{
    if (!isNumericOrNAs(start) || !isNumericOrNAs(end))
        stop("'start' and 'end' must be numerics")
    if (!is.integer(start))
        start <- as.integer(start)
    start[is.na(start)] <- as.integer(1)
    if (!is.integer(end))
        end <- as.integer(end)
    end[is.na(end)] <- subject@length
    if (length(start) < length(end))
        start <- recycleVector(start, length(end))
    else if (length(end) < length(start))
        end <- recycleVector(end, length(start))
    if (!all(start <= end))
        stop("'start' and 'end' must verify 'start <= end'")
    new("Views", start=start, width=end-start+1L, check=FALSE)
}

### Typical use:
###   dna <- DNAString("AA-CC-GG-TT")
### Just one view:
###   dnav1 <- views(dna, 2, 7)
### 9 views, 3 are out of limits:
###   dnav2 <- views(dna, 6:-2, 6:14)
### 5 out of limits views, all have a width of 6:
###   dnav3 <- views(dna, -5:-1, 0:4)
### Same as doing views(dna, 1, length(dna)):
###   dnav4 <- views(dna)
### An XStringViews object with no view:
###   dnav5 <- views(dna, integer(0), integer(0))
### TODO: Use same args for this function as for the "subviews" function i.e.:
###   - add a 'width' arg (default to NA)
###   - add a 'check.limits' arg (default to TRUE) for raising an error if
###     some views are "out of limits"
views <- function(subject, start=NA, end=NA)
{
    if (!is(subject, "XString"))
        subject <- BString(subject)
    ans <- new("XStringViews", subject, check=FALSE)
    ranges <- .safeMakeViews(subject(ans), start, end)
    update(ans, start=start(ranges), width=width(ranges))
}

### 'width' is the vector of view widths.
### 'gapwidth' is the vector of inter-view widths (recycled).
### TODO: Use intToAdjacentRanges() in adjacentViews(), this
### will be MUCH faster (intToAdjacentRanges() is written in C).
### But first, the 'gapwidth' arg must be added to it...
adjacentViews <- function(subject, width, gapwidth=0)
{
    ONE <- as.integer(1)
    if (!is(subject, "XString"))
        subject <- BString(subject)
    if (!is.numeric(width) || !isTRUE(all(width >= ONE))) # NA-proof
        stop("'width' must be numerics >= 1")
    if (!is.numeric(gapwidth) || !isTRUE(all(gapwidth >= 0))) # NA-proof
        stop("'gapwidth' must be numerics >= 0")
    lw <- length(width)
    if (lw == 0)
        return(new("XStringViews", subject, check=FALSE))
    if (!is.integer(width))
        width <- as.integer(width)
    lg <- length(gapwidth)
    if (lw > 1 && lg == 0)
        stop("'gapwidth' is empty")
    if (!is.integer(gapwidth))
        gapwidth <- as.integer(gapwidth)
    ans_start <- integer(lw)
    ans_start[ONE] <- ONE
    j <- ONE
    for (i in seq_len(lw-1)) {
        ans_start[i+ONE] <- ans_start[i] + width[i] + gapwidth[j]
        if (j < lg) j <- j + ONE else j <- ONE
    }
    new("XStringViews", subject, start=ans_start, width=width, check=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "XStringViews" generic and methods.
###
### 'subjectClass' must be the name of one of the direct XString subtypes
### i.e. "BString", "DNAString", "RNAString" or "AAString".
###
### Benchmarks:
###   n <- 40000
###   x <- sapply(1:n, function(i)
###                    paste(sample(DNA_ALPHABET, 250, replace=TRUE), collapse=""))
###   v <- XStringViews(x, "DNAString")
### Comparing XStringViews() speed vs "old" vectorized DNAString() speed:
###       n  XStringViews  "old" DNAString
###   -----  ------------  ---------------
###    5000        0.26 s           4.15 s
###   10000        0.51 s          16.29 s
###   20000        0.99 s          64.85 s
###   40000        1.69 s         488.43 s
### The quadratic behaviour of "old" DNAString() (in Biostrings 1) was first
### reported by Wolfgang.
###

setGeneric("XStringViews", signature="x",
    function(x, subjectClass, collapse="") standardGeneric("XStringViews")
)

### The main purpose of this "XStringViews" method is to work on a character
### vector. However we use "ANY" instead of "character" for the signature
### because we've found some weird objects that look very much like character
### vectors but break the dispatch mechanism.
### For example:
###   library(hgu95av2probe)
###   x <- hgu95av2probe$sequence
###   is.character(x) # TRUE
###   is(x, "character") # FALSE
### Welcome to R object model mess!
setMethod("XStringViews", "ANY",
    function(x, subjectClass, collapse="")
    {
        if (!is.character(collapse))
            collapse <- toString(collapse)
        seq <- paste(x, collapse=collapse)
        subject <- XString(subjectClass, seq)
        ans <- adjacentViews(subject, nchar(x), gapwidth=nchar(collapse))
        names(ans) <- names(x)
        ans
    }
)

setMethod("XStringViews", "XString",
    function(x, subjectClass, collapse="")
    {
        if (!missing(collapse)) {
            ## The semantic is: views are delimited by the occurences of
            ## 'collapse' in 'x' (a kind of strsplit() for XString objects).
            ## Uncomment when normalize() and ! method are ready (see TODO file):
            #return(!normalize(matchPattern(collapse, b, fixed=TRUE)))
            stop("'collapse' not yet supported when 'x' is an XString object")
        }
        if (!missing(subjectClass) && subjectClass != class(x))
            x <- XString(subjectClass, x)
        new("XStringViews", x, start=1L, width=nchar(x), check=FALSE)
    }
)

### The 'collapse' arg is ignored.
setMethod("XStringViews", "XStringViews",
    function(x, subjectClass, collapse="")
    {
        if (!missing(collapse))
            stop("'collapse' not supported when 'x' is an XStringViews object")
        x@subject <- XString(subjectClass, subject(x))
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "restrict" method and the "trim" function.
###

setMethod("restrict", "XStringViews",
    function(x, start, end, keep.all.ranges=FALSE, use.names=TRUE)
    {
        if (!missing(keep.all.ranges))
            stop("'keep.all.ranges' is not supported for XStringViews objects")
        callNextMethod(x, start, end, use.names=use.names)
    }
)

trim <- function(x, use.names=TRUE)
{
    if (!is(x, "XStringViews"))
        stop("'x' must be an XStringViews object")
    y <- restrict(x, 1L, nchar(subject(x)), use.names=use.names)
    if (length(y) != length(x))
        stop("some views are not overlapping with the subject, cannot trim them")
    y
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "narrow" method and the "subviews" function.
###

setMethod("narrow", "XStringViews",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE)
    {
        y <- callNextMethod()
        if (any(width(y) == 0))
            stop("some views would have a null width after narrowing")
        y
    }
)

subviews <- function(x, start=NA, end=NA, width=NA, use.names=TRUE)
{
    if (!is(x, "XStringViews"))
        stop("'x' must be an XStringViews object")
    narrow(x, start=start, end=end, width=width, use.names=use.names)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Deprecated generics and methods.
###

setGeneric("BStringViews", signature="src",
    function(src, subjectClass, collapse="") standardGeneric("BStringViews")
)

.redirect.BStringViews.to.XStringViews <- function(src, subjectClass, collapse="")
{
    .Deprecated("XStringViews")
    XStringViews(src, subjectClass, collapse=collapse)
}

setMethod("BStringViews", "ANY", .redirect.BStringViews.to.XStringViews)

setMethod("BStringViews", "file",
    function(src, subjectClass, collapse="")
    {
        .Deprecated("read.XStringViews")
        read.XStringViews(src, "fasta", subjectClass, collapse=collapse)
    }
)

setMethod("BStringViews", "XString", .redirect.BStringViews.to.XStringViews)

setMethod("BStringViews", "XStringViews", .redirect.BStringViews.to.XStringViews)

