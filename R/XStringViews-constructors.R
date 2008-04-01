### =========================================================================
### Constructor-like functions and generics for BStringViews objects
### -------------------------------------------------------------------------

### The "views" and "adjacentViews" functions below share the following
### properties:
###   - They are exported (and safe).
###   - First argument is 'subject'. It must be a character vector or an
###     XString object.
###   - Passing something else to 'subject' provokes an error.
###   - They return a BStringViews object whose 'subject' slot is the object
###     passed in the 'subject' argument.

.safeMakeIRanges <- function(subject, start, end)
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
    ## The NA-proof version of 'if (any(end < start))'
    if (!isTRUE(all(start <= end)))
        stop("'start' and 'end' must verify 'start <= end'")
    new("IRanges", start=start, width=end-start+1L)
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
### A BStringViews object with no view:
###   dnav5 <- views(dna, integer(0), integer(0))
### TODO: Use same args for this function as for the "subviews" function i.e.:
###   - add a 'width' arg (default to NA)
###   - add a 'check.limits' arg (default to TRUE) for raising an error if
###     some views are "out of limits"
views <- function(subject, start=NA, end=NA)
{
    if (!is(subject, "XString"))
        subject <- BString(subject)
    ans <- new("BStringViews", subject, check=FALSE)
    ranges <- .safeMakeIRanges(subject(ans), start, end)
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
    if (!is.numeric(width) || !isTRUE(all(width >= ONE))) # NA-proof
        stop("'width' must be numerics >= 1")
    if (!is.numeric(gapwidth) || !isTRUE(all(gapwidth >= 0))) # NA-proof
        stop("'gapwidth' must be numerics >= 0")
    lw <- length(width)
    if (lw == 0)
        return(new("BStringViews", subject, check=FALSE))
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
    new("BStringViews", subject, start=ans_start, width=width, check=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "BStringViews" generic and methods.
###

setGeneric("BStringViews", signature="src",
    function(src, subjectClass, collapse="") standardGeneric("BStringViews")
)

### 'subjectClass' must be the name of one of the direct XString subtypes i.e.
### "BString", "DNAString", "RNAString" or "AAString".
###
### Benchmarks:
###   n <- 40000
###   src <- sapply(1:n, function(i)
###                      paste(sample(DNA_ALPHABET, 250, replace=TRUE), collapse=""))
###   v <- BStringViews(src, "DNAString")
### Comparing BStringViews() speed vs "old" vectorized DNAString() speed:
###       n  BStringViews  "old" DNAString
###   -----  ------------  ---------------
###    5000        0.26 s           4.15 s
###   10000        0.51 s          16.29 s
###   20000        0.99 s          64.85 s
###   40000        1.69 s         488.43 s
### The quadratic behaviour of "old" DNAString() was first reported
### by Wolfgang.

### This is the "BStringViews" for "character" but we use the "ANY" signature
### anyway. This because we've found some weird objects that look very much
### like "character" vectors but break the dispatch mechanism.
### For example:
###   library(hgu95av2probe)
###   src <- hgu95av2probe$sequence
###   is.character(src) # TRUE
###   is(src, "character") # FALSE
### Welcome to R object model mess!
setMethod("BStringViews", "ANY",
    function(src, subjectClass, collapse="")
    {
        if (!is.character(collapse))
            collapse <- toString(collapse)
        seq <- paste(src, collapse=collapse)
        subject <- XString(subjectClass, seq)
        ans <- adjacentViews(subject, nchar(src), nchar(collapse))
        desc(ans) <- names(src)
        ans
    }
)

setMethod("BStringViews", "file",
    function(src, subjectClass, collapse="")
    {
        .Deprecated("read.BStringViews")
        read.BStringViews(src, "fasta", subjectClass, collapse)
    }
)

### Called when 'src' is an XString object.
### When not missing, 'subjectClass' must be the name of one of the direct
### XString subtypes i.e. "BString", "DNAString", "RNAString" or "AAString".
setMethod("BStringViews", "XString",
    function(src, subjectClass, collapse="")
    {
        if (!missing(collapse)) {
            ## The semantic is: views are delimited by the occurences of
            ## 'collapse' in 'src' (a kind of strsplit() for XString objects).
            ## Uncomment when normalize() and ! method are ready (see TODO file):
            #return(!normalize(matchPattern(collapse, b, fixed=TRUE)))
            stop("'collapse' not yet supported when 'src' is an XString object")
        }
        if (!missing(subjectClass) && subjectClass != class(src))
            src <- XString(subjectClass, src)
        new("BStringViews", src, start=1L, width=nchar(src), check=FALSE)
    }
)

### Called when 'src' is a BStringViews object.
### 'subjectClass' must be the name of one of the direct XString subtypes i.e.
### "BString", "DNAString", "RNAString" or "AAString".
### The 'collapse' arg is ignored.
setMethod("BStringViews", "BStringViews",
    function(src, subjectClass, collapse="")
    {
        if (!missing(collapse))
            stop("'collapse' not supported when 'src' is a BStringViews object")
        src@subject <- XString(subjectClass, subject(src))
        src
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "restrict" method and the "trim" function.
###

setMethod("restrict", "BStringViews",
    function(x, start, end, keep.all.ranges=FALSE, use.names=TRUE)
    {
        if (!missing(keep.all.ranges))
            stop("'keep.all.ranges' is not supported for BStringViews objects")
        callNextMethod(x, start, end, use.names=use.names)
    }
)

trim <- function(x, use.names=TRUE)
{
    if (!is(x, "BStringViews"))
        stop("'x' must be a BStringViews object")
    y <- restrict(x, 1L, nchar(subject(x)), use.names=use.names)
    if (length(y) != length(x))
        stop("some views are not overlapping with the subject, cannot trim them")
    y
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "narrow" method and the "subviews" function.
###

setMethod("narrow", "BStringViews",
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
    if (!is(x, "BStringViews"))
        stop("'x' must be a BStringViews object")
    narrow(x, start=start, end=end, width=width, use.names=use.names)
}

