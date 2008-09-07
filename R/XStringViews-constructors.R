### =========================================================================
### Constructor-like functions and generics for XStringViews objects
### -------------------------------------------------------------------------

### 'width' is the vector of view widths.
### 'gapwidth' is the vector of inter-view widths (recycled).
### TODO: Use intToAdjacentRanges() in adjacentViews(), this
### will be MUCH faster (intToAdjacentRanges() is written in C).
### But first, the 'gapwidth' arg must be added to it...
adjacentViews <- function(subject, width, gapwidth=0)
{
    ONE <- as.integer(1)
    if (!is(subject, "XString"))
        subject <- XString(NULL, subject)
    if (!is.numeric(width) || !isTRUE(all(width >= ONE))) # NA-proof
        stop("'width' must be numerics >= 1")
    if (!is.numeric(gapwidth) || !isTRUE(all(gapwidth >= 0))) # NA-proof
        stop("'gapwidth' must be numerics >= 0")
    lw <- length(width)
    if (lw == 0)
        return(new("XStringViews", subject))
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
    new("XStringViews", subject, start=ans_start, end=ans_start+width-1)
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
        new("XStringViews", x, start=1L, end=nchar(x))
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

