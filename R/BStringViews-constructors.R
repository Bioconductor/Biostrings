### =========================================================================
### Constructor-like functions and generics for BStringViews objects
### -------------------------------------------------------------------------

### WARNING: This function is unsafe! (it doesn't check its arguments)
### Only 2 valid ways to use it:
###   new("BStringViews", subject)
###   new("BStringViews", subject, start, end)
### where 'subject' is a BString (or derived) object,
### and 'start' and 'end' are integer vectors of the same length
### such that 'start <= end'.
setMethod("initialize", "BStringViews",
    function(.Object, subject, start, end)
    {
        .Object@subject <- subject
        if (!missing(start))
            .Object@views <- data.frame(start=start, end=end)
        .Object
    }
)

### The 3 functions below share the following properties:
###   - They are exported (and safe).
###   - First argument is 'subject'. It must be a character vector or a BString
###     (or derived) object.
###   - Passing something else to 'subject' provokes an error.
###   - They return a BStringViews object whose 'subject' slot is the object
###     passed in the 'subject' argument.

.makeViews <- function(subject, start, end)
{
    if (!isLooseNumeric(start) || !isLooseNumeric(end))
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
    data.frame(start=start, end=end)
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
views <- function(subject, start=NA, end=NA)
{
    if (class(subject) == "character")
        subject <- BString(subject)
    ans <- new("BStringViews", subject)
    ans@views <- .makeViews(subject, start, end)
    ans
}

### 'width' is the vector of view widths.
### 'gapwidth' is the vector of inter-view widths (recycled).
adjacentViews <- function(subject, width, gapwidth=0)
{
    if (class(subject) == "character")
        subject <- BString(subject)
    ans <- new("BStringViews", subject)
    ONE <- as.integer(1)
    if (!is.numeric(width) || !isTRUE(all(width >= ONE))) # NA-proof
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
    start <- integer(lw)
    end <- integer(lw)
    start[ONE] <- ONE
    end[ONE] <- width[ONE]
    if (lw >= 2) {
        j <- ONE
        for (i in 2:lw) {
            start[i] <- end[i-ONE] + ONE + gapwidth[j]
            end[i] <- start[i] + width[i] - ONE
            if (j < lg) j <- j + ONE else j <- ONE
        }
    }
    ans@views <- data.frame(start=start, end=end)
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "Normalized" BStringViews objects
### 
### Definition: BStringViews object x is "normalized" means
###  (1) no "out of limits" views,
###  (2) views are sorted from left to right (start(x) is ascending),
###  (3) views don't overlap and can't be adjacents i.e. there is at least
###      1 letter between 2 any given views.
### If length(x) >= 2, then the 3 above conditions are equivalent to:
###   1 <= start(x)[i] <= end(x)[i] <
###        start(x)[i+1] <= end(x)[i+1] <= nchar(subject(x))
### for every 1 <= i < length(x).
### If length(x) == 0, then x is normalized.
### If length(x) == 1, then x is normalized <=> view is not "out of limits".
###
### "Normalizing" a BStringViews object:
### For any given BStringViews object x, let's call S(x) the subset of
### integers defined by:
###   (union of all [start(x)[i],end(x)[i]]) inter [1,nchar(subject(x))]
### We can see that there is a unique "normalized" BStringViews object x0
### (with same subject as x) such that S(x0) = S(x).
### So we can define the normalize() function: x -> x0 = normalize(x).
### An interesting property of this function is that, for any x,
### x0 = normalize(x) is the shortest BStringViews object (with same
### subject than x) such that S(x0) = S(x).
### It can also been shown that length(x0) <= (nchar(subject(x)) + 1) / 2.
###
### Some basic operations on "normalized" BStringViews objects with same
### subject (the results of these operations are "normalized" too):
###   x2 <- !x: (unary operation) S(x2) =  [1:nchar(subject(x))] - S(x)
###   x3 <- x | y: s(x3) = (S(x) union S(y))
###   x4 <- x & y: s(x4) = (S(x) inter S(y))
### So we end up having the equivalent of the fundamental set operations
### (complementary, union, intersection).
### No need to define an equivalent for the set difference (S(x) - S(y)),
### or the set symetric difference ((S(x) -S(y)) union (S(y) - S(x)))
### since they can be achieved by using the fundamental operations.
###
### Notes:
### - !x is actually obtained with mask(x) (see below)
### - The results of the | and & operations are undefined if one of
###   the operand is not "normalized" or if the 2 operands don't have the
###   same subject! (no need to do any check of any sort, or to try to fail
###   gracefully, this will just result in slowing down the operators).
###   The result has always the same subject as the operands.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "mask" generic and methods.
###
### For any given BStringViews object x, y <- mask(x) is the shortest (in
### term of number of views) BStringViews object such that:
###   (a) subject(y) == subject(x)
###   (b) y views cover the parts of the subject that are not covered by x
###       views
###   (c) y views are sorted from left to right
###   (d) y views are not "out of limits"
### Relationship with the concept of "normalized" BStringViews objects (see
### above):
###   - mask(x) is !x
###   - mask(x) is normalized.
###   - mask(mask(x)) is normalize(x).
###

BStringViews.mask <- function(x)
{
    ii <- order(start(x))
    new_start <- new_end <- integer(0)
    start0 <- 1L
    for (i in ii) {
        end0 <- start(x)[i] - 1L
        start1 <- end(x)[i] + 1L
        if (end0 >= start0) {
            new_start <- c(new_start, start0)
            new_end <- c(new_end, end0)
        }
        if (start1 > start0)
            start0 <- start1
    }
    if (start0 <= length(subject(x))) {
        new_start <- c(new_start, start0)
        new_end <- c(new_end, length(subject(x)))
    }
    x@views <- data.frame(start=new_start, end=new_end)
    x
}

setGeneric("mask", function(x, ...) standardGeneric("mask"))

setMethod("mask", "BString",
    function(x, start=NA, end=NA)
    {
        BStringViews.mask(views(x, start, end))
    }
)

setMethod("mask", "BStringViews",
    function(x)
    {
        BStringViews.mask(x)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "BStringViews" generic and methods.
###

setGeneric(
    "BStringViews", signature="src",
    function(src, subjectClass, sep="") standardGeneric("BStringViews")
)

### 'subjectClass' must be "BString" or one of its derivations ("DNAString",
### "RNAString" or "AAString").
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
    function(src, subjectClass, sep="")
    {
        if (!is.character(sep))
            sep <- toString(sep)
        collapsed <- paste(src, collapse=sep)
        subject <- new(subjectClass, collapsed)
        adjacentViews(subject, nchar(src), nchar(sep))
    }
)

### Only FASTA files are supported for now.
### Typical use:
###   file <- system.file("Exfiles", "someORF.fsa", package="Biostrings")
###   v <- BStringViews(file(file), "DNAString")
setMethod("BStringViews", "file",
    function(src, subjectClass, sep="")
    {
        fasta <- readFASTA(src)
        src <- sapply(fasta, function(rec) rec$seq)
        desc <- sapply(fasta, function(rec) rec$desc)
        ans <- BStringViews(src, subjectClass, sep)
        ans@views$desc <- desc
        ans
    }
)

### Called when 'src' is a BString (or derived) object.
### When not missing, 'subjectClass' must be "BString" or one of its
### derivations ("DNAString", "RNAString" or "AAString").
setMethod("BStringViews", "BString",
    function(src, subjectClass, sep="")
    {
        if (!missing(sep)) {
            ## The semantic is: views are delimited by the occurences of 'sep'
            ## in 'src' (a kind of strsplit() for BString objects).
            ## Uncomment when normalize() and ! method are ready (see TODO file):
            #return(!normalize(matchPattern(sep, b, fixed=TRUE)))
            stop("'sep' not yet supported when 'src' is a \"BString\" object")
        }
        if (!missing(subjectClass) && subjectClass != class(src))
            src <- new(subjectClass, src)
        new("BStringViews", src, as.integer(1), src@length)
    }
)

### Called when 'src' is a BStringViews object.
### 'subjectClass' must be "BString" or one of its derivations ("DNAString",
### "RNAString" or "AAString").
### The 'sep' arg is ignored.
setMethod("BStringViews", "BStringViews",
    function(src, subjectClass, sep="")
    {
        if (!missing(sep))
            stop("'sep' not supported when 'src' is a \"BStringViews\" object")
        src@subject <- new(subjectClass, src@subject)
        src
    }
)

setGeneric(
    "writeFASTA", signature="x",
    function(x, file, width=80) standardGeneric("writeFASTA")
)

setMethod("writeFASTA", "BStringViews",
    function(x, file, width=80)
    {
        y <- list()
        for (i in seq_len(length(x))) {
            y[[i]] <- list(desc=desc(x)[i], seq=as.character(x[[i]]))
        }
        writeFASTA(y, file, width)
    }
)

