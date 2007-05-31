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

### The 2 functions above share the following properties:
###   - They are exported (and safe).
###   - First argument is 'subject'. It must be a character vector or a BString
###     (or derived) object.
###   - Passing something else to 'subject' provokes an error.
###   - They return a BStringViews object whose 'subject' slot is the object
###     passed in the 'subject' argument.

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
    ## Integrity checking
    if (!isLooseNumeric(start) || !isLooseNumeric(end))
        stop("'start' and 'end' must be numerics")
    #if (length(start) != length(end))
    #    stop("'start' and 'end' must have the same length")
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
    ans@views <- data.frame(start=start, end=end)
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


### The BStringViews() generic function
### -----------------------------------
### 'src' should typically be a character vector but the function will also
### work for other kind of input like numeric or even logical vectors.
### 'subjectClass' must be "BString" or one of its derivations ("DNAString",
### "RNAString" or "AAString").
###
### Benchmarks:
###   n <- 40000
###   src <- sapply(1:n, function(i) {paste(sample(DNA_ALPHABET, 250, replace=TRUE), collapse="")})
###   v <- BStringViews(src, "DNAString")
### Comparing BStringViews() speed vs "old" vectorized DNAString() speed:
###       n  BStringViews  "old" DNAString
###   -----  ------------  ---------------
###    5000        0.26 s           4.15 s
###   10000        0.51 s          16.29 s
###   20000        0.99 s          64.85 s
###   40000        1.69 s         488.43 s
### The quadratic behaviour of "old" DNAString() was start reported
### by Wolfgang.
BStringViews <- function(src, subjectClass, sep="")
{
    if (!is.character(sep))
        sep <- toString(sep)
    collapsed <- paste(src, collapse=sep)
    subject <- new(subjectClass, collapsed)
    adjacentViews(subject, nchar(src), nchar(sep))
}

setGeneric(
    "BStringViews",
    function(src, subjectClass, sep="") standardGeneric("BStringViews")
)

### Only FASTA files are supported for now.
### Typical use:
###   file <- system.file("Exfiles", "someORF.fsa", package="Biostrings")
###   v <- BStringViews(file(file), "DNAString")
setMethod("BStringViews", "file",
    function(src, subjectClass, sep)
    {
        fasta <- readFASTA(src)
        src <- sapply(fasta, function(rec) rec$seq)
        desc <- sapply(fasta, function(rec) rec$desc)
        ans <- BStringViews(src, subjectClass, sep) # call the default method
        ans@views$desc <- desc
        ans
    }
)

### Called when 'src' is a BString (or derived) object.
### When not missing, 'subjectClass' must be "BString" or one of its
### derivations ("DNAString", "RNAString" or "AAString").
setMethod("BStringViews", "BString",
    function(src, subjectClass, sep)
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
    function(src, subjectClass, sep)
    {
        if (!missing(sep))
            stop("'sep' not supported when 'src' is a \"BStringViews\" object")
        src@subject <- new(subjectClass, src@subject)
        src
    }
)

setGeneric(
    "writeFASTA",
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

