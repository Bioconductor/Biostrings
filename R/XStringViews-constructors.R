### =========================================================================
### Constructor-like functions and generics for XStringViews objects
### -------------------------------------------------------------------------

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "XStringViews" generic and methods.
###
### 'subjectClass' must be the name of one of the direct XString subclasses
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
        ## drop the "String" suffix
        subject_basetype <- substr(subjectClass, 1, nchar(subjectClass)-6)
        subject <- XString(subject_basetype, seq)
        ans <- successiveViews(subject, nchar(x), gapwidth=nchar(collapse))
        names(ans) <- names(x)
        ans
    }
)

setMethod("XStringViews", "XString",
    function(x, subjectClass, collapse="")
    {
        if (!missing(collapse)) {
            ## The semantic is: views are delimited by the occurrences of
            ## 'collapse' in 'x' (a kind of strsplit() for XString objects).
            ## Uncomment when normalize() and ! method are ready (see TODO file):
            #return(!normalize(matchPattern(collapse, b, fixed=TRUE)))
            stop("'collapse' not yet supported when 'x' is an XString object")
        }
        if (!missing(subjectClass) && subjectClass != class(x)) {
            ## drop the "String" suffix
            subject_basetype <- substr(subjectClass, 1, nchar(subjectClass)-6)
            x <- XString(subject_basetype, x)
        }
        unsafe.newXStringViews(x, 1L, nchar(x))
    }
)

### The 'collapse' arg is ignored.
setMethod("XStringViews", "XStringViews",
    function(x, subjectClass, collapse="")
    {
        if (!missing(collapse))
            stop("'collapse' not supported when 'x' is an XStringViews object")
        ## drop the "String" suffix
        subject_basetype <- substr(subjectClass, 1, nchar(subjectClass)-6)
        x@subject <- XString(subject_basetype, subject(x))
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Old stuff (Deprecated or Defunct).
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

adjacentViews <- function(subject, width, gapwidth=0)
{
    .Deprecated("successiveViews")
    successiveViews(subject, width, gapwidth=gapwidth)
}

