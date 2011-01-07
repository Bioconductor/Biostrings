### =========================================================================
### The XStringViews constructor
### -------------------------------------------------------------------------
###
### TODO: This is old and ugly stuff which has very little value. It needs to
### be revisited just to confirm that, most of the times, the functionalities
### it provides can be achieved in better and more elegant ways. So, at some
### point, it should be deprecated.

### 'subjectClass' must be the name of one of the direct XString subclasses
### i.e. "BString", "DNAString", "RNAString" or "AAString".
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
        if (!missing(collapse))
            stop("'collapse' is not supported ",
                 "when 'x' is an XString object")
        if (!missing(subjectClass) && subjectClass != class(x)) {
            ## drop the "String" suffix
            subject_basetype <- substr(subjectClass, 1, nchar(subjectClass)-6)
            x <- XString(subject_basetype, x)
        }
        as(x, "Views")
    }
)

setMethod("XStringViews", "XStringViews",
    function(x, subjectClass, collapse="")
    {
        if (!missing(collapse))
            stop("'collapse' is not supported ",
                 "when 'x' is an XStringViews object")
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

setMethod("BStringViews", "ANY", .Defunct("XStringViews"))

setMethod("BStringViews", "file",
    function(src, subjectClass, collapse="")
        .Defunct("read.DNAStringSet")
)

adjacentViews <- function(subject, width, gapwidth=0)
    .Defunct("successiveViews", package="IRanges")

