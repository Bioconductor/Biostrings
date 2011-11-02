### =========================================================================
### The XStringViews constructor
### -------------------------------------------------------------------------
###
### Stuff in this file is defunct.
###
### 'subjectClass' must be the name of one of the direct XString subclasses
### i.e. "BString", "DNAString", "RNAString" or "AAString".
###

setGeneric("XStringViews", signature="x",
    function(x, subjectClass, collapse="") standardGeneric("XStringViews")
)

setMethod("XStringViews", "ANY",
    function(x, subjectClass, collapse="")
    {
        if (!isSingleString(subjectClass))
            stop("'subjectClass' must be a single string")
        msg <- c("  Using XStringViews() on a character vector is ",
                 "defunct.\n  Please use instead something like:\n",
                 "      as(", subjectClass, "Set(x)), \"Views\")\n",
                 "  if you really want views, otherwise just:\n",
                 "      ", subjectClass, "Set(x)")
        .Defunct(msg=msg)
    }
)

setMethod("XStringViews", "XString",
    function(x, subjectClass, collapse="")
    {
        msg <- c("  Using XStringViews() on a ", class(x), " object is ",
                 "defunct.\n  Please use 'as(x, \"Views\")' instead.")
        .Defunct(msg=msg)
    }
)

setMethod("XStringViews", "XStringViews",
    function(x, subjectClass, collapse="")
    {
        if (!isSingleString(subjectClass))
            stop("'subjectClass' must be a single string")
        if (!missing(collapse))
            stop("'collapse' is not supported ",
                 "when 'x' is an XStringViews object")
        ## drop the "String" suffix
        basetype <- substr(subjectClass, 1, nchar(subjectClass)-6)
        msg <-  c("  Using XStringViews(..., subjectClass=\"",
                  subjectClass, "\") on an XStringViews\n  object is ",
                  "defunct. Please use 'xsbasetype(x) <- \"",
                  basetype, "\"' instead.")
        .Defunct(msg=msg)
    }
)

