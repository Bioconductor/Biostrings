### =========================================================================
### The BStringPartialMatches class
### -------------------------------------------------------------------------
### A BStringPartialMatches object contains a set of partial matches
### on the same BString object, the subject string.

setClass(
    "BStringPartialMatches",
    contains="BStringViews",
    representation(
        subpatterns="BStringViews"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods

setGeneric("subpatterns", function(x) standardGeneric("subpatterns"))
setMethod("subpatterns", "BStringPartialMatches", function(x) x@subpatterns)

setGeneric("pattern", function(x) standardGeneric("pattern"))
setMethod("pattern", "BStringPartialMatches", function(x) x@subpatterns@subject)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method
###

setMethod("show", "BStringPartialMatches",
    function(object)
    {
        subject <- subject(object)
        lsub <- length(subject)
        cat("  Views on a ", lsub, "-letter ",
            class(subject), " subject", sep="")
        #if (!is.null(subject@codec))
        #    cat(" with alphabet:", toString(subject@codec@letters))
        cat("\nSubject:", BString.get_snippet(subject, 70))
        BStringViews.show_vframe(object, 5L)

        pattern <- pattern(object)
        lpat <- length(pattern)
        cat("  Views on a ", lpat, "-letter ",
            class(pattern), " pattern", sep="")
        #if (!is.null(pattern@codec))
        #    cat(" with alphabet:", toString(pattern@codec@letters))
        cat("\nPattern:", BString.get_snippet(pattern, 70))
        BStringViews.show_vframe(subpatterns(object), 5L)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting
###

setMethod("[", "BStringPartialMatches",
    function(x, i, j, ..., drop)
    {
        ans <- callNextMethod()
        ans@subpatterns <- ans@subpatterns[i]
        ans
    }
)

