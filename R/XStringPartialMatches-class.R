### =========================================================================
### XStringPartialMatches objects
### -------------------------------------------------------------------------
### A XStringPartialMatches object contains a set of partial matches
### on the same XString object, the subject string.

setClass("XStringPartialMatches",
    contains="XStringViews",
    representation(
        subpatterns="XStringViews"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods

setGeneric("subpatterns",
    function(x)
    {
        .Defunct()
        standardGeneric("subpatterns")
    }
)

setMethod("subpatterns", "XStringPartialMatches", function(x) x@subpatterns)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method
###

setMethod("show", "XStringPartialMatches",
    function(object)
    {
        msg <- "XStringPartialMatches objects are defunct"
        .Defunct(msg=wmsg(msg))
        subject <- subject(object)
        lsub <- length(subject)
        cat("  Views on a ", lsub, "-letter ",
            class(subject), " subject", sep="")
        #if (!is.null(subject@codec))
        #    cat(" with alphabet:", toString(subject@codec@letters))
        cat("\nSubject:", toSeqSnippet(subject, 70))
        XStringViews.show_vframe(object)

        pattern <- pattern(object)
        lpat <- length(pattern)
        cat("  Views on a ", lpat, "-letter ",
            class(pattern), " pattern", sep="")
        #if (!is.null(pattern@codec))
        #    cat(" with alphabet:", toString(pattern@codec@letters))
        cat("\nPattern:", toSeqSnippet(pattern, 70))
        XStringViews.show_vframe(subpatterns(object))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting
###

setMethod("[", "XStringPartialMatches",
    function(x, i, j, ..., drop)
    {
        msg <- "XStringPartialMatches objects are defunct"
        .Defunct(msg=wmsg(msg))
        ans <- callNextMethod()
        ans@subpatterns <- ans@subpatterns[i]
        ans
    }
)

