### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "chartr" generic function and methods.
###

BString.tr <- function(x, lkup=NULL, reverse=FALSE)
{
    lx <- length(x)
    data <- XRaw(lx)
    if (reverse) {
        XRaw.reverseCopy(data, x@offset + 1, x@offset + lx, src=x@data, lkup=lkup)
    } else {
        XRaw.copy(data, x@offset + 1, x@offset + lx, src=x@data, lkup=lkup)
    }
    new(class(x), data, 0L, length(data), check=FALSE)
}

### This setGeneric() statement will unfortunately cause the following message
### at installation time:
###   New generic for "chartr" does not agree with implicit generic from
###   package "base"; a new generic will be assigned with package "Biostrings"
### But we need this setGeneric() statement anyway otherwise the dispatch would
### happen on the 'old' argument and not on the 'x' argument.
setGeneric("chartr", signature="x",
    function(old, new, x) standardGeneric("chartr")
)

setMethod("chartr", "BString",
    function(old, new, x)
    {
        if (class(old) != class(x))
            old <- mkBString(class(x), old)
        if (class(new) != class(x))
            new <- mkBString(class(x), new)
        if (nchar(old) != nchar(new))
            stop("'old' and 'new' must have the same length")
        old_codes <- BString.readInts(old, 1, nchar(old))
        new_codes <- BString.readInts(new, 1, nchar(new))
        lkup <- buildLookupTable(codes(x), codes(x))
        lkup[1 + old_codes] <- new_codes
        BString.tr(x, lkup=lkup)
    }
)

setMethod("chartr", "BStringViews",
    function(old, new, x)
        error("not ready yet, sorry!")
)

setMethod("chartr", "BStringSet",
    function(old, new, x)
        error("not ready yet, sorry!")
)

