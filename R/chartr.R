### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "chartr" generic function and methods.
###

.mkOldToNewLkup <- function(old, new, x)
{
    if (class(old) != class(x))
        old <- XString(class(x), old)
    if (class(new) != class(x))
        new <- XString(class(x), new)
    if (nchar(old) != nchar(new))
        stop("'old' and 'new' must have the same length")
    old_codes <- XString.readCodes(old, 1, nchar(old))
    new_codes <- XString.readCodes(new, 1, nchar(new))
    lkup <- buildLookupTable(codes(x), codes(x))
    lkup[1 + old_codes] <- new_codes
    lkup
}

XString.tr <- function(x, lkup=NULL, reverse=FALSE)
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

XStringSet.tr <- function(x, lkup=NULL, reverse=FALSE, use.names=TRUE)
{
    frame <- reduce(x, with.inframe.attrib=TRUE)
    super <- .Call("XStringSet_char_translate",
                   frame, lkup, reverse,
                   PACKAGE="Biostrings")
    ranges <- attr(frame, "inframe")
    newXStringSet(class(x), super, ranges, x, use.names)
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

setMethod("chartr", "XString",
    function(old, new, x)
    {
        lkup <- .mkOldToNewLkup(old, new, x)
        XString.tr(x, lkup=lkup)
    }
)

setMethod("chartr", "XStringSet",
    function(old, new, x)
    {
        lkup <- .mkOldToNewLkup(old, new, super(x))
        XStringSet.tr(x, lkup=lkup, use.names=TRUE)
    }
)

setMethod("chartr", "BStringViews",
    function(old, new, x)
    {
        x@subject <- chartr(old, new, subject(x))
        x
    }
)

