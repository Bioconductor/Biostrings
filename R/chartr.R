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
    xdata <- XRaw(lx)
    if (reverse) {
        XRaw.reverseCopy(xdata, x@offset + 1, x@offset + lx, src=x@xdata, lkup=lkup)
    } else {
        XRaw.copy(xdata, x@offset + 1, x@offset + lx, src=x@xdata, lkup=lkup)
    }
    new(class(x), xdata=xdata, length=length(xdata))
}

XStringSet.tr <- function(x, lkup=NULL, reverse=FALSE, use.names=TRUE)
{
    frame <- reduce(x, with.inframe.attrib=TRUE)
    super <- .Call("XStringSet_char_translate",
                   frame, lkup, reverse,
                   PACKAGE="Biostrings")
    ranges <- attr(frame, "inframe")
    newXStringSet(class(x), super, ranges, use.names=use.names, names=names(x))
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

setMethod("chartr", "XStringViews",
    function(old, new, x)
    {
        x@subject <- chartr(old, new, subject(x))
        x
    }
)

setMethod("chartr", "MaskedXString",
    function(old, new, x)
    {
        if (any(active(masks(x))))
            stop("\"chartr\" method for MaskedXString objects ",
                 "with active masks not ready yet\nPlease complain!")
        ans <- chartr(old, new, unmasked(x))
        masks(ans) <- masks(x)
        ans
    }
)

