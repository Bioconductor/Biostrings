### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "chartr" generic function and methods.
###

.mkOldToNewLkup <- function(old, new, x)
{
    basetype <- xsbasetype(x)
    if (!is(old, "XString") || xsbasetype(old) != basetype)
        old <- XString(basetype, old)
    if (!is(new, "XString") || xsbasetype(new) != basetype)
        new <- XString(basetype, new)
    if (nchar(old) != nchar(new))
        stop("'old' and 'new' must have the same length")
    old_codes <- XString.readCodes(old, 1, nchar(old))
    new_codes <- XString.readCodes(new, 1, nchar(new))
    lkup <- buildLookupTable(xscodes(x), xscodes(x))
    lkup[1 + old_codes] <- new_codes
    lkup
}

XString.tr <- function(x, lkup=NULL, reverse=FALSE)
{
    lx <- length(x)
    xdata <- RawPtr(lx)
    if (reverse) {
        RawPtr.reverseCopy(xdata, x@offset + 1, x@offset + lx, src=x@xdata, lkup=lkup)
    } else {
        RawPtr.copy(xdata, x@offset + 1, x@offset + lx, src=x@xdata, lkup=lkup)
    }
    new(class(x), xdata=xdata, length=length(xdata))
}

XStringSet.tr <- function(x, lkup=NULL, reverse=FALSE, use.names=TRUE)
{
    x@ranges <- reduce(x@ranges, with.inframe.attrib=TRUE)
    ans_super <- .Call("XStringSet_char_translate",
                       x, lkup, reverse,
                       PACKAGE="Biostrings")
    ans_ranges <- attr(x@ranges, "inframe")
    unsafe.newXStringSet(ans_super, ans_ranges, use.names=use.names, names=names(x))
}

setMethod("chartr", c(old = "ANY", new = "ANY", x = "XString"),
    function(old, new, x)
    {
        lkup <- .mkOldToNewLkup(old, new, x)
        XString.tr(x, lkup=lkup)
    }
)

setMethod("chartr", c(old = "ANY", new = "ANY", x = "XStringSet"),
    function(old, new, x)
    {
        lkup <- .mkOldToNewLkup(old, new, x)
        XStringSet.tr(x, lkup=lkup, use.names=TRUE)
    }
)

setMethod("chartr", c(old = "ANY", new = "ANY", x = "XStringViews"),
    function(old, new, x)
    {
        x@subject <- chartr(old, new, subject(x))
        x
    }
)

setMethod("chartr", c(old = "ANY", new = "ANY", x = "MaskedXString"),
    function(old, new, x)
    {
        if (any(active(masks(x))))
            stop("\"chartr\" method for MaskedXString objects ",
                 "with active masks not ready yet\n  Please complain!")
        ans <- chartr(old, new, unmasked(x))
        masks(ans) <- masks(x)
        ans
    }
)
