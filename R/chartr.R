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

### TODO: Make this a "copy" method for XVector objects and move it to the
### IRanges package.
XString.tr <- function(x, lkup=NULL, reverse=FALSE)
{
    ans_length <- length(x)
    ans_shared <- copy(x@shared, start=x@offset+1L, width=ans_length,
                       lkup=lkup, reverse=reverse)
    new(class(x), shared=ans_shared, length=ans_length)
}

### FIXME: Currently broken if 'length(x@pool)) != 1'. This is because the
### "copy" method for SharedVector_Pool objects (defined in IRanges) is
### itself broken. Fix it!
### NOTE: Memory footprint could be reduced by copying only the regions in
### each x@pool element that are actually used by 'x'. See old code (commented
### out) for how this was done at the time of the old XStringSet container.
### TODO: Make this a "copy" method for XVectorList objects and move it to the
### IRanges package.
XStringSet.tr <- function(x, lkup=NULL, reverse=FALSE, use.names=TRUE)
{
    #x@ranges <- reduce(x@ranges, with.inframe.attrib=TRUE)
    #ans_super <- .Call("XStringSet_char_translate",
    #                   x, lkup, reverse,
    #                   PACKAGE="Biostrings")
    #ans_ranges <- attr(x@ranges, "inframe")
    #unsafe.newXStringSet(ans_super, ans_ranges, use.names=use.names, names=names(x))
    x@pool <- copy(x@pool, lkup=lkup)
    x
}

setMethod("chartr", c(old="ANY", new="ANY", x="XString"),
    function(old, new, x)
    {
        lkup <- .mkOldToNewLkup(old, new, x)
        XString.tr(x, lkup=lkup)
    }
)

setMethod("chartr", c(old="ANY", new="ANY", x="XStringSet"),
    function(old, new, x)
    {
        lkup <- .mkOldToNewLkup(old, new, x)
        XStringSet.tr(x, lkup=lkup, use.names=TRUE)
    }
)

setMethod("chartr", c(old="ANY", new="ANY", x="XStringViews"),
    function(old, new, x)
    {
        x@subject <- chartr(old, new, subject(x))
        x
    }
)

setMethod("chartr", c(old="ANY", new="ANY", x="MaskedXString"),
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

