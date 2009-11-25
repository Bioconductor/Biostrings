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

setMethod("chartr", c(old="ANY", new="ANY", x="XString"),
    function(old, new, x) xvcopy(x, lkup=.mkOldToNewLkup(old, new, x))
)

setMethod("chartr", c(old="ANY", new="ANY", x="XStringSet"),
    function(old, new, x) xvcopy(x, lkup=.mkOldToNewLkup(old, new, x))
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

