### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "chartr" generic function and methods.
###

.mkOldToNewLkup <- function(old, new, x)
{
    x_seqtype <- seqtype(x)
    if (!is(old, "XString") || seqtype(old) != x_seqtype)
        old <- XString(x_seqtype, old)
    if (!is(new, "XString") || seqtype(new) != x_seqtype)
        new <- XString(x_seqtype, new)
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

