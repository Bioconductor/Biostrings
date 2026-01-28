### =========================================================================
### Find Longest Common Prefix/Suffix
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Longest Common Prefix
###

### 's1' and 's2' must be XString objects containing sequences of the same
### type. Return the length (integer) of the Longest Common Prefix.
XString.lcprefix <- function(s1, s2)
{
    .Call2("lcprefix", s1@shared@xp, s1@offset, s1@length,
                      s2@shared@xp, s2@offset, s2@length,
                      PACKAGE="Biostrings")
}

setGeneric("lcprefix", signature=c("s1", "s2"),
    function(s1, s2) standardGeneric("lcprefix")
)
setMethod("lcprefix", signature(s1="character", s2="character"),
    function(s1, s2)
        XString.lcprefix(BString(s1), BString(s2))
)
setMethod("lcprefix", signature(s1="character", s2="XString"),
    function(s1, s2)
        XString.lcprefix(XString(seqtype(s2), s1), s2)
)
setMethod("lcprefix", signature(s1="XString", s2="character"),
    function(s1, s2)
        XString.lcprefix(s1, XString(seqtype(s1), s2))
)
setMethod("lcprefix", signature(s1="XString", s2="XString"),
    function(s1, s2)
    {
        if (seqtype(s1) != seqtype(s2))
            stop("'s1' and 's2' must be XString objects containing ",
                 "sequences of the same type")
        XString.lcprefix(s1, s2)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Longest Common Suffix
###

### 's1' and 's2' must be XString objects containing sequences of the same
### type. Return the length (integer) of the Longest Common Suffix.
XString.lcsuffix <- function(s1, s2)
{
    .Call2("lcsuffix", s1@shared@xp, s1@offset, s1@length,
                      s2@shared@xp, s2@offset, s2@length,
                      PACKAGE="Biostrings")
}

setGeneric("lcsuffix", signature=c("s1", "s2"),
    function(s1, s2) standardGeneric("lcsuffix")
)
setMethod("lcsuffix", signature(s1="character", s2="character"),
    function(s1, s2)
        XString.lcsuffix(BString(s1), BString(s2))
)
setMethod("lcsuffix", signature(s1="character", s2="XString"),
    function(s1, s2)
        XString.lcsuffix(XString(seqtype(s2), s1), s2)
)
setMethod("lcsuffix", signature(s1="XString", s2="character"),
    function(s1, s2)
        XString.lcsuffix(s1, XString(seqtype(s1), s2))
)
setMethod("lcsuffix", signature(s1="XString", s2="XString"),
    function(s1, s2)
    {
        if (seqtype(s1) != seqtype(s2))
            stop("'s1' and 's2' must be XString objects containing ",
                 "sequences of the same type")
        XString.lcsuffix(s1, s2)
    }
)

