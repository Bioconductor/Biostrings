### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Longest Common Prefix: the "lcprefix" new generic
###

### 's1' and 's2' must be BString (or derived) objects of the same class.
### Return the length of the Longest Common Prefix.
BString.lcprefix <- function(s1, s2)
{
    stop("coming soon...")
}

setGeneric(
    "lcprefix", signature=c("s1", "s2"),
    function(s1, s2) standardGeneric("lcprefix")
)
setMethod(
    "lcprefix", signature(s1="character", s2="character"),
    function(s1, s2)
        BString.lcprefix(BString(s1), BString(s2))
)
setMethod(
    "lcprefix", signature(s1="character", s2="BString"),
    function(s1, s2)
        BString.lcprefix(new(class(s2), s1), s2)
)
setMethod(
    "lcprefix", signature(s1="BString", s2="character"),
    function(s1, s2)
        BString.lcprefix(s1, new(class(s1), s2))
)
setMethod(
    "lcprefix", signature(s1="BString", s2="BString"),
    function(s1, s2)
    {
        if (class(s1) != class(s2))
            stop("'s1' and 's2' are not of the same class")
        BString.lcprefix(s1, s2)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Longest Common Suffix: the "lcsuffix" new generic
###

### 's1' and 's2' must be BString (or derived) objects of the same class.
### Return the length of the Longest Common Suffix.
BString.lcsuffix <- function(s1, s2)
{
    stop("coming soon...")
}

setGeneric(
    "lcsuffix", signature=c("s1", "s2"),
    function(s1, s2) standardGeneric("lcsuffix")
)
setMethod(
    "lcsuffix", signature(s1="character", s2="character"),
    function(s1, s2)
        BString.lcsuffix(BString(s1), BString(s2))
)
setMethod(
    "lcsuffix", signature(s1="character", s2="BString"),
    function(s1, s2)
        BString.lcsuffix(new(class(s2), s1), s2)
)
setMethod(
    "lcsuffix", signature(s1="BString", s2="character"),
    function(s1, s2)
        BString.lcsuffix(s1, new(class(s1), s2))
)
setMethod(
    "lcsuffix", signature(s1="BString", s2="BString"),
    function(s1, s2)
    {
        if (class(s1) != class(s2))
            stop("'s1' and 's2' are not of the same class")
        BString.lcsuffix(s1, s2)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Longest Common Substring: the "lcsubstr" new generic
###

### 's1' and 's2' must be BString (or derived) objects of the same class.
### Return a BStringPartialMatches object.
BString.lcsubstr <- function(s1, s2)
{
    stop("coming soon...")
}

setGeneric(
    "lcsubstr", signature=c("s1", "s2"),
    function(s1, s2) standardGeneric("lcsubstr")
)
setMethod(
    "lcsubstr", signature(s1="character", s2="character"),
    function(s1, s2)
        BString.lcsubstr(BString(s1), BString(s2))
)
setMethod(
    "lcsubstr", signature(s1="character", s2="BString"),
    function(s1, s2)
        BString.lcsubstr(new(class(s2), s1), s2)
)
setMethod(
    "lcsubstr", signature(s1="BString", s2="character"),
    function(s1, s2)
        BString.lcsubstr(s1, new(class(s1), s2))
)
setMethod(
    "lcsubstr", signature(s1="BString", s2="BString"),
    function(s1, s2)
    {
        if (class(s1) != class(s2))
            stop("'s1' and 's2' are not of the same class")
        BString.lcsubstr(s1, s2)
    }
)

