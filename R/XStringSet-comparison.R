### =========================================================================
### Comparing and ordering the elements in one or more XStringSet objects
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### compare().
###

setMethod("compare", c("XStringSet", "XStringSet"),
    function(x, y)
    {
        if (!comparable_seqtypes(seqtype(x), seqtype(y)))
            stop("comparison between a \"", class(x), "\" instance ",
                 "and a \"", class(y), "\" instance ",
                 "is not supported")
        callNextMethod()  # call method for XRawList objects
    }
)

setMethod("compare", c("XStringSet", "character"),
    function(x, y)
    {
        y <- try(XStringSet(seqtype(x), y))
        if (is(y, "try-error"))
            stop("could not turn 'y' into a ", seqtype(x), " instance")
        callGeneric()  # call method for XStringSet,XStringSet
    }
)

setMethod("compare", c("character", "XStringSet"),
    function(x, y)
    {
        x <- try(XStringSet(seqtype(y), x))
        if (is(x, "try-error"))
            stop("could not turn 'x' into a ", seqtype(y), " instance")
        callGeneric()  # call method for XStringSet,XStringSet
    }
)

setMethod("compare", c("XStringSet", "XString"),
    function(x, y)
    {
        y <- as(y, "XStringSet")
        callGeneric()  # call method for XStringSet,XStringSet
    }
)

setMethod("compare", c("XString", "XStringSet"),
    function(x, y)
    {
        x <- as(x, "XStringSet")
        callGeneric()  # call method for XStringSet,XStringSet
    }
)

setMethod("compare", c("character", "XString"),
    function(x, y)
    {
        y <- as(y, "XStringSet")
        callGeneric()  # call method for character,XStringSet
    }
)

setMethod("compare", c("XString", "character"),
    function(x, y)
    {
        x <- as(x, "XStringSet")
        callGeneric()  # call method for XStringSet,character
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### match().
###

setMethod("match", c("XStringSet", "XStringSet"),
    function(x, table, nomatch=NA_integer_, incomparables=NULL)
    {
        if (!comparable_seqtypes(seqtype(x), seqtype(table)))
            stop("match() between a \"", class(x), "\" instance ",
                 "and a \"", class(table), "\" instance ",
                 "is not supported")
        callNextMethod()  # call method for XRawList objects
    }
)

setMethod("match", c("XStringSet", "character"),
    function(x, table, nomatch=NA_integer_, incomparables=NULL)
    {
        table <- try(XStringSet(seqtype(x), table))
        if (is(table, "try-error"))
            stop("could not turn 'table' into a ", seqtype(x), " instance")
        callGeneric()  # call method for XStringSet,XStringSet
    }
)

setMethod("match", c("character", "XStringSet"),
    function(x, table, nomatch=NA_integer_, incomparables=NULL)
    {
        x <- try(XStringSet(seqtype(table), x))
        if (is(x, "try-error"))
            stop("could not turn 'x' into a ", seqtype(table), " instance")
        callGeneric()  # call method for XStringSet,XStringSet
    }
)

setMethod("match", c("XStringSet", "XString"),
    function(x, table, nomatch=NA_integer_, incomparables=NULL)
    {
        table <- as(table, "XStringSet")
        callGeneric()  # call method for XStringSet,XStringSet
    }
)

setMethod("match", c("XString", "XStringSet"),
    function(x, table, nomatch=NA_integer_, incomparables=NULL)
    {
        x <- as(x, "XStringSet")
        callGeneric()  # call method for XStringSet,XStringSet
    }
)

setMethod("match", c("character", "XString"),
    function(x, table, nomatch=NA_integer_, incomparables=NULL)
    {
        table <- as(table, "XStringSet")
        callGeneric()  # call method for character,XStringSet
    }
)

setMethod("match", c("XString", "character"),
    function(x, table, nomatch=NA_integer_, incomparables=NULL)
    {
        x <- as(x, "XStringSet")
        callGeneric()  # call method for XStringSet,character
    }
)

