### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "injectHardMask" generic function and methods.
###

setGeneric("injectHardMask", signature="x",
    function(x, letter="+") standardGeneric("injectHardMask")
)

setMethod("injectHardMask", "XStringViews",
    function(x, letter="+")
    {
        base_class <- baseXStringSubtype(x)
        if (class(letter) == base_class) {
            if (length(letter) != 1)
                stop("'letter' must be a single letter")
        } else {
            if (!isSingleString(letter) || nchar(letter) != 1)
                stop("'letter' must be a single letter")
            letter <- XString(base_class, letter)
        }
        code <- XString.readCodes(letter, 1L)
        y <- gaps(x)
        ## Because y is obtained with the "gaps" method for XStringViews
        ## objects, then the set of ranges defined by start(y) and width(y)
        ## is normal and is guaranteed to be within the limits of y.
        ## Hence start(y) and width(y) can be considered safe.
        .Call("inject_code",
              subject(y), start(y), width(y), code,
              PACKAGE="Biostrings")
    }
)

setMethod("injectHardMask", "MaskedDNAString",
    function(x, letter="+")
    {
        y <- as(x, "XStringViews")
        injectHardMask(y, letter=letter)
    }
)

