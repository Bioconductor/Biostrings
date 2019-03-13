### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "injectHardMask" generic function and methods.
###

setGeneric("injectHardMask", signature="x",
    function(x, letter="+") standardGeneric("injectHardMask")
)

setMethod("injectHardMask", "XStringViews",
    function(x, letter="+")
    {
        x_seqtype <- seqtype(x)
        if (is(letter, "XString") && seqtype(letter) == x_seqtype) {
            if (length(letter) != 1)
                stop("'letter' must be a single letter")
        } else {
            if (!isSingleString(letter) || nchar(letter) != 1)
                stop("'letter' must be a single letter")
            letter <- XString(x_seqtype, letter)
        }
        code <- XString.readCodes(letter, 1L)
        y <- gaps(x)
        ## Because y is obtained with the "gaps" method for XStringViews
        ## objects, then the set of ranges defined by start(y) and width(y)
        ## is normal and is guaranteed to be within the limits of y.
        ## Hence start(y) and width(y) can be considered safe.
        .Call2("XString_inject_code",
              subject(y), start(y), width(y), code,
              PACKAGE="Biostrings")
    }
)

setMethod("injectHardMask", "MaskedXString",
    function(x, letter="+")
    {
        y <- as(x, "XStringViews")
        injectHardMask(y, letter=letter)
    }
)

