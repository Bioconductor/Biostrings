### =========================================================================
### The findPalindromes() generic & related functions
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "findPalindromes" and "findComplementedPalindromes" generics and
### methods.
###

### Return a list with the "start" and the "end" components.
.find.palindromes <- function(subject, min.armlength,
                              max.looplength, min.looplength,
                              max.mismatch, L2R_lkup)
{
    ## check min.armlength
    if (!isSingleNumber(min.armlength))
        stop("'min.armlength' must be a single integer")
    min.armlength <- as.integer(min.armlength)
    if (min.armlength < 2)
        stop("'min.armlength' must be >= 2")
    ## check max.looplength
    if (!isSingleNumber(max.looplength))
        stop("'max.looplength' must be a single integer")
    max.looplength <- as.integer(max.looplength)
    if (max.looplength < 0)
        stop("'max.looplength' must be a non-negative integer")
    ## check min.looplength
    if (!isSingleNumber(min.looplength))
        stop("'min.looplength' must be a single integer")
    min.looplength <- as.integer(min.looplength)
    if (min.looplength > max.looplength)
        stop("'min.looplength' must be <= 'max.looplength'")
    if (min.looplength < 0)
        stop("'min.looplength' must be a non-negative integer")
    if (min.looplength >= 1)
        stop("'min.looplength' >= 1 not yet supported (will be very soon)")
    ## check max.mismatch
    max.mismatch <- normargMaxMismatch(max.mismatch)
    if (max.mismatch != 0)
        stop("'max.mismatch' != 0 not yet supported (will be very soon)")
    C_ans <- .Call("find_palindromes",
                   subject@xdata@xp, subject@offset, subject@length,
                   min.armlength, max.looplength, L2R_lkup,
                   PACKAGE="Biostrings")
    ans_start <- C_ans$start
    ans_width <- C_ans$end - ans_start + 1L
    new("XStringViews", subject,
        start=ans_start, width=ans_width, check=FALSE)
}

setGeneric("findPalindromes", signature="subject",
    function(subject, min.armlength=4,
                      max.looplength=1, min.looplength=0,
                      max.mismatch=0)
        standardGeneric("findPalindromes")
)
setGeneric("findComplementedPalindromes", signature="subject",
    function(subject, min.armlength=4,
                      max.looplength=1, min.looplength=0,
                      max.mismatch=0)
        standardGeneric("findComplementedPalindromes")
)

setMethod("findPalindromes", "XString",
    function(subject, min.armlength=4,
                      max.looplength=1, min.looplength=0,
                      max.mismatch=0)
    {
        .find.palindromes(subject, min.armlength,
                          max.looplength, min.looplength,
                          max.mismatch, NULL)
    }
)
setMethod("findComplementedPalindromes", "DNAString",
    function(subject, min.armlength=4,
                      max.looplength=1, min.looplength=0,
                      max.mismatch=0)
    {
        .find.palindromes(subject, min.armlength,
                          max.looplength, min.looplength,
                          max.mismatch, getDNAComplementLookup())
    }
)

### WARNING: Unlike with the "findPalindromes" method for XString objects, the
### XStringViews object returned by this method is not guaranteed to have its
### views ordered from left to right! One important particular case where this
### is guaranteed though is when 'isNormal(subject)' is TRUE (i.e. 'subject' is
### a normal XStringViews object).
setMethod("findPalindromes", "XStringViews",
    function(subject, min.armlength=4,
                      max.looplength=1, min.looplength=0,
                      max.mismatch=0)
    {
        ans_start <- ans_width <- integer(0)
        for (i in seq_len(length(subject))) {
            pals <- findPalindromes(subject[[i]],
                                    min.armlength=min.armlength,
                                    max.looplength=max.looplength,
                                    min.looplength=min.looplength,
                                    max.mismatch=max.mismatch)
            offset <- start(subject)[i] - 1L
            ans_start <- c(ans_start, offset + start(pals))
            ans_width <- c(ans_width, width(pals))
        }
        new("XStringViews", subject(subject),
            start=ans_start, width=ans_width, check=FALSE)
    }
)
setMethod("findPalindromes", "MaskedXString",
    function(subject, min.armlength=4,
                      max.looplength=1, min.looplength=0,
                      max.mismatch=0)
    {
        findPalindromes(toXStringViewsOrXString(subject),
                        min.armlength=min.armlength,
                        max.looplength=max.looplength,
                        min.looplength=min.looplength,
                        max.mismatch=max.mismatch)
    }
)
setMethod("findComplementedPalindromes", "XStringViews",
    function(subject, min.armlength=4,
                      max.looplength=1, min.looplength=0,
                      max.mismatch=0)
    {
        ans_start <- ans_width <- integer(0)
        for (i in seq_len(length(subject))) {
            pals <- findComplementedPalindromes(subject[[i]],
                                                min.armlength=min.armlength,
                                                max.looplength=max.looplength,
                                                min.looplength=min.looplength,
                                                max.mismatch=max.mismatch)
            offset <- start(subject)[i] - 1L
            ans_start <- c(ans_start, offset + start(pals))
            ans_width <- c(ans_width, width(pals))
        }
        new("XStringViews", subject(subject),
            start=ans_start, width=ans_width, check=FALSE)
    }
)
setMethod("findComplementedPalindromes", "MaskedXString",
    function(subject, min.armlength=4,
                      max.looplength=1, min.looplength=0,
                      max.mismatch=0)
    {
        findComplementedPalindromes(toXStringViewsOrXString(subject),
                                    min.armlength=min.armlength,
                                    max.looplength=max.looplength,
                                    min.looplength=min.looplength,
                                    max.mismatch=max.mismatch)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "palindromeArmLength" and "complementedPalindromeArmLength" generics
### and methods.
###

setGeneric("palindromeArmLength", signature="x",
    function(x, max.mismatch=0, ...)
        standardGeneric("palindromeArmLength")
)
setGeneric("complementedPalindromeArmLength", signature="x",
    function(x, max.mismatch=0, ...)
        standardGeneric("complementedPalindromeArmLength")
)

setMethod("palindromeArmLength", "XString",
    function(x, max.mismatch=0, ...)
    {
        max.mismatch <- normargMaxMismatch(max.mismatch)
        revx <- reverse(x)
        armlength <- lcprefix(x, revx)
        if (armlength == 0L)
            stop("'x' is not a palindrome (no arms found)")
        armlength
    }
)
setMethod("complementedPalindromeArmLength", "DNAString",
    function(x, max.mismatch=0, ...)
    {
        max.mismatch <- normargMaxMismatch(max.mismatch)
        revx <- reverseComplement(x)
        armlength <- lcprefix(x, revx)
        if (armlength == 0L)
            stop("'x' is not a complemented palindrome (no arms found)")
        armlength
    }
)

setMethod("palindromeArmLength", "XStringViews",
    function(x, max.mismatch=0, ...)
    {
        if (length(x) == 0)
            return(integer(0))
        sapply(seq_len(length(x)),
               function(i)
                   palindromeArmLength(x[[i]],
                       max.mismatch=max.mismatch, ...))
    }
)
setMethod("complementedPalindromeArmLength", "XStringViews",
    function(x, max.mismatch=0, ...)
    {
        if (length(x) == 0)
            return(integer(0))
        sapply(seq_len(length(x)),
               function(i)
                   complementedPalindromeArmLength(x[[i]],
                       max.mismatch=max.mismatch, ...))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "palindromeLeftArm" and "complementedPalindromeLeftArm" generics and
### methods.
###

setGeneric("palindromeLeftArm", signature="x",
    function(x, max.mismatch=0, ...)
        standardGeneric("palindromeLeftArm")
)
setGeneric("complementedPalindromeLeftArm", signature="x",
    function(x, max.mismatch=0, ...)
        standardGeneric("complementedPalindromeLeftArm")
)

setMethod("palindromeLeftArm", "XString",
    function(x, max.mismatch=0, ...)
    {
        XString.substr(x, 1L,
            palindromeArmLength(x, max.mismatch=max.mismatch, ...))
    }
)
setMethod("complementedPalindromeLeftArm", "DNAString",
    function(x, max.mismatch=0, ...)
    {
        XString.substr(x, 1L,
            complementedPalindromeArmLength(x, max.mismatch=max.mismatch, ...))
    }
)

setMethod("palindromeLeftArm", "XStringViews",
    function(x, max.mismatch=0, ...)
    {
        ans_start <- start(x)
        ans_width <- palindromeArmLength(x, max.mismatch=max.mismatch, ...)
        new("XStringViews", subject(x),
            start=ans_start, width=ans_width, check=FALSE)
    }
)
setMethod("complementedPalindromeLeftArm", "XStringViews",
    function(x, max.mismatch=0, ...)
    {
        ans_start <- start(x)
        ans_width <- complementedPalindromeArmLength(x,
                         max.mismatch=max.mismatch, ...)
        new("XStringViews", subject(x),
            start=ans_start, width=ans_width, check=FALSE)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "palindromeRightArm" and "complementedPalindromeRightArm" generics
### and methods.
###

setGeneric("palindromeRightArm", signature="x",
    function(x, max.mismatch=0, ...)
        standardGeneric("palindromeRightArm")
)
setGeneric("complementedPalindromeRightArm", signature="x",
    function(x, max.mismatch=0, ...)
        standardGeneric("complementedPalindromeRightArm")
)

setMethod("palindromeRightArm", "XString",
    function(x, max.mismatch=0, ...)
    {
        start <- nchar(x) - palindromeArmLength(x,
                                max.mismatch=max.mismatch, ...) + 1L
        XString.substr(x, start, nchar(x))
    }
)
setMethod("complementedPalindromeRightArm", "DNAString",
    function(x, max.mismatch=0, ...)
    {
        start <- nchar(x) - complementedPalindromeArmLength(x,
                                max.mismatch=max.mismatch, ...) + 1L
        XString.substr(x, start, nchar(x))
    }
)

setMethod("palindromeRightArm", "XStringViews",
    function(x, max.mismatch=0, ...)
    {
        ans_width <- palindromeArmLength(x,
                         max.mismatch=max.mismatch, ...)
        ans_start <- end(x) - ans_width + 1L
        new("XStringViews", subject(x),
            start=ans_start, width=ans_width, check=FALSE)
    }
)
setMethod("complementedPalindromeRightArm", "XStringViews",
    function(x, max.mismatch=0, ...)
    {
        ans_width <- complementedPalindromeArmLength(x,
                         max.mismatch=max.mismatch, ...)
        ans_start <- end(x) - ans_width + 1L
        new("XStringViews", subject(x),
            start=ans_start, width=ans_width, check=FALSE)
    }
)

