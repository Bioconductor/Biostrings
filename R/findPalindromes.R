### =========================================================================
### The findPalindromes() generic & related functions
### -------------------------------------------------------------------------


.normalize.max.mismatch <- function(max.mismatch)
{
    if (!is.numeric(max.mismatch) || length(max.mismatch) != 1 || is.na(max.mismatch))
        stop("'max.mismatch' must be a single integer")
    max.mismatch <- as.integer(max.mismatch)
    if (max.mismatch < 0)
        stop("'max.mismatch' must be a non-negative integer")
    if (max.mismatch >= 1)
        stop("'max.mismatch' >= 1 not yet supported (will be very soon)")
    max.mismatch
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "findPalindromes" and "findComplementedPalindromes" generics and
### methods.
###

debug_find_palindromes <- function()
{
    invisible(.Call("find_palindromes_debug", PACKAGE="Biostrings"))
}

### Return a list with the "start" and the "end" components.
.find.palindromes <- function(subject, min.armlength,
                              max.looplength, min.looplength,
                              max.mismatch, L2R_lkup)
{
    ## check min.armlength
    if (!is.numeric(min.armlength) || length(min.armlength) != 1 || is.na(min.armlength))
        stop("'min.armlength' must be a single integer")
    min.armlength <- as.integer(min.armlength)
    if (min.armlength < 2)
        stop("'min.armlength' must be >= 2")
    ## check max.looplength
    if (!is.numeric(max.looplength) || length(max.looplength) != 1 || is.na(max.looplength))
        stop("'max.looplength' must be a single integer")
    max.looplength <- as.integer(max.looplength)
    if (max.looplength < 0)
        stop("'max.looplength' must be a non-negative integer")
    ## check min.looplength
    if (!is.numeric(min.looplength) || length(min.looplength) != 1 || is.na(min.looplength))
        stop("'min.looplength' must be a single integer")
    min.looplength <- as.integer(min.looplength)
    if (min.looplength > max.looplength)
        stop("'min.looplength' must be <= 'max.looplength'")
    if (min.looplength < 0)
        stop("'min.looplength' must be a non-negative integer")
    if (min.looplength >= 1)
        stop("'min.looplength' >= 1 not yet supported (will be very soon)")
    ## check max.mismatch
    max.mismatch <- .normalize.max.mismatch(max.mismatch)
    views <- .Call("find_palindromes",
                   subject@data@xp, subject@offset, subject@length,
                   min.armlength, max.looplength, L2R_lkup,
                   PACKAGE="Biostrings")
    new("BStringViews", subject=subject,
        start=views$start, end=views$end, check.views=FALSE)
}

setGeneric("findPalindromes", signature="subject",
    function(subject, min.armlength=4, max.looplength=1, min.looplength=0, max.mismatch=0)
        standardGeneric("findPalindromes")
)
setGeneric("findComplementedPalindromes", signature="subject",
    function(subject, min.armlength=4, max.looplength=1, min.looplength=0, max.mismatch=0)
        standardGeneric("findComplementedPalindromes")
)

setMethod("findPalindromes", "BString",
    function(subject, min.armlength=4, max.looplength=1, min.looplength=0, max.mismatch=0)
    {
        .find.palindromes(subject, min.armlength,
                          max.looplength, min.looplength,
                          max.mismatch, NULL)
    }
)
setMethod("findComplementedPalindromes", "DNAString",
    function(subject, min.armlength=4, max.looplength=1, min.looplength=0, max.mismatch=0)
    {
        .find.palindromes(subject, min.armlength,
                          max.looplength, min.looplength,
                          max.mismatch, getDNAComplementLookup())
    }
)

### WARNING: Unlike with the "findPalindromes" method for BString objects, the
### BStringViews object returned by this method is not guaranteed to have its
### views ordered from left to right! One important particular case where this
### is guaranteed though is when 'subject' is a normalized BStringViews object.
setMethod("findPalindromes", "BStringViews",
    function(subject, min.armlength=4, max.looplength=1, min.looplength=0, max.mismatch=0)
    {
        ans_start <- ans_end <- integer(0)
        for (i in seq_len(length(subject))) {
            pals <- findPalindromes(subject[[i]],
                                    min.armlength=min.armlength,
                                    max.looplength=max.looplength,
                                    min.looplength=min.looplength,
                                    max.mismatch=max.mismatch)
            offset <- start(subject)[i] - 1L
            ans_start <- c(ans_start, offset + start(pals))
            ans_end <- c(ans_end, offset + end(pals))
        }
        new("BStringViews", subject=subject(subject),
            start=ans_start, end=ans_end, check.views=FALSE)
    }
)
setMethod("findComplementedPalindromes", "BStringViews",
    function(subject, min.armlength=4, max.looplength=1, min.looplength=0, max.mismatch=0)
    {
        ans_start <- ans_end <- integer(0)
        for (i in seq_len(length(subject))) {
            pals <- findComplementedPalindromes(subject[[i]],
                                                min.armlength=min.armlength,
                                                max.looplength=max.looplength,
                                                min.looplength=min.looplength,
                                                max.mismatch=max.mismatch)
            offset <- start(subject)[i] - 1L
            ans_start <- c(ans_start, offset + start(pals))
            ans_end <- c(ans_end, offset + end(pals))
        }
        new("BStringViews", subject=subject(subject),
            start=ans_start, end=ans_end, check.views=FALSE)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "palindromeArmLength" and "complementedPalindromeArmLength" generics
### and methods.
###

setGeneric("palindromeArmLength", signature="x",
    function(x, max.mismatch=0, ...) standardGeneric("palindromeArmLength")
)
setGeneric("complementedPalindromeArmLength", signature="x",
    function(x, max.mismatch=0, ...) standardGeneric("complementedPalindromeArmLength")
)

setMethod("palindromeArmLength", "BString",
    function(x, max.mismatch=0, ...)
    {
        max.mismatch <- .normalize.max.mismatch(max.mismatch)
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
        max.mismatch <- .normalize.max.mismatch(max.mismatch)
        revx <- reverseComplement(x)
        armlength <- lcprefix(x, revx)
        if (armlength == 0L)
            stop("'x' is not a complemented palindrome (no arms found)")
        armlength
    }
)

setMethod("palindromeArmLength", "BStringViews",
    function(x, max.mismatch=0, ...)
    {
        if (length(x) == 0)
            return(integer(0))
        sapply(seq_len(length(x)), function(i)
               palindromeArmLength(x[[i]], max.mismatch=max.mismatch, ...))
    }
)
setMethod("complementedPalindromeArmLength", "BStringViews",
    function(x, max.mismatch=0, ...)
    {
        if (length(x) == 0)
            return(integer(0))
        sapply(seq_len(length(x)), function(i)
               complementedPalindromeArmLength(x[[i]], max.mismatch=max.mismatch, ...))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "palindromeLeftArm" and "complementedPalindromeLeftArm" generics and
### methods.
###

setGeneric("palindromeLeftArm", signature="x",
    function(x, max.mismatch=0, ...) standardGeneric("palindromeLeftArm")
)
setGeneric("complementedPalindromeLeftArm", signature="x",
    function(x, max.mismatch=0, ...) standardGeneric("complementedPalindromeLeftArm")
)

setMethod("palindromeLeftArm", "BString",
    function(x, max.mismatch=0, ...)
    {
        BString.substr(x, 1L, palindromeArmLength(x, max.mismatch=max.mismatch, ...))
    }
)
setMethod("complementedPalindromeLeftArm", "DNAString",
    function(x, max.mismatch=0, ...)
    {
        BString.substr(x, 1L, complementedPalindromeArmLength(x, max.mismatch=max.mismatch, ...))
    }
)

setMethod("palindromeLeftArm", "BStringViews",
    function(x, max.mismatch=0, ...)
    {
        arm_start <- start(x)
        arm_end <- arm_start + palindromeArmLength(x, max.mismatch=max.mismatch, ...) - 1L
        new("BStringViews", subject=subject(x),
            start=arm_start, end=arm_end, check.views=FALSE)
    }
)
setMethod("complementedPalindromeLeftArm", "BStringViews",
    function(x, max.mismatch=0, ...)
    {
        arm_start <- start(x)
        arm_end <- arm_start + complementedPalindromeArmLength(x, max.mismatch=max.mismatch, ...) - 1L
        new("BStringViews", subject=subject(x),
            start=arm_start, end=arm_end, check.views=FALSE)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "palindromeRightArm" and "complementedPalindromeRightArm" generics
### and methods.
###

setGeneric("palindromeRightArm", signature="x",
    function(x, max.mismatch=0, ...) standardGeneric("palindromeRightArm")
)
setGeneric("complementedPalindromeRightArm", signature="x",
    function(x, max.mismatch=0, ...) standardGeneric("complementedPalindromeRightArm")
)

setMethod("palindromeRightArm", "BString",
    function(x, max.mismatch=0, ...)
    {
        start <- nchar(x) - palindromeArmLength(x, max.mismatch=max.mismatch, ...) + 1L
        BString.substr(x, start, nchar(x))
    }
)
setMethod("complementedPalindromeRightArm", "DNAString",
    function(x, max.mismatch=0, ...)
    {
        start <- nchar(x) - complementedPalindromeArmLength(x, max.mismatch=max.mismatch, ...) + 1L
        BString.substr(x, start, nchar(x))
    }
)

setMethod("palindromeRightArm", "BStringViews",
    function(x, max.mismatch=0, ...)
    {
        arm_end <- end(x)
        arm_start <- arm_end - palindromeArmLength(x, max.mismatch=max.mismatch, ...) + 1L
        new("BStringViews", subject=subject(x),
            start=arm_start, end=arm_end, check.views=FALSE)
    }
)
setMethod("complementedPalindromeRightArm", "BStringViews",
    function(x, max.mismatch=0, ...)
    {
        arm_end <- end(x)
        arm_start <- arm_end - complementedPalindromeArmLength(x, max.mismatch=max.mismatch, ...) + 1L
        new("BStringViews", subject=subject(x),
            start=arm_start, end=arm_end, check.views=FALSE)
    }
)

