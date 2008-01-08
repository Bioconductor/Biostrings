### =========================================================================
### The findPalindromes() generic & related functions
### -------------------------------------------------------------------------


.normalize.anti <- function(anti, class)
{
    if (class %in% c("DNAString", "RNAString")) {
        if (is.null(anti))
            return(TRUE)
        if (!is.logical(anti) || length(anti) != 1 || is.na(anti))
            stop("'anti' can only be 'NULL', 'TRUE' or 'FALSE'")
        return(anti)
    }
    if (is.null(anti))
        return(FALSE)
    stop("'anti' must be 'NULL' with a non-DNAString (or non-RNAString) sequence")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "findPalindromes" generic and methods.
###

debug_find_palindromes <- function()
{
    invisible(.Call("find_palindromes_debug", PACKAGE="Biostrings"))
}

### Return a list with the "start" and the "end" components.
.find.palindromes <- function(subject, min.armlength, max.ngaps, L2R_lkup)
{
    
    views <- .Call("find_palindromes",
                   subject@data@xp, subject@offset, subject@length,
                   min.armlength, max.ngaps, L2R_lkup,
                   PACKAGE="Biostrings")
    new("BStringViews", subject=subject,
        start=views$start, end=views$end, check.views=FALSE)
}

setGeneric("findPalindromes", signature="subject",
    function(subject, min.armlength=4, max.ngaps=1, anti=NULL)
        standardGeneric("findPalindromes")
)

setMethod("findPalindromes", "BString",
    function(subject, min.armlength=4, max.ngaps=1, anti=NULL)
    {
        if (!is.numeric(min.armlength) || length(min.armlength) != 1 || is.na(min.armlength))
            stop("'min.armlength' must be a single integer")
        min.armlength <- as.integer(min.armlength)
        if (min.armlength < 2)
            stop("'min.armlength' must be >= 2")
        if (!is.numeric(max.ngaps) || length(max.ngaps) != 1 || is.na(max.ngaps))
            stop("'max.ngaps' must be a single integer")
        max.ngaps <- as.integer(max.ngaps)
        if (max.ngaps < 0)
            stop("'max.ngaps' must be a non-negative integer")
        if (.normalize.anti(anti, class(subject)))
            L2R_lkup <- getDNAComplementLookup()
        else
            L2R_lkup <- NULL
        .find.palindromes(subject, min.armlength, max.ngaps, L2R_lkup)
    }
)

### WARNING: Unlike with the "findPalindromes" method for BString objects, the
### BStringViews object returned by this method is not guaranteed to have its
### views ordered from left to right! One important particular case where this
### is guaranteed though is when 'subject' is a normalized BStringViews object.
setMethod("findPalindromes", "BStringViews",
    function(subject, min.armlength=4, max.ngaps=1, anti=NULL)
    {
        ans_start <- ans_end <- integer(0)
        for (i in seq_len(length(subject))) {
            pals <- findPalindromes(subject[[i]], min.armlength, max.ngaps, anti)
            offset <- start(subject)[i] - 1L
            ans_start <- c(ans_start, offset + start(pals))
            ans_end <- c(ans_end, offset + end(pals))
        }
        new("BStringViews", subject=subject(subject),
            start=ans_start, end=ans_end, check.views=FALSE)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "palindromeArmLength" generic and methods.
###

setGeneric("palindromeArmLength", signature="x",
    function(x, anti=NULL, ...) standardGeneric("palindromeArmLength")
)

setMethod("palindromeArmLength", "BString",
    function(x, anti=NULL, ...)
    {
        if (.normalize.anti(anti, class(x)))
            revx <- reverseComplement(x)
        else
            revx <- reverse(x)
        armlength <- lcprefix(x, revx)
        if (armlength == 0L)
            stop("'x' is not a palindrome (no arms found)")
        armlength
    }
)

setMethod("palindromeArmLength", "BStringViews",
    function(x, anti=NULL, ...)
    {
        if (length(x) == 0)
            return(integer(0))
        sapply(seq_len(length(x)), function(i) palindromeArmLength(x[[i]], anti=anti))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "palindromeLeftArm" generic and methods.
###

setGeneric("palindromeLeftArm", signature="x",
    function(x, anti=NULL, ...) standardGeneric("palindromeLeftArm")
)

setMethod("palindromeLeftArm", "BString",
    function(x, anti=NULL, ...)
    {
        BString.substr(x, 1L, palindromeArmLength(x, anti=anti))
    }
)

setMethod("palindromeLeftArm", "BStringViews",
    function(x, anti=NULL, ...)
    {
        arm_start <- start(x)
        arm_end <- arm_start + palindromeArmLength(x, anti=anti) - 1L
        new("BStringViews", subject=subject(x),
            start=arm_start, end=arm_end, check.views=FALSE)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "palindromeRightArm" generic and methods.
###

setGeneric("palindromeRightArm", signature="x",
    function(x, anti=NULL, ...) standardGeneric("palindromeRightArm")
)

setMethod("palindromeRightArm", "BString",
    function(x, anti=NULL, ...)
    {
        BString.substr(x, nchar(x) - palindromeArmLength(x, anti=anti) + 1L, nchar(x))
    }
)

setMethod("palindromeRightArm", "BStringViews",
    function(x, anti=NULL, ...)
    {
        arm_end <- end(x)
        arm_start <- arm_end - palindromeArmLength(x, anti=anti) + 1L
        new("BStringViews", subject=subject(x),
            start=arm_start, end=arm_end, check.views=FALSE)
    }
)

