### =========================================================================
### findPalindromes() and family
### -------------------------------------------------------------------------


.get_DNAorRNA_palindrome_L2R_lkup <- function()
{
    keys <- unname(DNA_CODES[c(DNA_BASES, "-")])
    vals <- unname(DNA_CODES[c(rev(DNA_BASES), "-")])
    buildLookupTable(keys, vals)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The findPalindromes() generic and methods
###

### Return an IRanges object.
.find_palindromes <- function(subject, min.armlength,
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
    C_ans <- .Call2("find_palindromes",
                    subject,
                    min.armlength, max.looplength, max.mismatch,
                    L2R_lkup,
                    PACKAGE="Biostrings")
    unsafe.newXStringViews(subject, start(C_ans), width(C_ans))
}

setGeneric("findPalindromes", signature="subject",
    function(subject, min.armlength=4,
                      max.looplength=1, min.looplength=0,
                      max.mismatch=0)
        standardGeneric("findPalindromes")
)

setMethod("findPalindromes", "XString",
    function(subject, min.armlength=4,
                      max.looplength=1, min.looplength=0,
                      max.mismatch=0)
    {
        .find_palindromes(subject, min.armlength,
                          max.looplength, min.looplength,
                          max.mismatch, NULL)
    }
)

setMethod("findPalindromes", "DNAString",
    function(subject, min.armlength=4,
                      max.looplength=1, min.looplength=0,
                      max.mismatch=0)
    {
        L2R_lkup <- .get_DNAorRNA_palindrome_L2R_lkup()
        .find_palindromes(subject, min.armlength,
                          max.looplength, min.looplength,
                          max.mismatch, L2R_lkup)
    }
)

setMethod("findPalindromes", "RNAString",
    function(subject, min.armlength=4,
                      max.looplength=1, min.looplength=0,
                      max.mismatch=0)
    {
        L2R_lkup <- .get_DNAorRNA_palindrome_L2R_lkup()
        .find_palindromes(subject, min.armlength,
                          max.looplength, min.looplength,
                          max.mismatch, L2R_lkup)
    }
)

setMethod("findPalindromes", "XStringViews",
    function(subject, min.armlength=4,
                      max.looplength=1, min.looplength=0,
                      max.mismatch=0)
    {
        tmp <- vector(mode="list", length=length(subject))
        offsets <- start(subject) - 1L
        for (i in seq_along(subject)) {
            pals <- findPalindromes(subject[[i]],
                                    min.armlength=min.armlength,
                                    max.looplength=max.looplength,
                                    min.looplength=min.looplength,
                                    max.mismatch=max.mismatch)
            tmp[[i]] <- shift(ranges(pals), shift=offsets[i])
        }
        ans_ranges <- do.call("c", tmp)
        unsafe.newXStringViews(subject(subject),
                               start(ans_ranges), width(ans_ranges))
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The palindromeArmLength() generic and methods
###

.palindrome_arm_length <- function(x, max.mismatch, L2R_lkup)
{
    max.mismatch <- normargMaxMismatch(max.mismatch)
    armlength <- .Call2("palindrome_arm_length",
                        x, max.mismatch, L2R_lkup,
                        PACKAGE="Biostrings")
    if (armlength == 0L)
        stop("'x' is not a palindrome (no arms found)")
    armlength
}

setGeneric("palindromeArmLength", signature="x",
    function(x, max.mismatch=0, ...)
        standardGeneric("palindromeArmLength")
)

setMethod("palindromeArmLength", "XString",
    function(x, max.mismatch=0, ...)
    {
        .palindrome_arm_length(x, max.mismatch, NULL)
    }
)

setMethod("palindromeArmLength", "DNAString",
    function(x, max.mismatch=0, ...)
    {
        L2R_lkup <- .get_DNAorRNA_palindrome_L2R_lkup()
        .palindrome_arm_length(x, max.mismatch, L2R_lkup)
    }
)

setMethod("palindromeArmLength", "RNAString",
    function(x, max.mismatch=0, ...)
    {
        L2R_lkup <- .get_DNAorRNA_palindrome_L2R_lkup()
        .palindrome_arm_length(x, max.mismatch, L2R_lkup)
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The palindromeLeftArm() and palindromeRightArm() generics and methods
###

setGeneric("palindromeLeftArm", signature="x",
    function(x, max.mismatch=0, ...)
        standardGeneric("palindromeLeftArm")
)

setGeneric("palindromeRightArm", signature="x",
    function(x, max.mismatch=0, ...)
        standardGeneric("palindromeRightArm")
)

setMethod("palindromeLeftArm", "XString",
    function(x, max.mismatch=0, ...)
        subseq(x,
            start=1L,
            end=palindromeArmLength(x, max.mismatch=max.mismatch, ...)
        )
)

setMethod("palindromeRightArm", "XString",
    function(x, max.mismatch=0, ...)
    {
        start <- nchar(x) - palindromeArmLength(x,
                                max.mismatch=max.mismatch, ...) + 1L
        subseq(x, start=start, end=nchar(x))
    }
)

setMethod("palindromeLeftArm", "XStringViews",
    function(x, max.mismatch=0, ...)
    {
        ans_start <- start(x)
        ans_width <- palindromeArmLength(x, max.mismatch=max.mismatch, ...)
        unsafe.newXStringViews(subject(x), ans_start, ans_width)
    }
)

setMethod("palindromeRightArm", "XStringViews",
    function(x, max.mismatch=0, ...)
    {
        ans_width <- palindromeArmLength(x,
                         max.mismatch=max.mismatch, ...)
        ans_start <- end(x) - ans_width + 1L
        unsafe.newXStringViews(subject(x), ans_start, ans_width)
    }
)

