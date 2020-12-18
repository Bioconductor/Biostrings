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
                              max.mismatch, allow.wobble, L2R_lkup)
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
    if (!isTRUEorFALSE(allow.wobble))
        stop("'allow.wobble' must be a logical")
    if (allow.wobble && !(seqtype(subject) %in% c("DNA", "RNA")))
        stop("subject must be DNA or RNA if 'allow.wobble' is TRUE")
    ## check max.mismatch
    max.mismatch <- normargMaxMismatch(max.mismatch)
    C_ans <- .Call2("find_palindromes",
                    subject,
                    min.armlength, max.looplength, max.mismatch,
                    min.looplength, allow.wobble, L2R_lkup,
                    PACKAGE="Biostrings")
    unsafe.newXStringViews(subject, start(C_ans), width(C_ans))
}

setGeneric("findPalindromes", signature="subject",
    function(subject, min.armlength=4,
                      max.looplength=1, min.looplength=0,
                      max.mismatch=0, allow.wobble=FALSE)
        standardGeneric("findPalindromes")
)

setMethod("findPalindromes", "XString",
    function(subject, min.armlength=4,
                      max.looplength=1, min.looplength=0,
                      max.mismatch=0, allow.wobble=FALSE)
    {
        .find_palindromes(subject, min.armlength,
                          max.looplength, min.looplength,
                          max.mismatch, allow.wobble,
                          NULL)
    }
)

setMethod("findPalindromes", "DNAString",
    function(subject, min.armlength=4,
                      max.looplength=1, min.looplength=0,
                      max.mismatch=0, allow.wobble=FALSE)
    {
        L2R_lkup <- .get_DNAorRNA_palindrome_L2R_lkup()
        .find_palindromes(subject, min.armlength,
                          max.looplength, min.looplength,
                          max.mismatch, allow.wobble,
                          L2R_lkup)
    }
)

setMethod("findPalindromes", "RNAString",
    function(subject, min.armlength=4,
                      max.looplength=1, min.looplength=0,
                      max.mismatch=0, allow.wobble=FALSE)
    {
        L2R_lkup <- .get_DNAorRNA_palindrome_L2R_lkup()
        .find_palindromes(subject, min.armlength,
                          max.looplength, min.looplength,
                          max.mismatch, allow.wobble,
                          L2R_lkup)
    }
)

setMethod("findPalindromes", "XStringViews",
    function(subject, min.armlength=4,
                      max.looplength=1, min.looplength=0,
                      max.mismatch=0, allow.wobble=FALSE)
    {
        tmp <- vector(mode="list", length=length(subject))
        offsets <- start(subject) - 1L
        for (i in seq_along(subject)) {
            pals <- findPalindromes(subject[[i]],
                                    min.armlength=min.armlength,
                                    max.looplength=max.looplength,
                                    min.looplength=min.looplength,
                                    max.mismatch=max.mismatch,
                                    allow.wobble=allow.wobble)
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
                      max.mismatch=0, allow.wobble=FALSE)
    {
        findPalindromes(toXStringViewsOrXString(subject),
                        min.armlength=min.armlength,
                        max.looplength=max.looplength,
                        min.looplength=min.looplength,
                        max.mismatch=max.mismatch,
                        allow.wobble=allow.wobble)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The palindromeArmLength() generic and methods
###

.palindrome_arm_length <- function(x, max.mismatch, allow.wobble, L2R_lkup)
{
    if (!isTRUEorFALSE(allow.wobble))
        stop("'allow.wobble' must be a logical")
    if (allow.wobble && !(seqtype(x) %in% c("DNA", "RNA")))
        stop("x must be DNA or RNA if 'allow.wobble' is TRUE")
    max.mismatch <- normargMaxMismatch(max.mismatch)
    .Call2("palindrome_arm_length",
           x, max.mismatch, allow.wobble, L2R_lkup,
           PACKAGE="Biostrings")
}

setGeneric("palindromeArmLength", signature="x",
    function(x, max.mismatch=0, allow.wobble=FALSE) standardGeneric("palindromeArmLength")
)

setMethod("palindromeArmLength", "XString",
    function(x, max.mismatch=0, allow.wobble=FALSE)
    {
        x <- as(x, "XStringSet")
        .palindrome_arm_length(x, max.mismatch, allow.wobble, NULL)
    }
)

setMethod("palindromeArmLength", "DNAString",
    function(x, max.mismatch=0, allow.wobble=FALSE)
    {
        L2R_lkup <- .get_DNAorRNA_palindrome_L2R_lkup()
        x <- as(x, "DNAStringSet")
        .palindrome_arm_length(x, max.mismatch, allow.wobble, L2R_lkup)
    }
)

setMethod("palindromeArmLength", "RNAString",
    function(x, max.mismatch=0, allow.wobble=FALSE)
    {
        L2R_lkup <- .get_DNAorRNA_palindrome_L2R_lkup()
        x <- as(x, "RNAStringSet")
        .palindrome_arm_length(x, max.mismatch, allow.wobble, L2R_lkup)
    }
)

setMethod("palindromeArmLength", "XStringViews",
    function(x, max.mismatch=0, allow.wobble=FALSE)
    {
        if (length(x) == 0)
            return(integer(0))
        if (seqtype(x) %in% c("DNA", "RNA")) {
            L2R_lkup <- .get_DNAorRNA_palindrome_L2R_lkup()
        } else {
            L2R_lkup <- NULL
        }
        max.mismatch <- normargMaxMismatch(max.mismatch)
        x <- fromXStringViewsToStringSet(x, out.of.limits="error")
        .palindrome_arm_length(x, max.mismatch, allow.wobble, L2R_lkup)
    }
)

setMethod("palindromeArmLength", "XStringSet",
    function(x, max.mismatch=0, allow.wobble=FALSE)
    {
        if (length(x) == 0)
            return(integer(0))
        if (seqtype(x) %in% c("DNA", "RNA")) {
            L2R_lkup <- .get_DNAorRNA_palindrome_L2R_lkup()
        } else {
            L2R_lkup <- NULL
        }
        max.mismatch <- normargMaxMismatch(max.mismatch)
        .palindrome_arm_length(x, max.mismatch, allow.wobble, L2R_lkup)
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The palindromeLeftArm() and palindromeRightArm() generics and methods
###

setGeneric("palindromeLeftArm", signature="x",
    function(x, max.mismatch=0, allow.wobble=FALSE)
        standardGeneric("palindromeLeftArm")
)

setGeneric("palindromeRightArm", signature="x",
    function(x, max.mismatch=0, allow.wobble=FALSE)
        standardGeneric("palindromeRightArm")
)

setMethod("palindromeLeftArm", "XString",
    function(x, max.mismatch=0, allow.wobble=FALSE)
        subseq(x,
            start=1L,
            end=palindromeArmLength(x, max.mismatch=max.mismatch, allow.wobble=allow.wobble)
        )
)

setMethod("palindromeRightArm", "XString",
    function(x, max.mismatch=0, allow.wobble=FALSE)
    {
        start <- nchar(x) -
                 palindromeArmLength(x, max.mismatch=max.mismatch, allow.wobble=allow.wobble) +
                 1L
        subseq(x, start=start, end=nchar(x))
    }
)

setMethod("palindromeLeftArm", "XStringViews",
    function(x, max.mismatch=0, allow.wobble=FALSE)
    {
        ans_start <- start(x)
        ans_width <- palindromeArmLength(x, max.mismatch=max.mismatch, allow.wobble=allow.wobble)
        unsafe.newXStringViews(subject(x), ans_start, ans_width)
    }
)

setMethod("palindromeRightArm", "XStringViews",
    function(x, max.mismatch=0, allow.wobble=FALSE)
    {
        ans_width <- palindromeArmLength(x, max.mismatch=max.mismatch, allow.wobble=allow.wobble)
        ans_start <- end(x) - ans_width + 1L
        unsafe.newXStringViews(subject(x), ans_start, ans_width)
    }
)

