### =========================================================================
### The matchPWM() generic & related functions
### -------------------------------------------------------------------------


.normargPwm <- function(pwm)
{
    if (!is.numeric(pwm) || !is.matrix(pwm))
        stop("'pwm' must be a numeric matrix")
    if (!identical(rownames(pwm), DNA_BASES))
        stop("'rownames(pwm)' must be the DNA bases ('DNA_BASES')")
    if (!is.double(pwm))
        storage.mode(pwm) <- "double"
    if (any(is.na(pwm)))
        stop("'pwm' contains NAs")
    pwm
}

.normargSubject <- function(subject)
{
    if (is(subject, "MaskedDNAString")) {
        if (any(active(masks(subject))))
            stop("active masks are not supported yet, please complain!")
        subject <- unmasked(subject)
    } else if (!is(subject, "DNAString"))
        stop("'subject' must be a DNAString or MaskedDNAString object ",
             "with no active masks")
    subject
}

### Utility function for getting the max weight for each position in the PWM
maxWeights <- function(pwm)
{
    pwm <- .normargPwm(pwm)
    #sapply(seq_len(ncol(pwm)), function(i) max(pwm[ , i]))
    ## This will be faster than the above on large matrices
    do.call(pmax, lapply(seq_len(nrow(pwm)), function(i) pwm[i, ]))
}

### Utility function for getting the highest possible score
maxScore <- function(pwm)
{
    sum(maxWeights(pwm))
}

.normargMinScore <- function(min.score, pwm)
{
    if (!isSingleNumber(min.score) && !isSingleString(min.score))
        stop("'min.score' must be a single number or string")
    if (is.numeric(min.score)) {
        if (!is.double(min.score))
            storage.mode(min.score) <- "double"
        return(min.score)
    }
    nc <- nchar(min.score)
    if (substr(min.score, nc, nc) == "%")
        min.score <- substr(min.score, 1L, nc-1L)
    maxScore(pwm) * as.double(min.score) / 100.00
}

### Method needed for searching the minus strand of a chromosome like
### this:
###   > matchPWM(reverseComplement(pwm), chr1)
### Note that the generic function is defined in Biostrings.
setMethod("reverseComplement", "matrix",
    function(x, ...)
    {
        ans <- rev(x)
        attributes(ans) <- attributes(x)
        ans
    }
)

PWMscoreStartingAt <- function(pwm, subject, starting.at=1)
{
    ## checking 'pwm'
    pwm <- .normargPwm(pwm)
    ## checking 'subject'
    subject <- .normargSubject(subject)
    ## checking 'starting.at'
    if (!is.numeric(starting.at))
        stop("'starting.at' must be a vector of integers")
    if (!is.integer(starting.at))
        starting.at <- as.integer(starting.at)

    base_codes <- xscodes(subject, baseOnly=TRUE)
    .Call("PWM_score_starting_at", pwm, subject, base_codes, starting.at, PACKAGE="Biostrings")
}

PWMscore <- function(pwm, subject, start=1)
{
    .Deprecated('PWMscoreStartingAt')
    PWMscoreStartingAt(pwm, subject, starting.at=start)
}

### pwm: the Position Weight Matrix (numeric matrix with row names A, C, G and T)
### subject: a DNAString object containing the subject sequence
### min.score: given as a percentage (e.g. "90%") of the highest possible
###            score or as a single number
.XString.matchPWM <- function(pwm, subject, min.score, count.only=FALSE)
{
    ## checking 'pwm'
    pwm <- .normargPwm(pwm)
    ## checking 'min.score'
    min.score <- .normargMinScore(min.score, pwm)
    ## no need to check 'count.only' (not a user controlled argument)

    base_codes <- xscodes(subject, baseOnly=TRUE)
    C_ans <- .Call("XString_match_PWM",
                   pwm, subject, base_codes, min.score, count.only,
                   PACKAGE="Biostrings")
    if (count.only)
        return(C_ans)
    unsafe.newXStringViews(subject, start(C_ans), width(C_ans))
}

.XStringViews.matchPWM <- function(pwm, subject, min.score, count.only=FALSE)
{
    ## checking 'pwm'
    pwm <- .normargPwm(pwm)
    ## checking 'min.score'
    min.score <- .normargMinScore(min.score, pwm)
    ## no need to check 'count.only' (not a user controlled argument)

    base_codes <- xscodes(subject, baseOnly=TRUE)
    C_ans <- .Call("XStringViews_match_PWM",
                   pwm,
                   subject(subject), start(subject), width(subject),
                   base_codes, min.score, count.only,
                   PACKAGE="Biostrings")
    if (count.only)
        return(C_ans)
    unsafe.newXStringViews(subject(subject), start(C_ans), width(C_ans))
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "matchPWM" generic and methods.
###

setGeneric("matchPWM", signature="subject",
    function(pwm, subject, min.score="80%")
        standardGeneric("matchPWM")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("matchPWM", "character",
    function(pwm, subject, min.score="80%")
        matchPWM(pwm, DNAString(subject), min.score)
)

### Dispatch on 'subject' (see signature of generic).
setMethod("matchPWM", "DNAString",
    function(pwm, subject, min.score="80%")
        .XString.matchPWM(pwm, subject, min.score)
)

### Dispatch on 'subject' (see signature of generic).
### WARNING: Unlike the other "matchPWM" methods, the XStringViews object
### returned by this method is not guaranteed to have its views ordered from
### left to right in general! One important particular case where this is
### guaranteed though is when 'isNormal(subject)' is TRUE (i.e. 'subject' is
### a normal XStringViews object) and 'max.mismatch=0' (no "out of limits"
### matches).
setMethod("matchPWM", "XStringViews",
    function(pwm, subject, min.score="80%")
        .XStringViews.matchPWM(pwm, subject, min.score)
)

### Dispatch on 'subject' (see signature of generic).
setMethod("matchPWM", "MaskedDNAString",
    function(pwm, subject, min.score="80%")
        matchPWM(pwm, toXStringViewsOrXString(subject), min.score)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "countPWM" generic and methods.
###

setGeneric("countPWM", signature="subject",
    function(pwm, subject, min.score="80%")
        standardGeneric("countPWM")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("countPWM", "character",
    function(pwm, subject, min.score="80%")
        countPWM(pwm, DNAString(subject), min.score)
)

### Dispatch on 'subject' (see signature of generic).
setMethod("countPWM", "DNAString",
    function(pwm, subject, min.score="80%")
        .XString.matchPWM(pwm, subject, min.score, count.only=TRUE)
)

### Dispatch on 'subject' (see signature of generic).
setMethod("countPWM", "XStringViews",
    function(pwm, subject, min.score="80%")
        .XStringViews.matchPWM(pwm, subject, min.score, count.only=TRUE)
)

### Dispatch on 'subject' (see signature of generic).
setMethod("countPWM", "MaskedDNAString",
    function(pwm, subject, min.score="80%")
        countPWM(pwm, toXStringViewsOrXString(subject), min.score)
)
