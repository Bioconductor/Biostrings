### =========================================================================
### Position Weight Matrix (PWM) functions
### -------------------------------------------------------------------------

### A Position Weight Matrix (PWM) is represented as an ordinary matrix.
### We don't use an S4 class for this, not even an S3 class.
.normargPwm <- function(pwm, argname="pwm")
{
    if (!is.numeric(pwm) || !is.matrix(pwm))
        stop("'", argname, "' must be a numeric matrix")
    if (!identical(rownames(pwm), DNA_BASES))
        stop("'rownames(", argname, ")' must be the 4 DNA bases ('DNA_BASES')")
    if (!is.double(pwm))
        storage.mode(pwm) <- "double"
    if (any(is.na(pwm)))
        stop("'", argname, "' contains NAs")
    pwm
}

### A Position Frequency Matrix (PFM) is also represented as an ordinary
### matrix.
.normargPfm <- function(x)
{
    if (!is.numeric(x) || !is.matrix(x))
        stop("'x' must be a numeric matrix")
    ## Check the row names.
    if (is.null(rownames(x)))
        stop("'x' must have row names")
    if (!all(rownames(x) %in% DNA_ALPHABET))
        stop("'rownames(x)' must be a subset of 'DNA_ALPHABET'")
    if (!all(DNA_BASES %in% rownames(x)))
        stop("'rownames(x)' must contain the 4 DNA bases ('DNA_BASES')")
    if (any(duplicated(rownames(x))))
        stop("'x' has duplicated row names")
    ## Check the values.
    if (any(is.na(x)) || any(x < 0))
        stop("frequencies cannot be NA or negative")
    if (any(x[!(rownames(x) %in% DNA_BASES), ] != 0))
        stop("frequencies of IUPAC ambiguity letters must be 0")
    x <- x[DNA_BASES, , drop=FALSE]
    if (!is.double(x))
        storage.mode(x) <- "double"
    csums <- colSums(x)
    ## TODO: Make IRanges::isConstant() a generic function and define a method
    ## for "double" to use in this test.
    if (length(csums) != 0L && !isTRUE(all.equal(min(csums), max(csums))))
        stop("all columns in 'x' must sum to the same value")
    x
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Some utilities to operate on a PWM.
###

### Extracts the max weight for each position (i.e. column) in a PWM.
setGeneric("maxWeights", function(x) standardGeneric("maxWeights"))

setMethod("maxWeights", "matrix",
    function(x)
    {
        x <- .normargPwm(x, argname="x")
        #sapply(seq_len(ncol(x)), function(i) max(x[ , i]))
        ## This will be faster than the above on large matrices
        do.call(pmax, lapply(seq_len(nrow(x)), function(i) x[i, ]))
    })

### Extracts the min weight for each position (i.e. column) in a PWM.
setGeneric("minWeights", function(x) standardGeneric("minWeights"))

setMethod("minWeights", "matrix",
    function(x)
    {
        x <- .normargPwm(x, argname="x")
        #sapply(seq_len(ncol(x)), function(i) max(x[ , i]))
        ## This will be faster than the above on large matrices
        do.call(pmin, lapply(seq_len(nrow(x)), function(i) x[i, ]))
    })

### Computes the highest possible score that can be obtained with a PWM.
setGeneric("maxScore", function(x) standardGeneric("maxScore"))
setMethod("maxScore", "ANY", function(x) sum(maxWeights(x)))

### Computes the lowest possible score that can be obtained with a PWM.
setGeneric("minScore", function(x) standardGeneric("minScore"))
setMethod("minScore", "ANY", function(x) sum(minWeights(x)))

### TODO: There is no reason to treat this differently than the above
### utilities. So either this should be implemented as a generic+method
### or the above utilities should be implemented as ordinary functions.
unitScale <- function(x)
{
    minS <- minScore(x)
    (x - minS/ncol(x)) / (maxScore(x) - minS)
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "PWM" generic and methods.
###

setGeneric("PWM", signature="x",
    function(x, type = c("log2probratio", "prob"),
             prior.params = c(A=0.25, C=0.25, G=0.25, T=0.25))
        standardGeneric("PWM")
)

setMethod("PWM", "character",
    function(x, type = c("log2probratio", "prob"),
             prior.params = c(A=0.25, C=0.25, G=0.25, T=0.25))
        PWM(DNAStringSet(x), type = type, prior.params = prior.params)
)

setMethod("PWM", "DNAStringSet",
    function(x, type = c("log2probratio", "prob"),
             prior.params = c(A=0.25, C=0.25, G=0.25, T=0.25))
    {
        if (!isConstant(width(x)))
            stop("'x' must be rectangular (i.e. have a constant width)")
        cmat <- consensusMatrix(x)
        rsums <- rowSums(cmat)
        if (any(rsums[!(names(rsums) %in% DNA_BASES)] > 0))
            stop("'x' contains non 'DNA_BASES' letters")
        cmat <- cmat[DNA_BASES, , drop=FALSE]
        PWM(cmat, type = type, prior.params = prior.params)
    }
)

### Assumes 'x' is a Position *Frequency* Matrix (PFM) and computes the
### corresponding Position *Weight* Matrix (PWM).
setMethod("PWM", "matrix",
    function(x, type = c("log2probratio", "prob"),
             prior.params = c(A=0.25, C=0.25, G=0.25, T=0.25))
    {
        x <- .normargPfm(x)
        type <- match.arg(type)
        if (!is.numeric(prior.params))
            stop("'prior.params' must be a numeric vector")
        prior.params <- prior.params[DNA_BASES]
        if (any(is.na(prior.params) | prior.params < 0))
            stop("'prior.params' must have non-negative named elements corresponding to 'DNA_BASES'")
        priorN <- sum(prior.params)
        postProbs <- (x + prior.params) / (length(x) + priorN)
        if (type == "log2probratio") {
            priorProbs <- prior.params / priorN
            if (any(priorProbs == 0))
                stop("infinite values in PWM due to 0's in 'prior.params'")
            ans <- log2(postProbs / priorProbs)
        } else {
            ans <- postProbs
        }
        unitScale(ans)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "PWMscoreStartingAt" function.
###

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

### TODO: There is no reason to treat this differently than anything else
### in this file. So maybe this should be implemented as a generic+method
### too (or matchPWM/countPWM shouldn't).
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
    .Call("PWM_score_starting_at",
          pwm, subject, base_codes, starting.at,
          PACKAGE="Biostrings")
}

PWMscore <- function(pwm, subject, start=1)
    .Defunct('PWMscoreStartingAt')


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "matchPWM" generic and methods.
###

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

### pwm: the Position Weight Matrix (numeric matrix with row names A, C, G
###      and T);
### subject: a DNAString object containing the subject sequence;
### min.score: given as a percentage (e.g. "90%") of the highest possible
###            score or as a single number.
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

### Note the dispatch on 'subject'.
setGeneric("matchPWM", signature="subject",
    function(pwm, subject, min.score="80%", ...)
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
### a normal XStringViews object).
### matchPWM does not support "out of limits"  matches.
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

### Note the dispatch on 'subject'.
setGeneric("countPWM", signature="subject",
    function(pwm, subject, min.score="80%", ...)
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

