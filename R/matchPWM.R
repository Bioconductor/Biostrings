### =========================================================================
### Position Weight Matrix (PWM) functions
### -------------------------------------------------------------------------

### A Position Weight Matrix (PWM) is represented as an ordinary matrix.
### We don't use an S4 class for this, not even an S3 class.
.normargPwm <- function(pwm, argname="pwm")
{
    if (!is.matrix(pwm) || !is.numeric(pwm))
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
### matrix. Unlike a PWM, it must be of type integer (it will typically be
### the result of consensusMatrix()).
.normargPfm <- function(x)
{
    if (!is.matrix(x) || !is.integer(x))
        stop("invalid PFM 'x': not an integer matrix")
    ## Check the row names.
    if (is.null(rownames(x)))
        stop("invalid PFM 'x': no row names")
    if (!all(rownames(x) %in% DNA_ALPHABET))
        stop("invalid PFM 'x': row names must be in 'DNA_ALPHABET'")
    if (!all(DNA_BASES %in% rownames(x)))
        stop("invalid PFM 'x': row names must contain A, C, G and T")
    if (any(duplicated(rownames(x))))
        stop("invalid PFM 'x': duplicated row names")
    ## Check the nb of cols.
    if (ncol(x) == 0L)
        stop("invalid PFM 'x': no columns")
    ## Check the values.
    if (any(is.na(x)) || any(x < 0L))
        stop("invalid PFM 'x': values cannot be NA or negative")
    if (any(x[!(rownames(x) %in% DNA_BASES), ] != 0L))
        stop("invalid PFM 'x': IUPAC ambiguity letters are represented")
    x <- x[DNA_BASES, , drop=FALSE]
    if (!isConstant(colSums(x)))
        stop("invalid PFM 'x': all columns in 'x' must sum to the same ",
             "value.\n  If the PFM was obtained by calling consensusMatrix() ",
             "on a DNAStringSet\n  object, please make sure that this object ",
             "is rectangular (i.e. has a\n  constant width).")
    x
}

### Typical 'prior.params' vector: c(A=0.25, C=0.25, G=0.25, T=0.25)
.normargPriorParams <- function(prior.params)
{
    if (!is.numeric(prior.params))
        stop("'prior.params' must be a numeric vector")
    if (length(prior.params) != length(DNA_BASES) ||
        !setequal(names(prior.params), DNA_BASES))
        stop("'prior.params' elements must be named A, C, G and T")
    ## Re-order the elements.
    prior.params <- prior.params[DNA_BASES]
    if (any(is.na(prior.params)) || any(prior.params < 0))
        stop("'prior.params' contains NAs and/or negative values")
    prior.params
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
    {
        dnaset <- DNAStringSet(x)
        PWM(dnaset, type = type, prior.params = prior.params)
    }
)

setMethod("PWM", "DNAStringSet",
    function(x, type = c("log2probratio", "prob"),
             prior.params = c(A=0.25, C=0.25, G=0.25, T=0.25))
    {
        if (!isConstant(width(x)))
            stop("'x' must be rectangular (i.e. have a constant width)")
        pfm <- consensusMatrix(x)
        PWM(pfm, type = type, prior.params = prior.params)
    }
)

### Assumes 'x' is a Position *Frequency* Matrix (PFM) and computes the
### corresponding Position *Weight* Matrix (PWM).
setMethod("PWM", "matrix",
    function(x, type = c("log2probratio", "prob"),
             prior.params = c(A=0.25, C=0.25, G=0.25, T=0.25))
    {
        x <- .normargPfm(x)
        ## From here 'x' is guaranteed to have at least 1 column and to have
        ## all its columns sum to the same value.
        nseq <- sum(x[ , 1L])
        type <- match.arg(type)
        prior.params <- .normargPriorParams(prior.params)
        priorN <- sum(prior.params)
        ## NOTE (H.P.): What's the purpose of dividing by nseq + priorN here?
        ## It won't have any impact on the final result (because of the
        ## unitScale final step).
        postProbs <- (x + prior.params) / (nseq + priorN)
        if (type == "log2probratio") {
            if (any(prior.params == 0))
                stop("infinite values in PWM due to 0's in 'prior.params'")
            prior.probs <- prior.params / priorN
            ans <- log2(postProbs / prior.probs)
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
    if (is.character(subject)) {
        subject <- DNAString(subject)
    } else if (is(subject, "MaskedDNAString")) {
        if (any(active(masks(subject))))
            stop("active masks are not supported yet, please complain!")
        subject <- unmasked(subject)
    } else if (is(subject, "XStringViews")) {
        subject <- subject(subject)
    }
    if (!is(subject, "DNAString"))
        stop("'subject' must be a single character string, ",
             "or a DNAString object, or a MaskedDNAString object with ",
             "no active masks, or a Views object on a DNAString subject")
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
    .Call2("PWM_score_starting_at",
          pwm, subject, starting.at, base_codes,
          PACKAGE="Biostrings")
}


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
.XString.matchPWM <- function(pwm, subject, min.score,
                              with.score=FALSE, count.only=FALSE)
{
    ## checking 'pwm'
    pwm <- .normargPwm(pwm)
    ## checking 'min.score'
    min.score <- .normargMinScore(min.score, pwm)
    ## checking 'with.score'
    if (!isTRUEorFALSE(with.score))
        stop("'with.score' must be TRUE or FALSE")
    ## no need to check 'count.only' (not a user controlled argument)

    base_codes <- xscodes(subject, baseOnly=TRUE)
    C_ans <- .Call2("XString_match_PWM",
                   pwm, subject, min.score, count.only, base_codes,
                   PACKAGE="Biostrings")
    if (count.only)
        return(C_ans)
    ans_start <- start(C_ans)
    ans_width <- width(C_ans)
    ans <- unsafe.newXStringViews(subject, ans_start, ans_width)
    if (with.score) {
        score <- PWMscoreStartingAt(pwm, subject, starting.at=ans_start)
        mcols(ans) <- DataFrame(score=score)
    }
    ans
}

.XStringViews.matchPWM <- function(pwm, subject, min.score,
                                   with.score=FALSE, count.only=FALSE)
{
    ## checking 'pwm'
    pwm <- .normargPwm(pwm)
    ## checking 'min.score'
    min.score <- .normargMinScore(min.score, pwm)
    ## checking 'with.score'
    if (!isTRUEorFALSE(with.score))
        stop("'with.score' must be TRUE or FALSE")
    ## no need to check 'count.only' (not a user controlled argument)

    subject0 <- subject(subject)
    base_codes <- xscodes(subject0, baseOnly=TRUE)
    C_ans <- .Call2("XStringViews_match_PWM",
                   pwm,
                   subject0, start(subject), width(subject),
                   min.score, count.only, base_codes,
                   PACKAGE="Biostrings")
    if (count.only)
        return(C_ans)
    ans_start <- start(C_ans)
    ans_width <- width(C_ans)
    ans <- unsafe.newXStringViews(subject0, ans_start, ans_width)
    if (with.score) {
        score <- PWMscoreStartingAt(pwm, subject0, starting.at=ans_start)
        mcols(ans) <- DataFrame(score=score)
    }
    ans
}

### Note the dispatch on 'subject'.
setGeneric("matchPWM", signature="subject",
    function(pwm, subject, min.score="80%", with.score=FALSE, ...)
        standardGeneric("matchPWM")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("matchPWM", "character",
    function(pwm, subject, min.score="80%", with.score=FALSE)
        matchPWM(pwm, DNAString(subject),
                 min.score=min.score, with.score=with.score)
)

### Dispatch on 'subject' (see signature of generic).
setMethod("matchPWM", "DNAString",
    function(pwm, subject, min.score="80%", with.score=FALSE)
        .XString.matchPWM(pwm, subject, min.score, with.score=with.score)
)

### Dispatch on 'subject' (see signature of generic).
### WARNING: Unlike the other "matchPWM" methods, the XStringViews object
### returned by this method is not guaranteed to have its views ordered from
### left to right in general! One important particular case where this is
### guaranteed though is when 'isNormal(subject)' is TRUE (i.e. 'subject' is
### a normal XStringViews object).
### matchPWM does not support "out of limits"  matches.
setMethod("matchPWM", "XStringViews",
    function(pwm, subject, min.score="80%", with.score=FALSE)
        .XStringViews.matchPWM(pwm, subject, min.score, with.score=with.score)
)

### Dispatch on 'subject' (see signature of generic).
setMethod("matchPWM", "MaskedDNAString",
    function(pwm, subject, min.score="80%", with.score=FALSE)
        matchPWM(pwm, toXStringViewsOrXString(subject),
                 min.score=min.score, with.score=with.score)
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

