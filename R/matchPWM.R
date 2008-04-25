### =========================================================================
### The matchPWM() generic & related functions
### -------------------------------------------------------------------------


.normalize.pwm <- function(pwm)
{
    if (!is.numeric(pwm) || !is.matrix(pwm))
        stop("'pwm' must be a numeric matrix")
    if (!identical(rownames(pwm), c("A", "C", "G", "T")))
        stop("'rownames(pwm)' must be 'c(\"A\", \"C\", \"G\", \"T\")'")
    if (!is.integer(pwm))
        storage.mode(pwm) <- "integer"
    if (any(is.na(pwm)))
        stop("'pwm' contains NAs")
    pwm
}

### Utility function for getting the max weight for each position in the PWM
maxWeights <- function(pwm)
{
    pwm <- .normalize.pwm(pwm)
    #sapply(seq_len(ncol(pwm)), function(i) max(pwm[ , i]))
    ## This will be faster than the above on large matrices
    do.call("pmax", lapply(seq_len(nrow(pwm)), function(i) pwm[i, ]))
}

### Utility function for getting the highest possible score
maxScore <- function(pwm)
{
    sum(maxWeights(pwm))
}

.normalize.min.score <- function(min.score, pwm)
{
    if (!(is.numeric(min.score) || is.character(min.score))
        || length(min.score) != 1 || is.na(min.score))
        stop("'min.score' must be a single integer or string")
    if (is.numeric(min.score)) {
        if (!is.integer(min.score))
            min.score <- as.integer(min.score)
        return(min.score)
    }
    nc <- nchar(min.score)
    if (substr(min.score, nc, nc) == "%")
        min.score <- substr(min.score, 1, nc-1)
    as.integer(maxScore(pwm) * as.numeric(min.score) / 100)
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

PWMscore <- function(pwm, subject, start=1)
{
    ## checking 'pwm'
    pwm <- .normalize.pwm(pwm)
    ## checking 'subject'
    if (!is(subject, "DNAString"))
        stop("'subject' must be a DNAString object")
    ## checking 'start'
    if (!is.numeric(start))
        stop("'start' must be a vector of integers")
    if (!is.integer(start))
        start <- as.integer(start)
    .Call("PWM_score", pwm, subject, start, PACKAGE="Biostrings")
}

### pwm: the Position Weight Matrix (integer matrix with row names A, C, G and T)
### subject: a DNAString object containing the subject sequence
### min.score: given as a percentage (e.g. "90%") of the highest possible
###            score or as an integer
.matchPWM <- function(pwm, subject, min.score, count.only=FALSE)
{
    ## checking 'pwm'
    pwm <- .normalize.pwm(pwm)
    ## checking 'subject'
    if (!is(subject, "DNAString"))
        stop("'subject' must be a DNAString object")
    ## checking 'min.score'
    min.score <- .normalize.min.score(min.score, pwm)
    ## no need to check 'count.only' (not a user controlled argument)
    ans_start <- .Call("match_PWM",
                       pwm, subject, min.score, count.only,
                       PACKAGE="Biostrings")
    if (count.only)
        return(ans_start)
    ans_width <- rep.int(ncol(pwm), length(ans_start))
    new("XStringViews", subject=subject, start=ans_start, width=ans_width, check=FALSE)
}

### 
matchPWM <- function(pwm, subject, min.score="80%")
{
    .matchPWM(pwm, subject, min.score)
}

countPWM <- function(pwm, subject, min.score="80%")
{
    .matchPWM(pwm, subject, min.score, count.only=TRUE)
}

