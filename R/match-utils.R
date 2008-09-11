### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helper functions (not exported) used by matching functions from other
### files (like matchPattern(), matchPDict(), etc...) to check and normalize
### their arguments.
###

normargMaxMismatch <- function(max.mismatch)
{
    if (!isSingleNumber(max.mismatch))
        stop("'max.mismatch' must be a single integer")
    max.mismatch <- as.integer(max.mismatch)
    if (max.mismatch < 0)
        stop("'max.mismatch' must be a non-negative integer")
    max.mismatch
}

### Return a logical vector of length 2.
normargFixed <- function(fixed, subjectClass)
{
    if (!is.logical(fixed) && !is.character(fixed))
        stop("'fixed' not a logical or character vector")
    if (is.logical(fixed)) {
        if (any(is.na(fixed)))
            stop("'fixed' has NAs")
        fixed_names <- names(fixed)
        if (is.null(fixed_names)) {
            if (!(length(fixed) %in% 1:2))
                stop("when an unnamed logical vector, ",
                     "'fixed' fixed must be of length 1 or 2")
            if (length(fixed) == 1)
                fixed <- c(fixed, fixed)
        } else {
            if (length(fixed) != 2)
                stop("when a named logical vector, 'fixed' must be of length 2")
            if (!setequal(fixed_names, c("pattern", "subject")))
                stop("'fixed' names must be \"pattern\" and \"subject\"")
            fixed <- c(fixed["pattern"], fixed["subject"])
        }
    } else if (is.character(fixed)) {
        if (any(duplicated(fixed)) || !all(fixed %in% c("pattern", "subject")))
            stop("when a character vector, 'fixed' must be ",
                 "a subset of 'c(\"pattern\", \"subject\")' ",
                 "with no duplicated")
        fixed <- c("pattern" %in% fixed, "subject" %in% fixed)
    }
    if (!all(fixed) && !extends(subjectClass, "DNAString") && !extends(subjectClass, "RNAString"))
        stop("'fixed' value only supported for a DNAString or RNAString subject ",
             "(you can only use 'fixed=TRUE' with your subject)")
    fixed
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "nmismatchStartingAt" and "nmismatchEndingAt" generic and methods.
###
### 'starting.at' (or 'ending.at') must be integer vectors containing the
### starting (or ending) positions of the pattern relatively to the subject.
### The two functions return an integer vector of the same length as
### 'starting.at' (or 'ending.at').
###

setGeneric("nmismatchStartingAt", signature="subject",
    function(pattern, subject, starting.at=1, fixed=TRUE)
        standardGeneric("nmismatchStartingAt")
)

setGeneric("nmismatchEndingAt", signature="subject",
    function(pattern, subject, ending.at=1, fixed=TRUE)
        standardGeneric("nmismatchEndingAt")
)

### If 'starting=TRUE' then 'at' contains starting positions, otherwise it
### contains ending positions.
.nmismatchAt <- function(pattern, subject, starting, at, fixed)
{
    if (!is(subject, "XString"))
        subject <- XString(NULL, subject)
    if (class(pattern) != class(subject))
        pattern <- XString(class(subject), pattern)
    if (!is.numeric(at)) {
        what <- if (starting) "starting.at" else "ending.at"
        stop("'", what, "'  must be a vector of integers")
    }
    if (!is.integer(at))
        at <- as.integer(at)
    fixed <- normargFixed(fixed, class(subject))
    .Call("nmismatch_at", pattern, subject,
          starting, at, fixed,
          PACKAGE="Biostrings")
}

### Dispatch on 'subject' (see signature of generic).
setMethod("nmismatchStartingAt", "character",
    function(pattern, subject, starting.at=1, fixed=TRUE)
        .nmismatchAt(pattern, subject, TRUE, starting.at, fixed)
)
setMethod("nmismatchEndingAt", "character",
    function(pattern, subject, ending.at=1, fixed=TRUE)
        .nmismatchAt(pattern, subject, FALSE, ending.at, fixed)
)

### Dispatch on 'subject' (see signature of generic).
setMethod("nmismatchStartingAt", "XString",
    function(pattern, subject, starting.at=1, fixed=TRUE)
        .nmismatchAt(pattern, subject, TRUE, starting.at, fixed)
)
setMethod("nmismatchEndingAt", "XString",
    function(pattern, subject, ending.at=1, fixed=TRUE)
        .nmismatchAt(pattern, subject, FALSE, ending.at, fixed)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "isMatching" generic and methods.
###
### Return a logical vector of the same length as 'start'.
###

setGeneric("isMatching", signature="subject",
    function(pattern, subject, start=1, max.mismatch=0, fixed=TRUE)
        standardGeneric("isMatching")
)

.isMatching <- function(pattern, subject, start, max.mismatch, fixed)
{
    if (!is(subject, "XString"))
        subject <- XString(NULL, subject)
    if (class(pattern) != class(subject))
        pattern <- XString(class(subject), pattern)
    if (!is.numeric(start))
        stop("'start' must be a vector of integers")
    if (!is.integer(start))
        start <- as.integer(start)
    max.mismatch <- normargMaxMismatch(max.mismatch)
    fixed <- normargFixed(fixed, class(subject))
    .Call("is_matching", pattern, subject, start,
          max.mismatch, fixed,
          PACKAGE="Biostrings")
}

### Dispatch on 'subject' (see signature of generic).
setMethod("isMatching", "character",
    function(pattern, subject, start=1, max.mismatch=0, fixed=TRUE)
        .isMatching(pattern, subject, start, max.mismatch, fixed)
)

### Dispatch on 'subject' (see signature of generic).
setMethod("isMatching", "XString",
    function(pattern, subject, start=1, max.mismatch=0, fixed=TRUE)
        .isMatching(pattern, subject, start, max.mismatch, fixed)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### mismatch()
###

### Helper function used by .mismatch()
### Returns a vector of the positions of mismatches of 'pattern'
### in a view on 'subject' starting at 'start' and whose width is length(pattern).
.bsMismatch <- function(pattern, subject, start, fixed)
{
    mm <- integer(0)
    j0 <- start - as.integer(1)
    for (i in seq_len(length(pattern))) {
        j <- j0 + i
        if (j < 1 || j > length(subject)) {
            mm <- c(mm, i)
        } else {
            l <- XString.substr(pattern, i, i)
            cp <- isMatching(l, subject, j, max.mismatch=0, fixed=fixed)
            if (cp == 0)
                mm <- c(mm, i)
        }
    }
    mm
}

.mismatch <- function(pattern, x, fixed)
{
    if (length(x) == 0)
        return(list())
    if (any(width(x) != length(pattern)))
        warning("views in 'x' don't have a width equal to pattern length")
    lapply(1:length(x),
           function(i) .bsMismatch(pattern, subject(x), start(x)[i], fixed))
}

setGeneric("mismatch", signature=c("pattern", "x"),
    function(pattern, x, fixed=TRUE) standardGeneric("mismatch")
)

### Typical use:
###   mp <- matchPattern("TGA", DNAString("GTGACGTGCAT"), max.mismatch=2)
###   mismatch("TGA", mp)
### Dispatch on 'x' (see signature of generic).
setMethod("mismatch", c(pattern="ANY", x="XStringViews"),
    function(pattern, x, fixed)
    {
        if (class(pattern) != class(x@subject))
            pattern <- XString(class(x@subject), pattern)
        .mismatch(pattern, x, fixed)
    }
)

setMethod("mismatch", c(pattern = "AlignedXStringSet", x = "missing"),
    function(pattern, x, fixed)
        pattern@mismatch
)

setGeneric("nmatch", signature=c("pattern", "x"),
    function(pattern, x, fixed=TRUE) standardGeneric("nmatch")
)

setMethod("nmatch", c(pattern="ANY", x="XStringViews"),
    function(pattern, x, fixed)
    {
        funCall <- match.call()
        funCall[[1]] <- as.name("nmismatch")
        nchar(pattern) - eval(funCall, sys.parent())
    }
)

setMethod("nmatch", c(pattern = "PairwiseAlignment", x = "missing"),
    function(pattern, x, fixed)
        .Call("nmatch_PairwiseAlignment", nchar(pattern), nmismatch(pattern),
              nindel(subject(pattern))[,"WidthSum"], nindel(pattern(pattern))[,"WidthSum"],
              PACKAGE="Biostrings")
)

setMethod("nmatch", c(pattern = "PairwiseAlignmentSummary", x = "missing"),
    function(pattern, x, fixed)
        pattern@nmatch
)

setGeneric("nmismatch", signature=c("pattern", "x"),
    function(pattern, x, fixed=TRUE) standardGeneric("nmismatch")
)

setMethod("nmismatch", c(pattern="ANY", x="XStringViews"),
    function(pattern, x, fixed)
    {
        funCall <- match.call()
        funCall[[1]] <- as.name("mismatch")
        mismatches <- eval(funCall, sys.parent())
        sapplyLength(mismatches)
    }
)

setMethod("nmismatch", c(pattern = "AlignedXStringSet", x = "missing"),
    function(pattern, x, fixed)
    {
        mismatches <- mismatch(pattern)
        sapplyLength(mismatches)
    }
)

setMethod("nmismatch", c(pattern = "PairwiseAlignment", x = "missing"),
    function(pattern, x, fixed)
        nmismatch(pattern(pattern))
)

setMethod("nmismatch", c(pattern = "PairwiseAlignmentSummary", x = "missing"),
    function(pattern, x, fixed)
        pattern@nmismatch
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "mismatchTable" generic and methods.
###

setGeneric("mismatchTable", signature = "x",
    function(x, shiftLeft = 0L, shiftRight = 0L, ...)
    standardGeneric("mismatchTable")
)

setMethod("mismatchTable", "AlignedXStringSet",
    function(x, shiftLeft = 0L, shiftRight = 0L, prefixColNames = "")
    {
        if (!isSingleNumber(shiftLeft) || shiftLeft > 0)
            stop("'shiftLeft' must be a non-positive integer")
        if (!isSingleNumber(shiftRight) || shiftRight < 0)
            stop("'shiftRight' must be a non-negative integer")
        shiftLeft <- as.integer(shiftLeft)
        shiftRight <- as.integer(shiftRight)
        nMismatch <- nmismatch(x)
        id <- rep.int(seq_len(length(nMismatch)), nMismatch)
        singleString <- (length(unaligned(x)) == 1)
        if (singleString)
            subset <- unaligned(x)[[1]]
        else
            subset <- unaligned(x)[id]
        position <- unlist(mismatch(x))
        if (shiftLeft == 0L)
            start <- position
        else
            start <- pmax(position + shiftLeft, 1L)
        if (shiftRight == 0L)
            end <- position
        else
            end <- pmin(position + shiftRight, width(subset))
        if (singleString)
            substring <- Views(subset, start=start, end=end)
        else
            substring <- narrow(subset, start=start, end=end)
        output <-
          data.frame("Id" = id,
                     "Start" = start,
                     "End" = end,
                     "Substring" = as.character(substring))
        if (any(nchar(prefixColNames) > 0))
            names(output) <- paste(prefixColNames, names(output), sep = "")
        output
    }
)

setMethod("mismatchTable", "QualityAlignedXStringSet",
    function(x, shiftLeft = 0L, shiftRight = 0L, prefixColNames = "")
    {
        output <-
          callNextMethod(x, shiftLeft = shiftLeft, shiftRight = shiftRight,
                         prefixColNames = "")
        if (nrow(output) == 0) {
            output <- cbind(output, "Quality" = character(0))
		} else {
            if (length(quality(x)) == 1) {
                if (width(quality(x)) == 1)
                    quality <-
                      Views(quality(x)[[1]][rep.int(1L, max(output[["End"]]))],
                            start = output[["Start"]], end = output[["End"]])
                else
                    quality <-
                      Views(quality(x)[[1]], start = output[["Start"]],
                            end = output[["End"]])
            }
            else
                quality <-
                  narrow(quality(x)[output[["Id"]]], start = output[["Start"]],
                         end = output[["End"]])
            output <- cbind(output, "Quality" = as.character(quality))
        }
        if (any(nchar(prefixColNames) > 0))
            names(output) <- paste(prefixColNames, names(output), sep = "")
        output
    }
)

setMethod("mismatchTable", "PairwiseAlignment",
    function(x, shiftLeft = 0L, shiftRight = 0L)
    {
        cbind(mismatchTable(pattern(x), shiftLeft = shiftLeft,
                            shiftRight = shiftRight, prefixColNames = "Pattern"),
              mismatchTable(subject(x), shiftLeft = shiftLeft,
                            shiftRight = shiftRight, prefixColNames = "Subject")[,-1])
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "mismatchSummary" generic and methods.
###

setGeneric("mismatchSummary", signature = "x",
    function(x, ...) standardGeneric("mismatchSummary")
)

setMethod("mismatchSummary", "AlignedXStringSet",
    function(x, weight=1L, .mismatchTable=mismatchTable(x))
    {
        if (!is.numeric(weight) || !(length(weight) %in% c(1, length(x))))
            stop("'weight' must be an integer vector with length 1 or 'length(x)'")
        if (!is.integer(weight))
            weight <- as.integer(weight)
        coverageTable <- as.integer(coverage(x, weight = weight))
        n <- length(coverageTable)
        if (length(weight) == 1)
            endTable <- weight * table(.mismatchTable[["End"]])
        else
            endTable <- table(rep(.mismatchTable[["End"]], weight[.mismatchTable[["Id"]]]))
        countTable <- rep(0L, n)
        countTable[as.integer(names(endTable))] <- endTable
        list("position" =
             data.frame("Position" = seq_len(n),
                        "Count" = countTable,
                        "Probability" = countTable / coverageTable))
    }
)

setMethod("mismatchSummary", "QualityAlignedXStringSet",
    function(x, weight=1L, .mismatchTable=mismatchTable(x))
    {
        if (!is.numeric(weight) || !(length(weight) %in% c(1, length(x))))
            stop("'weight' must be an integer vector with length 1 or 'length(x)'")
        if (!is.integer(weight))
            weight <- as.integer(weight)
        if (is(quality(x), "PhredQuality")) {
            qualityValues <- 33:(33 + 99)
            qualityZero <- 33
        } else if (is(quality(x), "SolexaQuality")) {
            qualityValues <- 59:(59 + 104)
            qualityZero <- 64
        } else
            stop("unrecognized quality class")
        if ((length(quality(x)) == 1) && (nchar(quality(x)) == 1))
            qualityAll <-
              sum(as.numeric(weight) * width(x)) *
                alphabetFrequency(quality(x), collapse = TRUE)[qualityValues + 1]
        else {
            nonEmptyAlignment <- (width(x) > 0)
            if (length(weight) == 1)
                qualityAll <-
                  as.numeric(weight) *
                    alphabetFrequency(narrow(quality(x)[nonEmptyAlignment], start = start(x)[nonEmptyAlignment],
                                             end = end(x)[nonEmptyAlignment]), collapse = TRUE)[qualityValues + 1]
            else
                qualityAll <-
                  colSums(as.numeric(weight)[nonEmptyAlignment] *
                          alphabetFrequency(narrow(quality(x)[nonEmptyAlignment], start = start(x)[nonEmptyAlignment],
                                                   end = end(x)[nonEmptyAlignment]))[, qualityValues + 1, drop=FALSE])
        }
        names(qualityAll) <- sapply(as.raw(qualityValues), rawToChar)
        qualityAll <- qualityAll[qualityAll > 0]
        if (length(weight) == 1)
            qualityTable <- weight * table(.mismatchTable[["Quality"]])
        else
            qualityTable <-
              table(rep(.mismatchTable[["Quality"]], weight[.mismatchTable[["Id"]]]))
        qualityCounts <- rep(0L, length(qualityAll))
        names(qualityCounts) <- names(qualityAll)
        qualityCounts[names(qualityTable)] <- qualityTable
        c(callNextMethod(x, weight = weight, .mismatchTable = .mismatchTable),
          list("quality" =
               data.frame("Quality" = unlist(lapply(names(qualityAll), utf8ToInt)) - qualityZero,
                          "Count" = qualityCounts,
                          "Probability" = qualityCounts / qualityAll)))
    }
)

setMethod("mismatchSummary", "PairwiseAlignment",
    function(x, weight=1L)
    {
        if (!is.numeric(weight) || !(length(weight) %in% c(1, length(x))))
            stop("'weight' must be an integer vector with length 1 or 'length(x)'")
        if (!is.integer(weight))
            weight <- as.integer(weight)
        mismatchTable <- list("pattern" = mismatchTable(pattern(x)), "subject" = mismatchTable(subject(x)))
        combinedInfo <-
          paste(mismatchTable[["subject"]][["End"]], mismatchTable[["pattern"]][["Substring"]], sep = "\001")
        if (length(weight) == 1)
            subjectTable <- weight * table(combinedInfo)
        else
            subjectTable <- table(rep(combinedInfo, weight[mismatchTable[["pattern"]][["Id"]]]))
        if (length(subjectTable) == 0) {
            subjectTableLabels <- character(0)
            subjectPosition <- integer(0)
        } else {
            subjectTableLabels <- strsplit(names(subjectTable), split = "\001")
            subjectPosition <- as.integer(unlist(lapply(subjectTableLabels, "[", 1)))
        }
        output <-
          list("pattern" = mismatchSummary(pattern(x), weight = weight, .mismatchTable = mismatchTable[["pattern"]]),
               "subject" =
               data.frame("SubjectPosition" = subjectPosition,
                          "Subject" = safeExplode(letter(unaligned(subject(x))[[1]], subjectPosition)),
                          "Pattern" = unlist(lapply(subjectTableLabels, "[", 2)),
                          "Count" = as.vector(subjectTable),
                          "Probability" = as.vector(subjectTable) / coverage(subject(x), weight = weight)[subjectPosition]))
        output[["subject"]] <- output[["subject"]][order(output[["subject"]][[1]], output[["subject"]][[2]]),]
        rownames(output[["subject"]]) <- as.character(seq_len(nrow(output[["subject"]])))
        output
    }
)

setMethod("mismatchSummary", "PairwiseAlignmentSummary",
    function(x)
        x@mismatchSummary
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "coverage" methods.
###

setMethod("coverage", "XStringViews",
    function(x, start=NA, end=NA, weight=1L)
    {
        if (!isSingleNumberOrNA(start))
            stop("'start' must be a single integer or NA")
        if (!is.integer(start))
            start <- as.integer(start)
        if (is.na(start))
            start <- 1L
        if (!isSingleNumberOrNA(end))
            stop("'end' must be a single integer or NA")
        if (!is.integer(end))
            end <- as.integer(end)
        if (is.na(end))
            end <- length(subject(x))
        if (!is.numeric(weight) || !(length(weight) %in% c(1, length(x))))
            stop("'weight' must be an integer vector with length 1 or 'length(x)'")
        if (!is.integer(weight))
            weight <- as.integer(weight)
        callNextMethod(x, start=start, end=end, weight=weight)
    }
)

setMethod("coverage", "MaskedXString",
    function(x, start=NA, end=NA, weight=1L)
        coverage(masks(x), start=start, end=end, weight=weight)
)

setMethod("coverage", "MIndex",
    function(x, start=NA, end=NA)
    {
        if (!isSingleNumber(start))
            stop("'start' must be a single integer")
        if (!is.integer(start))
            start <- as.integer(start)
        if (!isSingleNumber(end))
            stop("'end' must be a single integer")
        if (!is.integer(end))
            end <- as.integer(end)
        width <- end - start + 1L
        if (width < 0)
            stop("'end' must be >= 'start' - 1")
		ans <- XInteger(end - start + 1L, initialize = TRUE)
        if (is(x, "ByPos_MIndex"))
            .Call("ByPos_MIndex_coverage",
                  endIndex(x), x@width, start, ans@xdata@xp,
                  PACKAGE="Biostrings")
        else if (is(x, "ByName_MIndex"))
            .Call("ByName_MIndex_coverage",
                  x@ends_envir, x@width, start, ans@xdata@xp,
                  PACKAGE="Biostrings")
        else
            stop("Biostrings internal error: unknown MIndex subtype ", class(x))
        ans
    }
)

setMethod("coverage", "AlignedXStringSet",
    function(x, start = NA, end = NA, weight = 1L)
    {
        if (any(is.na(start)))
            start <- 1
        if (any(is.na(end)))
            end <- max(nchar(unaligned(x)))
        coverage(x@range, start = start, end = end, weight = weight)
    }
)

setMethod("coverage", "PairwiseAlignment",
    function(x, start = NA, end = NA, weight = 1L)
        coverage(subject(x), start = NA, end = NA, weight = 1L)
)

setMethod("coverage", "PairwiseAlignmentSummary",
    function(x, start = NA, end = NA)
    {
        if (!isSingleNumberOrNA(start))
            stop("'start' must be a single integer or NA")
        if (!is.integer(start))
            start <- as.integer(start)
        if (is.na(start))
            start <- 1L
        if (!isSingleNumberOrNA(end))
            stop("'end' must be a single integer or NA")
        if (!is.integer(end))
            end <- as.integer(end)
        if (is.na(end))
            end <- length(x@coverage)
        if (start > 1 || end < length(x@coverage))
            ans <- XInteger(x@coverage[start:end])
        else
            ans <- x@coverage
        ans
	}
)
