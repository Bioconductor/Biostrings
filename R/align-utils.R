### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "mismatch", "nmatch" and "nmismatch" methods.
###

setMethod("mismatch", c(pattern = "AlignedXStringSet0", x = "missing"),
    function(pattern, x, fixed)
        pattern@mismatch
)

setMethod("nmatch", c(pattern = "PairwiseAlignments", x = "missing"),
    function(pattern, x, fixed)
        .Call2("PairwiseAlignments_nmatch", nchar(pattern), nmismatch(pattern),
              nindel(subject(pattern))[,"WidthSum"], nindel(pattern(pattern))[,"WidthSum"],
              PACKAGE="Biostrings")
)

setMethod("nmatch", c(pattern = "PairwiseAlignmentsSingleSubjectSummary", x = "missing"),
    function(pattern, x, fixed)
        pattern@nmatch
)

setMethod("nmismatch", c(pattern = "AlignedXStringSet0", x = "missing"),
    function(pattern, x, fixed)
        elementNROWS(mismatch(pattern))
)

setMethod("nmismatch", c(pattern = "PairwiseAlignments", x = "missing"),
    function(pattern, x, fixed)
        nmismatch(pattern(pattern))
)

setMethod("nmismatch", c(pattern = "PairwiseAlignmentsSingleSubjectSummary", x = "missing"),
    function(pattern, x, fixed)
        pattern@nmismatch
)

setGeneric("nedit", function(x) standardGeneric("nedit"))

setMethod("nedit", "PairwiseAlignments",
    function(x)
        nmismatch(x) + unname(nindel(subject(x))[,"WidthSum"]) +
          unname(nindel(pattern(x))[,"WidthSum"])
)

setMethod("nedit", "PairwiseAlignmentsSingleSubjectSummary",
    function(x)
        nmismatch(x) + unname(insertion(nindel(x))[,"WidthSum"]) +
          unname(deletion(nindel(x))[,"WidthSum"])
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "mismatchTable" generic and methods.
###

setGeneric("mismatchTable", signature = "x",
    function(x, shiftLeft = 0L, shiftRight = 0L, ...)
    standardGeneric("mismatchTable")
)

setMethod("mismatchTable", "AlignedXStringSet0",
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
            end <- pmin(position + shiftRight, nchar(subset))
        if (length(subset) == 0)
            substring <- character(0)
        else if (singleString)
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
            if (length(quality(unaligned(x))) == 1) {
                if (width(quality(unaligned(x))) == 1)
                    quality <-
                      Views(quality(unaligned(x))[[1]][rep.int(1L, max(output[["End"]]))],
                            start = output[["Start"]], end = output[["End"]])
                else
                    quality <-
                      Views(quality(unaligned(x))[[1]], start = output[["Start"]],
                            end = output[["End"]])
            }
            else
                quality <-
                  narrow(quality(unaligned(x))[output[["Id"]]],
                         start = output[["Start"]], end = output[["End"]])
            output <- cbind(output, "Quality" = as.character(quality))
        }
        if (any(nchar(prefixColNames) > 0))
            names(output) <- paste(prefixColNames, names(output), sep = "")
        output
    }
)

setMethod("mismatchTable", "PairwiseAlignments",
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

setMethod("mismatchSummary", "AlignedXStringSet0",
    function(x, weight=1L, .mismatchTable=mismatchTable(x))
    {
        if (!is.numeric(weight) || !(length(weight) %in% c(1, length(x))))
            stop("'weight' must be an integer vector with length 1 or 'length(x)'")
        if (!is.integer(weight))
            weight <- as.integer(weight)
        coverageTable <- as.vector(coverage(x, weight = weight))
        n <- length(coverageTable)
        if (length(weight) == 1)
            endTable <- weight * table(.mismatchTable[["End"]])
        else
            endTable <-
              table(rep(.mismatchTable[["End"]], weight[.mismatchTable[["Id"]]]))
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
        qualityValues <-
          (minQuality(quality(unaligned(x))) + offset(quality(unaligned(x)))):
          (maxQuality(quality(unaligned(x))) + offset(quality(unaligned(x))))
        qualityZero <- offset(quality(unaligned(x)))
        if ((length(quality(unaligned(x))) == 1) && (nchar(quality(unaligned(x))) == 1))
            qualityAll <-
              sum(as.numeric(weight) * width(x)) *
                alphabetFrequency(quality(unaligned(x)), collapse = TRUE)[qualityValues + 1]
        else {
            nonEmptyAlignment <- (width(x) > 0)
            if (length(weight) == 1)
                qualityAll <-
                  as.numeric(weight) *
                    alphabetFrequency(narrow(quality(unaligned(x))[nonEmptyAlignment],
                                             start = start(x)[nonEmptyAlignment],
                                             end = end(x)[nonEmptyAlignment]),
                                      collapse = TRUE)[qualityValues + 1]
            else
                qualityAll <-
                  colSums(as.numeric(weight)[nonEmptyAlignment] *
                          alphabetFrequency(narrow(quality(unaligned(x))[nonEmptyAlignment],
                                                   start = start(x)[nonEmptyAlignment],
                                                   end = end(x)[nonEmptyAlignment])
                                            )[, qualityValues + 1, drop=FALSE])
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
               data.frame("Quality" =
                          unlist(lapply(names(qualityAll), utf8ToInt)) - qualityZero,
                          "Count" = qualityCounts,
                          "Probability" = qualityCounts / qualityAll)))
    }
)

setMethod("mismatchSummary", "PairwiseAlignmentsSingleSubject",
    function(x, weight=1L)
    {
        if (!is.numeric(weight) || !(length(weight) %in% c(1, length(x))))
            stop("'weight' must be an integer vector with length 1 or 'length(x)'")
        if (!is.integer(weight))
            weight <- as.integer(weight)
        mismatchTable <-
          list("pattern" = mismatchTable(pattern(x)),
               "subject" = mismatchTable(subject(x)))
        combinedInfo <-
          paste(mismatchTable[["subject"]][["End"]],
                mismatchTable[["pattern"]][["Substring"]], sep = "\001")
        if (length(weight) == 1)
            subjectTable <- weight * table(combinedInfo)
        else
            subjectTable <-
              table(rep(combinedInfo, weight[mismatchTable[["pattern"]][["Id"]]]))
        if (length(subjectTable) == 0) {
            subjectTableLabels <- character(0)
            subjectPosition <- integer(0)
        } else {
            subjectTableLabels <- strsplit(names(subjectTable), split = "\001")
            subjectPosition <- as.integer(unlist(lapply(subjectTableLabels, "[", 1)))
        }
        output <-
          list("pattern" =
               mismatchSummary(pattern(x), weight = weight,
                               .mismatchTable = mismatchTable[["pattern"]]),
               "subject" =
               data.frame("SubjectPosition" = subjectPosition,
                          "Subject" =
                          safeExplode(letter(unaligned(subject(x))[[1]], subjectPosition)),
                          "Pattern" = unlist(lapply(subjectTableLabels, "[", 2)),
                          "Count" = as.vector(subjectTable),
                          "Probability" =
                           as.vector(subjectTable) /
                             coverage(subject(x), weight = weight)[subjectPosition, drop = TRUE]))
        output[["subject"]] <-
          output[["subject"]][order(output[["subject"]][[1]], output[["subject"]][[2]]),]
        rownames(output[["subject"]]) <- as.character(seq_len(nrow(output[["subject"]])))
        output
    }
)

setMethod("mismatchSummary", "PairwiseAlignmentsSingleSubjectSummary",
    function(x)
        x@mismatchSummary
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "coverage" methods.
###

setMethod("coverage", "AlignedXStringSet0",
    function(x, shift=0L, width=NULL, weight=1L)
    {
        shift <- recycleIntegerArg(shift, "shift", length(x@range))
        if (is.null(width))
            width <- max(nchar(unaligned(x))) + max(shift)
        coverage(x@range, shift=shift, width=width, weight=weight)
    }
)

setMethod("coverage", "PairwiseAlignmentsSingleSubject",
    function(x, shift=0L, width=NULL, weight=1L)
        coverage(subject(x), shift=shift, width=width, weight=weight)
)

setMethod("coverage", "PairwiseAlignmentsSingleSubjectSummary",
    function(x, shift=0L, width=NULL, weight=1L) {
        if (shift != 0L)
            stop("'shift' argument is not supported for 'PairwiseAlignmentsSingleSubjectSummary' objects")
        if (weight != 1L)
            stop("'weight' argument is not supported for 'PairwiseAlignmentsSingleSubjectSummary' objects")
        window(x@coverage, width=width)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Compare strings
###

setGeneric("compareStrings", signature = c("pattern", "subject"),
           function(pattern, subject)  standardGeneric("compareStrings"))
setMethod("compareStrings",
          signature = c(pattern = "character", subject = "character"),
          function(pattern, subject) {
              if (length(pattern) != length(subject))
                  stop("'pattern' and 'subject' must have the same length")
              ncharPattern <- nchar(pattern)
              if (any(ncharPattern != nchar(subject)))
                  stop("'pattern' and 'subject' must have the same number of characters")
              .Call2("align_compareStrings", pattern, subject, max(ncharPattern), "+", "-", "?",
                    PACKAGE="Biostrings")
          })
setMethod("compareStrings",
          signature = c(pattern = "XString", subject = "XString"),
          function(pattern, subject) {
              compareStrings(as.character(pattern), as.character(subject))
          })
setMethod("compareStrings",
          signature = c(pattern = "XStringSet", subject = "XStringSet"),
          function(pattern, subject) {
              compareStrings(as.character(pattern), as.character(subject))
          })
setMethod("compareStrings",
          signature = c(pattern = "AlignedXStringSet0", subject = "AlignedXStringSet0"),
          function(pattern, subject) {
              compareStrings(as.character(pattern), as.character(subject))
          })
setMethod("compareStrings",
          signature = c(pattern = "PairwiseAlignments", subject = "missing"),
          function(pattern, subject) {
              compareStrings(pattern@pattern, pattern@subject)
          })


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "consensusMatrix" method for PairwiseAlignmentsSingleSubject objects.
###

setMethod("consensusMatrix", "PairwiseAlignmentsSingleSubject",
          function(x, as.prob=FALSE, shift=0L, width=NULL,
                   baseOnly=FALSE, gapCode="-", endgapCode="-")
          {
              if (!identical(shift, 0L) || !identical(width, NULL))
                  stop("\"consensusMatrix\" method for PairwiseAlignmentsSingleSubject objects ",
                       "doesn't support the 'shift' and 'width' arguments")
              consensusMatrix(aligned(x, gapCode=gapCode, endgapCode=endgapCode),
                              as.prob=as.prob, baseOnly=baseOnly)
          })

