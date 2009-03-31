### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "mismatch", "nmatch" and "nmismatch" methods.
###

setMethod("mismatch", c(pattern = "AlignedXStringSet0", x = "missing"),
    function(pattern, x, fixed)
        pattern@mismatch
)

setMethod("nmatch", c(pattern = "PairwiseAlignedXStringSet", x = "missing"),
    function(pattern, x, fixed)
        .Call("PairwiseAlignedXStringSet_nmatch", nchar(pattern), nmismatch(pattern),
              nindel(subject(pattern))[,"WidthSum"], nindel(pattern(pattern))[,"WidthSum"],
              PACKAGE="Biostrings")
)

setMethod("nmatch", c(pattern = "PairwiseAlignedFixedSubjectSummary", x = "missing"),
    function(pattern, x, fixed)
        pattern@nmatch
)

setMethod("nmismatch", c(pattern = "AlignedXStringSet0", x = "missing"),
    function(pattern, x, fixed)
        elementLengths(mismatch(pattern))
)

setMethod("nmismatch", c(pattern = "PairwiseAlignedXStringSet", x = "missing"),
    function(pattern, x, fixed)
        nmismatch(pattern(pattern))
)

setMethod("nmismatch", c(pattern = "PairwiseAlignedFixedSubjectSummary", x = "missing"),
    function(pattern, x, fixed)
        pattern@nmismatch
)

setGeneric("nedit", function(x) standardGeneric("nedit"))

setMethod("nedit", "PairwiseAlignedXStringSet",
    function(x)
        nmismatch(x) + unname(nindel(subject(x))[,"WidthSum"]) + unname(nindel(pattern(x))[,"WidthSum"])
)

setMethod("nedit", "PairwiseAlignedFixedSubjectSummary",
    function(x)
        nmismatch(x) + unname(insertion(nindel(x))[,"WidthSum"]) + unname(deletion(nindel(x))[,"WidthSum"])
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
                  narrow(quality(unaligned(x))[output[["Id"]]], start = output[["Start"]],
                         end = output[["End"]])
            output <- cbind(output, "Quality" = as.character(quality))
        }
        if (any(nchar(prefixColNames) > 0))
            names(output) <- paste(prefixColNames, names(output), sep = "")
        output
    }
)

setMethod("mismatchTable", "PairwiseAlignedXStringSet",
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
                    alphabetFrequency(narrow(quality(unaligned(x))[nonEmptyAlignment], start = start(x)[nonEmptyAlignment],
                                             end = end(x)[nonEmptyAlignment]), collapse = TRUE)[qualityValues + 1]
            else
                qualityAll <-
                  colSums(as.numeric(weight)[nonEmptyAlignment] *
                          alphabetFrequency(narrow(quality(unaligned(x))[nonEmptyAlignment], start = start(x)[nonEmptyAlignment],
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

setMethod("mismatchSummary", "PairwiseAlignedFixedSubject",
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
                          "Probability" =
                           as.vector(subjectTable) /
                             coverage(subject(x), weight = weight)[subjectPosition, drop = TRUE]))
        output[["subject"]] <- output[["subject"]][order(output[["subject"]][[1]], output[["subject"]][[2]]),]
        rownames(output[["subject"]]) <- as.character(seq_len(nrow(output[["subject"]])))
        output
    }
)

setMethod("mismatchSummary", "PairwiseAlignedFixedSubjectSummary",
    function(x)
        x@mismatchSummary
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "coverage" methods.
###

setMethod("coverage", "AlignedXStringSet0",
    function(x, start = NA, end = NA, weight = 1L)
    {
        if (any(is.na(start)))
            start <- 1
        if (any(is.na(end)))
            end <- max(nchar(unaligned(x)))
        coverage(x@range, start = start, end = end, weight = weight)
    }
)

setMethod("coverage", "PairwiseAlignedFixedSubject",
    function(x, start = NA, end = NA, weight = 1L)
        coverage(subject(x), start = NA, end = NA, weight = 1L)
)

setMethod("coverage", "PairwiseAlignedFixedSubjectSummary",
    function(x, start = NA, end = NA)
        subseq(x@coverage, start=start, end=end)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Compare strings
###

setGeneric("compareStrings", signature = c("pattern", "subject"),
           function(pattern, subject)  standardGeneric("compareStrings"))
setMethod("compareStrings", signature = c(pattern = "character", subject = "character"),
          function(pattern, subject) {
              if (length(pattern) != length(subject))
                  stop("'pattern' and 'subject' must have the same length")
              ncharPattern <- nchar(pattern)
              if (any(ncharPattern != nchar(subject)))
                  stop("'pattern' and 'subject' must have the same number of characters")
              .Call("align_compareStrings", pattern, subject, max(ncharPattern), "+", "-", "?", PACKAGE="Biostrings")
          })
setMethod("compareStrings", signature = c(pattern = "XString", subject = "XString"),
          function(pattern, subject) {
              compareStrings(as.character(pattern), as.character(subject))
          })
setMethod("compareStrings", signature = c(pattern = "XStringSet", subject = "XStringSet"),
          function(pattern, subject) {
              compareStrings(as.character(pattern), as.character(subject))
          })
setMethod("compareStrings", signature = c(pattern = "AlignedXStringSet0", subject = "AlignedXStringSet0"),
          function(pattern, subject) {
              compareStrings(as.character(pattern), as.character(subject))
          })
setMethod("compareStrings", signature = c(pattern = "PairwiseAlignedXStringSet", subject = "missing"),
		  function(pattern, subject) {
			  compareStrings(pattern@pattern, pattern@subject)
		  })


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Alignment consensus matrix
###

setGeneric("consmat", signature="x", function(x, ...)  standardGeneric("consmat"))

setMethod("consmat", "ANY",
          function(x, ...)
          {
              .Deprecated("consensusMatrix")
              consensusMatrix(x, ...)
          })


setGeneric("consensusMatrix", signature="x", function(x, ...)  standardGeneric("consensusMatrix"))

setMethod("consensusMatrix", "character",
          function(x, freq=FALSE)
          {
              consensusMatrix(BStringSet(x), freq=freq)
          })

setMethod("consensusMatrix", "matrix",
          function(x, freq=FALSE)
          {
              consensusMatrix(BStringSet(apply(x, 1, paste, collapse="")), freq=freq)
          })

### 'x' must be a list of FASTA records as one returned by readFASTA()
setMethod("consensusMatrix", "list",
          function(x, freq=FALSE)
          {
              consensusMatrix(BStringSet(FASTArecordsToCharacter(x, use.names=FALSE)), freq=freq)
          })

setMethod("consensusMatrix", "XStringSet",
          function(x, baseOnly=FALSE, freq=FALSE)
          {
              codes <- xscodes(x, baseOnly=baseOnly)
              if (is.null(names(codes))) {
                  names(codes) <- intToUtf8(codes, multiple = TRUE)
                  removeUnused <- TRUE
              } else {
                  removeUnused <- FALSE
              }
              freq <- .normargFreq(freq)
              ans <- .Call("XStringSet_char_frequency_by_pos",
                      x, codes, baseOnly, PACKAGE="Biostrings")
              if (removeUnused) {
                  ans <- ans[rowSums(ans) > 0, , drop=FALSE]
              }
              if (freq) {
                  ans <- t(t(ans) / colSums(ans))
              }
              ans
          })

setMethod("consensusMatrix", "XStringViews",
          function(x, baseOnly=FALSE, freq=FALSE)
          {
              y <- XStringViewsToSet(x, use.names=FALSE, verbose=FALSE)
              consensusMatrix(y, baseOnly=baseOnly, freq=freq)
          })

setMethod("consensusMatrix", "PairwiseAlignedFixedSubject",
          function(x, baseOnly=FALSE, freq=FALSE, gapCode="-", endgapCode="-")
          {
              consensusMatrix(aligned(x, gapCode=gapCode, endgapCode=endgapCode),
                              baseOnly=baseOnly, freq=freq)
          })

setGeneric("consensusString", signature = "x", function(x)  standardGeneric("consensusString"))

setMethod("consensusString", "ANY",
          function(x)
          {
              mat <- consensusMatrix(x, freq=TRUE)
              letters <- rownames(mat)
              paste(apply(mat, 2, function(y) if(any(y > 0.5)) letters[which.max(y)] else "?"), collapse="")
          })
