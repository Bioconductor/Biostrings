### ==========================================================================
### PairwiseAlignmentsSingleSubject objects
### --------------------------------------------------------------------------
### A PairwiseAlignmentsSingleSubject object contains the result of the
### pairwise alignment of many patterns to one subject.
###

setClass("PairwiseAlignmentsSingleSubject",
    contains="PairwiseAlignments"
)

setClass("PairwiseAlignmentsSingleSubjectSummary",
    representation(
        type="character",
        score="numeric",
        nmatch="numeric",
        nmismatch="numeric",
        ninsertion="matrix",
        ndeletion="matrix",
        coverage="Rle",
        mismatchSummary="list"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructors.
###

setGeneric("PairwiseAlignmentsSingleSubject",
    function(pattern, subject, type = "global", substitutionMatrix = NULL,
             gapOpening = 0, gapExtension = 1, ...)
    standardGeneric("PairwiseAlignmentsSingleSubject"))
setMethod("PairwiseAlignmentsSingleSubject", signature(pattern = "XString", subject = "XString"),
    function(pattern, subject, type = "global", substitutionMatrix = NULL,
             gapOpening = 0, gapExtension = 1) {
        if (seqtype(pattern) != seqtype(subject))
            stop("'pattern' and 'subject' must contain ",
                 "sequences of the same type")
        PairwiseAlignmentsSingleSubject(as.character(pattern), as.character(subject),
                                    type = type, substitutionMatrix = substitutionMatrix,
                                    gapOpening = gapOpening, gapExtension = gapExtension,
                                    baseClass = xsbaseclass(pattern))
    }
)

setMethod("PairwiseAlignmentsSingleSubject", signature(pattern = "XStringSet", subject = "missing"),
    function(pattern, subject, type = "global", substitutionMatrix = NULL,
             gapOpening = 0, gapExtension = 1) {
        if (length(pattern) != 2)
            stop("'pattern' must be of length 2 when 'subject' is missing")
        if (diff(nchar(pattern)) != 0)
            stop("'pattern' elements must have the same number of characters")
        PairwiseAlignmentsSingleSubject(as.character(pattern[1]), as.character(pattern[2]),
                                    type = type, substitutionMatrix = substitutionMatrix,
                                    gapOpening = gapOpening, gapExtension = gapExtension,
                                    baseClass = xsbaseclass(pattern))
    }
)

setMethod("PairwiseAlignmentsSingleSubject", signature(pattern = "character", subject = "missing"),
    function(pattern, subject, type = "global", substitutionMatrix = NULL,
             gapOpening = 0, gapExtension = 1, baseClass = "BString") {
        if (length(pattern) != 2)
            stop("'pattern' must be of length 2 when 'subject' is missing")
        if (diff(nchar(pattern)) != 0)
            stop("'pattern' elements must have the same number of characters")
        PairwiseAlignmentsSingleSubject(pattern[1], pattern[2],
                                    type = type, substitutionMatrix = substitutionMatrix,
                                    gapOpening = gapOpening, gapExtension = gapExtension,
                                    baseClass = baseClass)
    }
)

setMethod("PairwiseAlignmentsSingleSubject", signature(pattern = "character", subject = "character"),
    function(pattern, subject, type = "global", substitutionMatrix = NULL,
             gapOpening = 0, gapExtension = 1, baseClass = "BString") {
        newPairwiseAlignments(pattern = pattern, subject = subject, type = type,
                              substitutionMatrix = substitutionMatrix,
                              gapOpening = gapOpening, gapExtension = gapExtension,
                              baseClass = baseClass, pwaClass = "PairwiseAlignmentsSingleSubject")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods.
###

### TODO: Support the 'width' argument.
setMethod("Views", signature = c(subject = "PairwiseAlignmentsSingleSubject"),
          function(subject, start=NULL, end=NULL, width=NULL, names=NULL)
          {
              if (!is.null(width))
                  stop("\"Views\" method for PairwiseAlignmentsSingleSubject objects ",
                       "does not support the 'width' argument yet, sorry!")
              if (is.null(start))
                  start <- NA
              if (all(is.na(start)))
                  start <- start(subject(subject))
              else if (!is.numeric(start) || length(start) > 1)
                  stop("'start' must be either NA or an integer vector of length 1")
              else
                  start <- as.integer(start) + start(subject(subject))
              if (is.null(end))
                  end <- NA
              if (all(is.na(end)))
                  end <- end(subject(subject))
              else if (!is.numeric(end) || length(end) > 1)
                  stop("'end' must be either NA or an integer vector of length 1")
              else
                  end <- as.integer(end) + start(subject(subject))
              tmp <- unaligned(subject(subject))
              if (length(tmp) != 1L)
                  stop("internal error: length(tmp) != 1")
              ans_subject <- tmp[[1L]]
              Views(ans_subject, start=start, end=end, names=names)
          })


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "summary" method.
###

setMethod("summary", "PairwiseAlignmentsSingleSubject", function(object, weight=1L, ...)
          {
              if (!is.numeric(weight) || !(length(weight) %in% c(1, length(object))))
                  stop("'weight' must be an integer vector with length 1 or 'length(object)'")
              if (!is.integer(weight))
                weight <- as.integer(weight)
              if (all(weight == 1))
                  new("PairwiseAlignmentsSingleSubjectSummary",
                      type = type(object),
                      score = score(object),
                      nmatch = nmatch(object),
                      nmismatch = nmismatch(object),
                      ninsertion = nindel(subject(object)),
                      ndeletion = nindel(pattern(object)),
                      coverage = coverage(object),
                      mismatchSummary = mismatchSummary(object))
              else
                  new("PairwiseAlignmentsSingleSubjectSummary",
                      type = type(object),
                      score = rep(score(object), weight),
                      nmatch = rep(nmatch(object), weight),
                      nmismatch = rep(nmismatch(object), weight),
                      ninsertion = nindel(subject(object))[rep(seq_len(length(object)), weight), , drop = FALSE],
                      ndeletion = nindel(pattern(object))[rep(seq_len(length(object)), weight), , drop = FALSE],
                      coverage = coverage(object, weight = weight),
                      mismatchSummary = mismatchSummary(object, weight = weight))
          })

setMethod("show", "PairwiseAlignmentsSingleSubjectSummary", function(object)
          {
              cat(switch(type(object), "global" = "Global", "overlap" = "Overlap",
                         "local" = "Local", "global-local" = "Global-Local",
                         "local-global" = "Local-Global"),
                  " Single Subject Pairwise Alignments\n", sep = "")
              cat("Number of Alignments:  ", length(score(object)), "\n", sep = "")
              cat("\nScores:\n")
              print(summary(score(object)))
              cat("\nNumber of matches:\n")
              print(summary(nmatch(object)))
              n <- min(nrow(mismatchSummary(object)[["subject"]]), 10)
              cat(paste("\nTop", n, "Mismatch Counts:\n"))
              mmtable <- 
                mismatchSummary(object)[["subject"]][
                  order(mismatchSummary(object)[["subject"]][["Count"]],
                        mismatchSummary(object)[["subject"]][["Probability"]],
                        decreasing = TRUE)[seq_len(n)],,drop=FALSE]
              rownames(mmtable) <- NULL
              print(mmtable)
          })

setMethod("type", "PairwiseAlignmentsSingleSubjectSummary", function(x) x@type)
setMethod("score", "PairwiseAlignmentsSingleSubjectSummary", function(x) x@score)
setMethod("nindel", "PairwiseAlignmentsSingleSubjectSummary",
          function(x) new("InDel", insertion = x@ninsertion, deletion = x@ndeletion))
setMethod("length", "PairwiseAlignmentsSingleSubjectSummary", function(x) length(score(x)))
setMethod("nchar", "PairwiseAlignmentsSingleSubjectSummary", 
          function(x, type="chars", allowNA=FALSE)
              unname(nmatch(x) + nmismatch(x) + x@ninsertion[,"WidthSum"] + x@ndeletion[,"WidthSum"]))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "as.character" method.
###

setMethod("aligned", "PairwiseAlignmentsSingleSubject",
          function(x, degap=FALSE, gapCode="-", endgapCode="-") {
              if (degap) {
                  value <- aligned(pattern(x), degap = degap)
              } else {
                  codecX <- xscodec(x)
                  if (is.null(codecX)) {
                      gapCode <- charToRaw(gapCode)
                      endgapCode <- charToRaw(endgapCode)
                  } else {
                      letters2codes <- codecX@codes
                      names(letters2codes) <- codecX@letters
                      gapCode <- as.raw(letters2codes[[gapCode]])
                      endgapCode <- as.raw(letters2codes[[endgapCode]])
                  }
                  value <-
                    .Call2("PairwiseAlignmentsSingleSubject_align_aligned", x, gapCode, endgapCode, PACKAGE="Biostrings")
              }
              value
          })

setMethod("as.character", "PairwiseAlignmentsSingleSubject",
          function(x)
          {
              as.character(aligned(x))
          })

setMethod("toString", "PairwiseAlignmentsSingleSubject", function(x, ...) toString(as.character(x), ...))

setMethod("as.matrix", "PairwiseAlignmentsSingleSubject",
          function(x) {
              as.matrix(aligned(x))
          })

