### ==========================================================================
### PairwiseAlignedFixedSubject objects
### --------------------------------------------------------------------------
### A PairwiseAlignedFixedSubject object contains the result of the pairwise
### alignment of many patterns to one subject.
###
### FIXME: The name of this class is confusing. In the Biostrings context, a
### "fixed subject" (or "fixed pattern") is a subject (or pattern) in which
### the IUPAC ambiguity letters are interpreted literally. See the 'fixed'
### argument of neditStartingAt(), isMatchingStartingAt(), matchPattern(),
### matchPDict(), etc... for the details.

setClass("PairwiseAlignedFixedSubject",
    contains="PairwiseAlignedXStringSet"
)

setClass("PairwiseAlignedFixedSubjectSummary",
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

setGeneric("PairwiseAlignedFixedSubject",
    function(pattern, subject, type = "global", substitutionMatrix = NULL,
             gapOpening = 0, gapExtension = -1, ...)
    standardGeneric("PairwiseAlignedFixedSubject"))
setMethod("PairwiseAlignedFixedSubject", signature(pattern = "XString", subject = "XString"),
    function(pattern, subject, type = "global", substitutionMatrix = NULL,
             gapOpening = 0, gapExtension = -1) {
        if (xsbasetype(pattern) != xsbasetype(subject))
            stop("'pattern' and 'subject' must have the same XString base type")
        PairwiseAlignedFixedSubject(as.character(pattern), as.character(subject),
                                    type = type, substitutionMatrix = substitutionMatrix,
                                    gapOpening = gapOpening, gapExtension = gapExtension,
                                    baseClass = xsbaseclass(pattern))
    }
)

setMethod("PairwiseAlignedFixedSubject", signature(pattern = "XStringSet", subject = "missing"),
    function(pattern, subject, type = "global", substitutionMatrix = NULL,
             gapOpening = 0, gapExtension = -1) {
        if (length(pattern) != 2)
            stop("'pattern' must be of length 2 when 'subject' is missing")
        if (diff(nchar(pattern)) != 0)
            stop("'pattern' elements must have the same number of characters")
        PairwiseAlignedFixedSubject(as.character(pattern[1]), as.character(pattern[2]),
                                    type = type, substitutionMatrix = substitutionMatrix,
                                    gapOpening = gapOpening, gapExtension = gapExtension,
                                    baseClass = xsbaseclass(pattern))
    }
)

setMethod("PairwiseAlignedFixedSubject", signature(pattern = "character", subject = "missing"),
    function(pattern, subject, type = "global", substitutionMatrix = NULL,
             gapOpening = 0, gapExtension = -1, baseClass = "BString") {
        if (length(pattern) != 2)
            stop("'pattern' must be of length 2 when 'subject' is missing")
        if (diff(nchar(pattern)) != 0)
            stop("'pattern' elements must have the same number of characters")
        PairwiseAlignedFixedSubject(pattern[1], pattern[2],
                                    type = type, substitutionMatrix = substitutionMatrix,
                                    gapOpening = gapOpening, gapExtension = gapExtension,
                                    baseClass = baseClass)
    }
)

setMethod("PairwiseAlignedFixedSubject", signature(pattern = "character", subject = "character"),
    function(pattern, subject, type = "global", substitutionMatrix = NULL,
             gapOpening = 0, gapExtension = -1, baseClass = "BString") {
        newPairwiseAlignedXStringSet(pattern = pattern, subject = subject, type = type,
                                     substitutionMatrix = substitutionMatrix,
                                     gapOpening = gapOpening, gapExtension = gapExtension,
                                     baseClass = baseClass, pwaClass = "PairwiseAlignedFixedSubject")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods.
###

setMethod("Views", signature = c(subject = "PairwiseAlignedFixedSubject"),
          function(subject, start=NA, end=NA, names=NULL)
          {
              if (all(is.na(start)))
                  start <- start(subject(subject))
              else if (!is.numeric(start) || length(start) > 1)
                  stop("'start' must be either NA or an integer vector of length 1")
              else
                  start <- as.integer(start) + start(subject(subject))
              if (all(is.na(end)))
                  end <- end(subject(subject))
              else if (!is.numeric(end) || length(end) > 1)
                  stop("'end' must be either NA or an integer vector of length 1")
              else
                  end <- as.integer(end) + start(subject(subject))
              Views(super(unaligned(subject(subject))), start=start, end=end, names=names)
          })


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "summary" method.
###

setMethod("summary", "PairwiseAlignedFixedSubject", function(object, weight=1L, ...)
          {
              if (!is.numeric(weight) || !(length(weight) %in% c(1, length(object))))
                  stop("'weight' must be an integer vector with length 1 or 'length(object)'")
              if (!is.integer(weight))
                weight <- as.integer(weight)
              if (all(weight == 1))
                  new("PairwiseAlignedFixedSubjectSummary",
                      type = type(object),
                      score = score(object),
                      nmatch = nmatch(object),
                      nmismatch = nmismatch(object),
                      ninsertion = nindel(subject(object)),
                      ndeletion = nindel(pattern(object)),
                      coverage = coverage(object),
                      mismatchSummary = mismatchSummary(object))
              else
                  new("PairwiseAlignedFixedSubjectSummary",
                      type = type(object),
                      score = rep(score(object), weight),
                      nmatch = rep(nmatch(object), weight),
                      nmismatch = rep(nmismatch(object), weight),
                      ninsertion = nindel(subject(object))[rep(seq_len(length(object)), weight), , drop = FALSE],
                      ndeletion = nindel(pattern(object))[rep(seq_len(length(object)), weight), , drop = FALSE],
                      coverage = coverage(object, weight = weight),
                      mismatchSummary = mismatchSummary(object, weight = weight))
          })

setMethod("show", "PairwiseAlignedFixedSubjectSummary", function(object)
          {
              cat(switch(type(object), "global" = "Global", "overlap" = "Overlap",
                         "local" = "Local", "global-local" = "Global-Local",
                         "local-global" = "Local-Global"),
                  " Fixed Subject Pairwise Alignment\n", sep = "")
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

setMethod("type", "PairwiseAlignedFixedSubjectSummary", function(x) x@type)
setMethod("score", "PairwiseAlignedFixedSubjectSummary", function(x) x@score)
setMethod("nindel", "PairwiseAlignedFixedSubjectSummary",
          function(x) new("InDel", insertion = x@ninsertion, deletion = x@ndeletion))
setMethod("length", "PairwiseAlignedFixedSubjectSummary", function(x) length(score(x)))
setMethod("nchar", "PairwiseAlignedFixedSubjectSummary", 
          function(x, type="chars", allowNA=FALSE)
              unname(nmatch(x) + nmismatch(x) + x@ninsertion[,"WidthSum"] + x@ndeletion[,"WidthSum"]))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "as.character" method.
###

setMethod("aligned", "PairwiseAlignedFixedSubject",
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
                    .Call("PairwiseAlignedFixedSubject_align_aligned", x, gapCode, endgapCode, PACKAGE="Biostrings")
              }
              value
          })

setMethod("as.character", "PairwiseAlignedFixedSubject",
          function(x)
          {
              as.character(aligned(x))
          })

setMethod("toString", "PairwiseAlignedFixedSubject", function(x, ...) toString(as.character(x), ...))

setMethod("as.matrix", "PairwiseAlignedFixedSubject",
          function(x) {
              as.matrix(aligned(x))
          })
