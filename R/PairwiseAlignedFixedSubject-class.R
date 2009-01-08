### ==========================================================================
### PairwiseAlignedFixedSubject objects
### --------------------------------------------------------------------------
### A PairwiseAlignedFixedSubject object contains the result of the pairwise
### alignment of many patterns to one subject.


setClass("PairwiseAlignedFixedSubject",
    representation(
        "PairwiseAlignedXStringSet"
    )
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
        if (baseXStringSubtype(pattern) != baseXStringSubtype(subject))
            stop("'pattern' and 'subject' must have the same XString base subtype")
        PairwiseAlignedFixedSubject(as.character(pattern), as.character(subject),
                                    type = type, substitutionMatrix = substitutionMatrix,
                                    gapOpening = gapOpening, gapExtension = gapExtension,
                                    baseClass = baseXStringSubtype(pattern))
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
                                    baseClass = baseXStringSubtype(pattern))
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
### Initialization.
###

setMethod("initialize", "PairwiseAlignedFixedSubject",
    function(.Object, pattern, subject, type, score, substitutionArray,
             gapOpening, gapExtension, check = TRUE)
    {
        if (!identical(class(unaligned(pattern)), class(unaligned(subject))))
            stop("'unaligned(pattern)' and 'unaligned(subject)' must be XString objects of the same subtype")
        if (length(type) != 1 || !(type %in% c("global", "local", "overlap", "patternOverlap", "subjectOverlap")))
            stop("'type' must be one of 'global', 'local', 'overlap', 'patternOverlap', or 'subjectOverlap'")
        if (length(pattern) != length(subject))
            stop("'length(pattern)' must equal 'length(subject)'")
        gapOpening <- as.double(- abs(gapOpening))
        if (length(gapOpening) != 1 || is.na(gapOpening))
            stop("'gapOpening' must be a non-positive numeric vector of length 1")
        gapExtension <- as.double(- abs(gapExtension))
        if (length(gapExtension) != 1 || is.na(gapExtension))
            stop("'gapExtension' must be a non-positive numeric vector of length 1")
        slot(.Object, "pattern", check = check) <- pattern
        slot(.Object, "subject", check = check) <- subject
        slot(.Object, "type", check = check) <- type
        slot(.Object, "score", check = check) <- score
        slot(.Object, "substitutionArray", check = check) <- substitutionArray
        slot(.Object, "gapOpening", check = check) <- gapOpening
        slot(.Object, "gapExtension", check = check) <- gapExtension
        .Object
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.PairwiseAlignedFixedSubject <- function(object)
{
    message <- character(0)
    if (!identical(class(unaligned(pattern(object))), class(unaligned(subject(object)))))
        message <- c(message, "'unaligned(pattern)' and 'unaligned(subject)' must be XString objects of the same subtype")
    if (length(object@type) != 1 || !(object@type %in% c("global", "local", "overlap", "patternOverlap", "subjectOverlap")))
        message <- c(message, "'type' must be one of 'global', 'local', 'overlap', 'patternOverlap', or 'subjectOverlap'")
    if (length(pattern) != length(subject))
        message <- c(message, "'length(pattern)' must equal 'length(subject)'")
    if (length(object@gapOpening) != 1 || is.na(object@gapOpening))
        message <- c(message, "'gapOpening' must be a non-positive numeric vector of length 1")
    if (length(object@gapExtension) != 1 || is.na(object@gapExtension))
        message <- c(message, "'gapExtension' must be a non-positive numeric vector of length 1")
    if (length(message) == 0)
        message <- NULL
    message
}

setValidity("PairwiseAlignedFixedSubject",
    function(object)
    {
        problems <- .valid.PairwiseAlignedFixedSubject(object)
        if (is.null(problems)) TRUE else problems
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
                         "patternOverlap" = "Pattern Overlap", "subjectOverlap" = "Subject Overlap",
                         "local" = "Local"), " Fixed Subject Pairwise Alignment\n", sep = "")
              cat("Number of Alignments:  ", length(score(object)), "\n", sep = "")
              cat("\nScores:\n")
              print(summary(score(object)))
              cat("\nNumber of matches:\n")
              print(summary(nmatch(object)))
              n <- min(nrow(mismatchSummary(object)[["subject"]]), 10)
              cat(paste("\nTop", n, "Mismatch Counts:\n"))
              print(mismatchSummary(object)[["subject"]][
                order(mismatchSummary(object)[["subject"]][["Count"]],
                      mismatchSummary(object)[["subject"]][["Probability"]],
                      decreasing = TRUE)[seq_len(n)],,drop=FALSE])
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
          function(x, degap = FALSE) {
              if (degap) {
                  value <- aligned(pattern(x), degap = degap)
              } else {
                  codecX <- codec(x)
                  if (is.null(codecX)) {
                      gapCode <- charToRaw("-")
                  } else {
                      letters2codes <- codecX@codes
                      names(letters2codes) <- codecX@letters
                      gapCode <- as.raw(letters2codes[["-"]])
                  }
                  value <-
                    .Call("PairwiseAlignedFixedSubject_align_aligned", x, gapCode, PACKAGE="Biostrings")
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
