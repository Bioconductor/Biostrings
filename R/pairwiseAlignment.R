### =========================================================================
### The pairwiseAlignment() generic & related functions
### -------------------------------------------------------------------------
###
### The pairwiseAligment() function provides optimal pairwise alignment of
### the following types:
### - Global alignment
### - Local alignment
### - Overlap alignment
###
### -------------------------------------------------------------------------


XString.pairwiseAlignment <-
function(string1,
         string2,
         quality1 = 1,
         quality2 = 1,
         matchScoring,
         gapOpening = -5,
         gapExtension = -2,
         type = "global")
{
  if (class(string1) != class(string2))
    stop("'string1' and 'string2' must have the same class")
  if (!(length(quality1) %in% c(1, nchar(string1))))
    stop("length(quality1) must be 1 or nchar(string1)")
  if (!(length(quality2) %in% c(1, nchar(string2))))
    stop("length(quality2) must be 1 or nchar(string2)")
  if (any(is.na(quality1)) || any(quality1 < 0) || any(quality1 > 1))
    stop("all elements in 'quality1' must be between 0 and 1")
  if (any(is.na(quality2)) || any(quality2 < 0) || any(quality2 > 1))
    stop("all elements in 'quality2' must be between 0 and 1")
  if (is.character(matchScoring)) {
    if (length(matchScoring) != 1)
      stop("'matchScoring' is a character vector of length != 1")
    tempMatrix <- matchScoring
    matchScoring <- try(getdata(tempMatrix), silent=TRUE)
    if (is(matchScoring, "try-error"))
      stop("unknown scoring matrix \"", tempMatrix, "\"")
  }
  if (!is.matrix(matchScoring) || !is.numeric(matchScoring))
    stop("'matchScoring' must be a numeric matrix")
  if (!identical(rownames(matchScoring), colnames(matchScoring)))
    stop("row and column names differ for matrix 'matchScoring'")
  if (is.null(rownames(matchScoring)))
    stop("matrix 'matchScoring' must have row and column names")
  if (any(duplicated(rownames(matchScoring))))
    stop("matrix 'matchScoring' has duplicated row names")
  matchScoring <-
    matrix(as.double(matchScoring),
           nrow = nrow(matchScoring), ncol = ncol(matchScoring),
           dimnames = dimnames(matchScoring))
  gapOpening <- as.double(- abs(gapOpening))
  if (is.na(gapOpening) || length(gapOpening) != 1)
    stop("'gapOpening' must be a non-positive numeric vector of length 1")
  gapExtension <- as.double(- abs(gapExtension))
  if (is.na(gapExtension) || length(gapExtension) != 1)
    stop("'gapExtension' must be a non-positive numeric vector of length 1")
  type <- match.arg(tolower(type), c("global", "local", "overlap"))
  typeCode <- c("global" = 1L, "local" = 2L, "overlap" = 3L)[[type]]
  if (is.null(codec(string1))) {
    codes <- as.integer(charToRaw(paste(rownames(matchScoring), collapse="")))
    gapCode <- charToRaw("-")
  } else {
    if (!all(rownames(matchScoring) %in% alphabet(string1)))
      stop("matrix 'matchScoring' is incompatible with 'string1' alphabet")
    lettersToCodes <- codec(string1)@codes
    names(lettersToCodes) <- codec(string1)@letters
    codes <- lettersToCodes[rownames(matchScoring)]
    gapCode <- as.raw(lettersToCodes[["-"]])
  }
  lookupTable <- buildLookupTable(codes, 0:(nrow(matchScoring) - 1))
  answer <- .Call("align_pairwiseAlignment",
                  string1,
                  string2,
                  quality1,
                  quality2,
                  gapCode,
                  typeCode,
                  lookupTable,
                  matchScoring,
                  dim(matchScoring),
                  gapOpening,
                  gapExtension,
                  PACKAGE="Biostrings")
  align1 <-
    new(class(string1),
        xdata = answer[["align1"]],
        length = length(answer[["align1"]]))
  align2 <-
    new(class(string2),
        xdata = answer[["align2"]],
        length = length(answer[["align2"]]))
  return(new("XStringAlign",
             align1 = align1,
             align2 = align2,
             quality1 = quality1,
             quality2 = quality2,
             type = type,
             matchScoring = matchScoring,
             gapOpening = gapOpening,
             gapExtension = gapExtension,
             score = answer[["score"]]))
}


setGeneric("pairwiseAlignment", signature = c("string1", "string2"),
           function(string1, string2, quality1 = 1, quality2 = 1,
                    matchScoring, gapOpening = -5, gapExtension = -2,
                    type = "global")
           standardGeneric("pairwiseAlignment"))

setMethod("pairwiseAlignment",
          signature(string1 = "character", string2 = "character"),
          function(string1, string2, quality1 = 1, quality2 = 1,
                   matchScoring, gapOpening = -5, gapExtension = -2,
                   type = "global")
          XString.pairwiseAlignment(BString(string1), BString(string2),
                                    quality1 = quality1,
                                    quality2 = quality2,
                                    matchScoring = matchScoring,
                                    gapExtension = gapExtension,
                                    gapOpening = gapOpening,
                                    type = type))

setMethod("pairwiseAlignment",
          signature(string1 = "character", string2 = "XString"),
          function(string1, string2, quality1 = 1, quality2 = 1,
                   matchScoring, gapOpening = -5, gapExtension = -2,
                   type = "global")
          XString.pairwiseAlignment(XString(class(string2), string1), string2,
                                    quality1 = quality1,
                                    quality2 = quality2,
                                    matchScoring = matchScoring,
                                    gapExtension = gapExtension,
                                    gapOpening = gapOpening,
                                    type = type))

setMethod("pairwiseAlignment",
          signature(string1 = "XString", string2 = "character"),
          function(string1, string2, quality1 = 1, quality2 = 1,
                   matchScoring, gapOpening = -5, gapExtension = -2,
                   type = "global")
          XString.pairwiseAlignment(string1, XString(class(string1), string2),
                                    quality1 = quality1,
                                    quality2 = quality2,
                                    matchScoring = matchScoring,
                                    gapExtension = gapExtension,
                                    gapOpening = gapOpening,
                                    type = type))

setMethod("pairwiseAlignment",
          signature(string1 = "XString", string2 = "XString"),
          function(string1, string2, quality1 = 1, quality2 = 1,
                   matchScoring, gapOpening = -5, gapExtension = -2,
                   type = "global")
          XString.pairwiseAlignment(string1, string2,
                                    quality1 = quality1,
                                    quality2 = quality2,
                                    matchScoring = matchScoring,
                                    gapExtension = gapExtension,
                                    gapOpening = gapOpening,
                                    type = type))
