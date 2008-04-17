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


### R-implementation of global, local, and overlap alignment.
### 'string1' and 'string2' must be character vectors of length 1.
.pairwiseAlignment <-
function(string1,
         string2,
         matchScores,
         gapOpening = -5L,
         gapExtension = -2L,
         type = "global")
{
  getIndex <- function(i, j, ncol) {
    return((i - 1) * ncol + j)
  }

  ## Step 0:  Check inputs
  if (length(string1) != 1 || length(string2) != 1) {
    stop("objects 'string1' and 'string2' each need to be vectors of length 1")
  }

  stringElements1 <- safeExplode(string1)
  stringElements2 <- safeExplode(string2)
  if (!all(unique(stringElements1) %in% rownames(matchScores)) &&
      !all(unique(stringElements2) %in% colnames(matchScores))) {
    stop(paste("not all symbols used in 'string1' and 'string2'",
               "are present in 'matchScores' matrix"))
  }

  ## Step 1:  Get information on input strings
  n1 <- length(stringElements1)
  n2 <- length(stringElements2)
  gapExtension <- as.integer(- abs(gapExtension))
  if(gapOpening != 0L) {
    stop("affine gaps are currently not supported")
  }
  type <- match.arg(tolower(type), c("global", "local", "overlap"))

  ## Step 2:  Create objects for scores and traceback values
  fMatrix <- matrix(0L, n1 + 1, n2 + 1)
  if(type == "global") {
    for(i in 1:(n1 + 1))
      fMatrix[i, 1] <- as.integer((i - 1) * gapExtension)
    for(j in 1:(n2 + 1))
      fMatrix[1, j] <- as.integer((j - 1) * gapExtension)
  }

  ## Traceback values:
  ##   1 = a match,
  ##   2 = a gap in string 2,
  ##   3 = a gap in string 1,
  ##   4 = 0 ("local" only)
  traceMatrix <- vector("list", (n1 + 1) * (n2 + 1))
  possibleTraceValues <- vector("list", 4)
  for(i in seq(along = possibleTraceValues)) {
    possibleTraceValues[[i]] <- new.env(parent = emptyenv())
    possibleTraceValues[[i]][["value"]] <- as.integer(i)
  }

  ## Step 3:  Generate scores and traceback values
  startIndex <- -1L
  startScore <- - .Machine[["integer.max"]]
  minPossibleScore <- ifelse(type == "local", 0, - .Machine[["integer.max"]])
  for (i in seq_len(n1)) {
    for (j in seq_len(n2)) {
      candidates <-
        c(fMatrix[i, j] + matchScores[stringElements1[i], stringElements2[j]],
          fMatrix[i, j + 1] + gapExtension,
          fMatrix[i + 1, j] + gapExtension,
          minPossibleScore)
      index <- getIndex(i + 1, j + 1, n2 + 1)
      score <- max(candidates)
      if(type == "local") {
        if(score == startScore) {
          startIndex <- c(startIndex, index)
        } else if(score > startScore) {
          startIndex <- index
          startScore <- score
        }
      }
      fMatrix[i + 1, j + 1] <- score
      traceMatrix[[index]] <- possibleTraceValues[[which.max(candidates)]]
    }
  }
  if(type == "global") {
    startIndex <- getIndex(n1 + 1, n2 + 1, n2 + 1)
    startScore <- fMatrix[n1 + 1, n2 + 1]
  } else if(type == "overlap") {
    for (i in seq_len(n1)) {
      index <- getIndex(i + 1, n2 + 1, n2 + 1)
      score <- fMatrix[i + 1, n2 + 1]
      if(score == startScore) {
        startIndex <- c(startIndex, index)
      } else if(score > startScore) {
        startIndex <- index
        startScore <- score
      }
    }
    for (j in seq_len(n2)) {
      index <- getIndex(n1 + 1, j + 1, n2 + 1)
      score <- fMatrix[n1 + 1, j + 1]
      if(score == startScore) {
        startIndex <- c(startIndex, index)
      } else if(score > startScore) {
        startIndex <- index
        startScore <- score
      }
    }
  }

  ## Step 4:  Get a starting location for the traceback
  startRow <- 1 + (startIndex[1] - 1) %/% (n2 + 1)
  startCol <- startIndex[1] - (startRow - 1) * (n2 + 1)
  align1 <- align2 <- character(0)
  if(type == "overlap") {
    if(startRow == n1 + 1) {
      align1 <- rep("-", n2 + 1 - startCol)
      if(startCol == n2 + 1)
        align2 <- character(0)
      else
        align2 <- stringElements2[startCol:n2]
    } else {
      if(startRow == n1 + 1)
        align1 <- character(0)
      else
        align1 <- stringElements1[startRow:n1]
      align2 <- rep("-", n1 + 1 - startRow)
    }
  }

  ## Step 5:  Traceback through the score matrix
  currentRow <- startRow
  currentCol <- startCol
  while(((type == "global" || type == "overlap") &&
         (currentRow > 1 && currentCol > 1)) ||
        ((type == "local") && (fMatrix[currentRow, currentCol] > 0))) {
    currentTrace <-
      traceMatrix[[getIndex(currentRow, currentCol, n2 + 1)]][["value"]]
    align1 <-
      c(ifelse(currentTrace == 3, "-", stringElements1[currentRow - 1]),
        align1)
    align2 <-
      c(ifelse(currentTrace == 2, "-", stringElements2[currentCol - 1]),
        align2)
    currentRow <- currentRow - as.integer(currentTrace != 3)
    currentCol <- currentCol - as.integer(currentTrace != 2)
  }
  if(type == "global" || type == "overlap") {
    if(currentCol > 1) {
      align1 <- c(rep("-", currentCol - 1), align1)
      align2 <- c(stringElements2[1:(currentCol - 1)], align2)
    }
    if(currentRow > 1) {
      align1 <- c(stringElements1[1:(currentRow - 1)], align1)
      align2 <- c(rep("-", currentRow - 1), align2)
    }
  }

  ## Step 6: Generate the return value
  answer <-
    c(align1 = paste(align1, collapse = ""),
      align2 = paste(align2, collapse = ""))
  attr(answer, "type") <- type
  attr(answer, "score") <- fMatrix[startRow, startCol]
  attr(answer, "fMatrix") <- fMatrix
  attr(answer, "traceMatrix") <- traceMatrix
  class(answer) <- "PairwiseAlignment"
  return(answer)
}


print.PairwiseAlignment <- function(x, ...)
{
  print(matrix(c(x["align1"], x["align2"]), nrow = 2, ncol = 1,
               dimnames = list(c("1", "2"), "Alignment1")))
  cat(paste("Type:  ",
            switch(attr(x, "type"),
                   "global" = "Global",
                   "overlap" = "Overlap",
                   "local" = "Local"),
            "\n", sep = ""))
  cat(paste("Score:  ", attr(x, "score"), sep = ""),"\n")
  return(invisible(x))
}



### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### C-implementation of the Needleman-Wunsch algo.
###
### The various "needwunsQS" methods below don't call
### .Call("align_needwunsQS", ...) directly but call XString.needwunsQS()
### instead which itself calls .Call("align_needwunsQS", ...).
### Some quick testing shows that, depending on the size of the strings to
### align, this C version is 100 to 1000 times faster than the above
### .needwunsQS().
### 's1' and 's2' must be XString objects of the same subtype.
### Return an XStringAlign object where the "al1" and "al2" slots contain the
### aligned versions of 's1' and 's2'.
XString.pairwiseAlignment <-
function(string1,
         string2,
         matchScores,
         gapOpening = -5L,
         gapExtension = -2L,
         type = "global")
{
  if (class(string1) != class(string2))
    stop("'string1' and 'string2' are not of the same class")
  ## The reason we test 'length(string1)' and 'length(string2)'
  ## separately is to avoid a "NAs produced by integer overflow"
  ## in the product of the two.
  if (length(string1) > 20000L || length(string2) > 20000L ||
      length(string1) * length(string2) > 20000L * 20000L)
    stop("'length(string1) * length(string2)' is too big (> 4e+08)")
  if (is.character(matchScores)) {
    if (length(matchScores) != 1)
      stop("'matchScores' is a character vector of length != 1")
    tempMatrix <- matchScores
    matchScores <- try(getdata(tempMatrix), silent=TRUE)
    if (is(matchScores, "try-error"))
      stop("unknown scoring matrix \"", tempMatrix, "\"")
  }
  if (!is.matrix(matchScores) || !is.integer(matchScores))
    stop("'matchScores' must be a matrix of integers")
  if (!identical(rownames(matchScores), colnames(matchScores)))
    stop("row and column names differ for matrix 'matchScores'")
  if (is.null(rownames(matchScores)))
    stop("matrix 'matchScores' must have row and column names")
  if (any(duplicated(rownames(matchScores))))
    stop("matrix 'matchScores' has duplicated row names")
  gapOpening <- as.integer(- abs(gapOpening))
  if (is.na(gapOpening) || length(gapOpening) != 1)
    stop("'gapOpening' must be a non-positive integer vector of length 1")
  gapExtension <- as.integer(- abs(gapExtension))
  if (is.na(gapExtension) || length(gapExtension) != 1)
    stop("'gapExtension' must be a non-positive integer vector of length 1")
  type <- match.arg(tolower(type), c("global", "local", "overlap"))
  typeCode <- c("global" = 1L, "local" = 2L, "overlap" = 3L)[[type]]
  if (is.null(codec(string1))) {
    codes <- as.integer(charToRaw(paste(rownames(matchScores), collapse="")))
    gapCode <- charToRaw("-")
  } else {
    if (!all(rownames(matchScores) %in% alphabet(string1)))
      stop("matrix 'matchScores' is incompatible with 'string1' alphabet")
    lettersToCodes <- codec(string1)@codes
    names(lettersToCodes) <- codec(string1)@letters
    codes <- lettersToCodes[rownames(matchScores)]
    gapCode <- as.raw(lettersToCodes[["-"]])
  }
  lookupTable <- buildLookupTable(codes, 0:(nrow(matchScores) - 1))
  answer <- .Call("R_pairwiseAlignment",
                  string1,
                  string2,
                  matchScores,
                  dim(matchScores),
                  lookupTable,
                  gapOpening,
                  gapExtension,
                  gapCode,
                  typeCode,
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
             type = type,
             score = answer[["score"]]))
}


setGeneric("pairwiseAlignment", signature = c("string1", "string2"),
           function(string1, string2, matchScores,
                    gapOpening = -5L, gapExtension = -2L,
                    type = "global")
           standardGeneric("pairwiseAlignment"))

setMethod("pairwiseAlignment",
          signature(string1 = "character", string2 = "character"),
          function(string1, string2, matchScores,
                   gapOpening = -5L, gapExtension = -2L,
                   type = "global")
          XString.pairwiseAlignment(BString(string1), BString(string2),
                                    matchScores = matchScores,
                                    gapExtension = gapExtension,
                                    gapOpening = gapOpening,
                                    type = type))

setMethod("pairwiseAlignment",
          signature(string1 = "character", string2 = "XString"),
          function(string1, string2, matchScores,
                   gapOpening = -5L, gapExtension = -2L,
                   type = "global")
          XString.pairwiseAlignment(XString(class(string2), string1), string2,
                                    matchScores = matchScores,
                                    gapExtension = gapExtension,
                                    gapOpening = gapOpening,
                                    type = type))

setMethod("pairwiseAlignment",
          signature(string1 = "XString", string2 = "character"),
          function(string1, string2, matchScores,
                   gapOpening = -5L, gapExtension = -2L,
                   type = "global")
          XString.pairwiseAlignment(string1, XString(class(string1), string2),
                                    matchScores = matchScores,
                                    gapExtension = gapExtension,
                                    gapOpening = gapOpening,
                                    type = type))

setMethod("pairwiseAlignment",
          signature(string1 = "XString", string2 = "XString"),
          function(string1, string2, matchScores,
                   gapOpening = -5L, gapExtension = -2L,
                   type = "global")
          XString.pairwiseAlignment(string1, string2,
                                    matchScores = matchScores,
                                    gapExtension = gapExtension,
                                    gapOpening = gapOpening,
                                    type = type))
