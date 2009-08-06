### =========================================================================
### Genome-oriented methods for RangedData classes
### -------------------------------------------------------------------------

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods.
###

setGeneric("genome", function(x, ...) standardGeneric("genome"))
setMethod("genome", "RangedData", function(x) universe(x))

setGeneric("genome<-", function(x, value) standardGeneric("genome<-"))
setReplaceMethod("genome", "RangedData",
                 function(x, value) {
                   genome(x@ranges) <- value
                   x
                 })

setGeneric("chrom", function(x, ...) standardGeneric("chrom"))
setMethod("chrom", "RangedData", function(x) {
  chrom(ranges(x))
})

## score: a common track column

setGeneric("score<-", function(x, ..., value) standardGeneric("score<-"))

setMethod("score", "RangedData", function(x) {
  score <- x[["score"]]
  if (is.null(score) && ncol(x) > 0 && is.numeric(x[[1]]))
    score <- x[[1]]
  score
})
setReplaceMethod("score", "RangedData", function(x, value) {
  if (!is.numeric(value))
    stop("score must be numeric")
  if (length(value) != nrow(x))
    stop("number of scores must equal the number of rows")
  x[["score"]] <- value
  x
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###

GenomicData <- function(ranges, ..., strand = NULL, chrom = NULL, genome = NULL)
{
  if (is(ranges, "data.frame") || is(ranges, "DataTable")) {
    colnames(ranges)[match("chrom", colnames(ranges))] <- "space"
  }
  if (!is(ranges, "Ranges")) {
    return(RangedData(ranges)) # direct coercion
  }
  if (length(chrom) > length(ranges))
    stop("length of 'chrom' greater than length of 'ranges'")
  if (length(chrom) > 0 && (length(ranges) %% length(chrom) != 0))
    stop("length of 'ranges' not a multiple of 'chrom' length")
  if (!is.null(genome) && (length(genome) != 1 || !is.character(genome)))
    stop("'genome' must be a single string")
  if (!is.null(strand)) {
    if (!all(strand[!is.na(strand)] %in% levels(strand())))
      stop("strand values should be 'NA', '-', '+' or '*'")
    strand <- factor(strand, levels(strand()))
    RangedData(ranges, ..., strand = strand, space = chrom, universe = genome)
  } else RangedData(ranges, ..., space = chrom, universe = genome)
}

### =========================================================================
### Genome-oriented methods for RangesList classes
### -------------------------------------------------------------------------

setMethod("genome", "RangesList", function(x) universe(x))

setReplaceMethod("genome", "RangesList",
                 function(x, value) {
                   universe(x) <- value
                   x
                 })

setMethod("chrom", "RangesList", function(x) {
  chrom <- names(x)
  if (!is.null(chrom))
    chrom <- rep(factor(chrom, chrom), unlist(lapply(x, length)))
  chrom
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

GenomicRanges <- function(start = integer(), end = integer(), chrom = NULL,
                          genome = NULL)
{
  ir <- IRanges(start, end)
  if (!is.null(chrom)) {
    if (!is.factor(chrom))
      chrom <- factor(chrom, unique(chrom))
    if (!length(ir))
      chrom <- chrom[FALSE]
    if (length(chrom) != length(ir))
      stop("length of 'chrom' must equal that of 'start' and 'end'",
           " unless 'start' and 'end' are zero length")
    rl <- split(ir, chrom)
  } else rl <- RangesList(ir)
  universe(rl) <- genome
  rl
}
