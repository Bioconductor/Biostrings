### =========================================================================
### Reading a mask from a file
### --------------------------
###
### From a RepeatMasker .out file (for chrM in ce2):
###   library(BSgenome.Celegans.UCSC.ce2)
###   file1 <- system.file("extdata", "ce2chrM.fa.out", package="Biostrings")
###   mask1 <- read.rmMask(file1, width=length(Celegans$chrM)) 
### From a Tandem Repeats Finder .bed file (for chrM in ce2):
###   file2 <- system.file("extdata", "ce2chrM.bed", package="Biostrings")
###   mask2 <- read.trfMask(file2, width=length(Celegans$chrM)) 
###
### -------------------------------------------------------------------------


read.rmMask <- function(file, width, with.ID=FALSE)
{
    if (!isSingleNumber(width))
        stop("'width' must be a single integer")
    if (!is.integer(width))
        width <- as.integer(width)
    if (!isTRUEorFALSE(with.ID))
        stop("'with.ID' must be 'TRUE' or 'FALSE'")
    ALLCOLS <- c(
        `SW_score`="integer",
        `perc_div`="numeric",
        `perc_del`="numeric",
        `perc_ins`="numeric",
        `query_sequence`="character",
        `begin_in_query`="integer",
        `end_in_query`="integer",
        `left_in_query`="character",
        `strand`="character",           # Looks like the strand but I'm not sure!
        `matching_repeat`="character",
        `repeat_class_or_family`="character",
        `begin_in_repeat`="integer",
        `end_in_repeat`="integer",
        `left_in_repeat`="character",
        `ID`="character"
    )
    COLS <- c("begin_in_query", "end_in_query", "ID")
    ALLCOLS[!(names(ALLCOLS) %in% COLS)] <- "NULL"
    data <- read.table(file,
                       col.names=names(ALLCOLS),
                       colClasses=ALLCOLS,
                       skip=3,
                       check.names=FALSE)
    ranges <- IRanges(start=data$begin_in_query, end=data$end_in_query)
    if (with.ID) {
        names(ranges) <- data$ID
        if (isNotStrictlySorted(start(ranges)))
            ranges <- ranges[order(start(ranges))]
        if (!isNormal(ranges))
            stop("cannot keep the repeat IDs when some repeats overlap")
        nir1 <- as(ranges, "NormalIRanges")
    } else {
        nir1 <- toNormalIRanges(ranges)
    }
    new("MaskCollection", nir_list=list(nir1), width=width, active=TRUE, NAMES="rm")
}

read.trfMask <- function(file, width)
{
    if (!isSingleNumber(width))
        stop("'width' must be a single integer")
    if (!is.integer(width))
        width <- as.integer(width)
    ALLCOLS <- c(
        `chrom`="character",
        `chromStart`="integer",
        `chromEnd`="integer",
        `name`="character",
        `period`="integer",
        `copyNum`="numeric",
        `consensusSize`="integer",
        `perMatch`="integer",
        `perIndel`="integer",
        `score`="integer",
        `A`="integer",
        `C`="integer",
        `G`="integer",
        `T`="integer",
        `entropy`="numeric",
        `sequence`="character"
    )
    COLS <- c("chromStart", "chromEnd")
    ALLCOLS[!(names(ALLCOLS) %in% COLS)] <- "NULL"
    data <- read.table(file,
                       col.names=names(ALLCOLS),
                       colClasses=ALLCOLS,
                       check.names=FALSE)
    ranges <- IRanges(start=data$chromStart+1, end=data$chromEnd)
    nir1 <- toNormalIRanges(ranges)
    #name1 <- "Tandem Repeats Finder (period<=12)"
    new("MaskCollection", nir_list=list(nir1), width=width, active=TRUE, NAMES="trf")
}

