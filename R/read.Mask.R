### =========================================================================
### Read a mask from a file
### -----------------------
###
### From an UCSC "lift" file (for hg18):
###   file1 <- system.file("extdata", "hg18liftAll.lft", package="Biostrings")
###   mask1 <- read.liftMask(file1, "chr1")
### From a RepeatMasker .out file (for chrM in ce2):
###   library(BSgenome.Celegans.UCSC.ce2)
###   file2 <- system.file("extdata", "ce2chrM.fa.out", package="Biostrings")
###   mask2 <- read.rmMask(file2, width=length(Celegans$chrM)) 
### From a Tandem Repeats Finder .bed file (for chrM in ce2):
###   file3 <- system.file("extdata", "ce2chrM.bed", package="Biostrings")
###   mask3 <- read.trfMask(file3, width=length(Celegans$chrM)) 
###
### -------------------------------------------------------------------------


read.liftMask <- function(file, seqname=NA, width=NA)
{
    if (!isSingleStringOrNA(seqname))
        stop("'seqname' must be a single string or 'NA'")
    if (!isSingleNumberOrNA(width))
        stop("'width' must be a single integer or 'NA'")
    if (!is.integer(width))
        width <- as.integer(width)
    ALLCOLS <- c(
        `offset`="integer",
        `xxxx`="NULL",  # not sure how to call this
        `width`="integer",
        `seqname`="character",
        `seqlen`="integer"
    )
    data <- read.table(file,
                       col.names=names(ALLCOLS),
                       colClasses=ALLCOLS,
                       check.names=FALSE)
    if (is.na(seqname)) {
        found_seqnames <- paste("\"", unique(data$seqname), "\"", sep="")
        found_seqnames <- paste(found_seqnames, collapse=", ")
        stop("seqnames found in this file: ", found_seqnames)
    }
    data <- data[data$seqname %in% seqname, ]
    if (nrow(data) == 0) {
        if (is.na(width))
            stop("unknown sequence \"", seqname, "\", ",
                 "please specify the width of the empty mask to return")
        warning("unknown sequence \"", seqname, "\", returning empty mask")
        ans <- Mask(width, start=integer(0), width=integer(0))
        names(ans) <- "inter-contig gaps (empty)"
        return(ans)
    }
    ## Sanity checks
    seqlen0 <- unique(data$seqlen)
    if (length(seqlen0) != 1)
        stop("broken \"lift\" file: contains different lengths ",
             "for sequence \"", seqname, "\"")
    if (!is.na(width) && width != seqlen0)
        stop("when specified, 'width' must match the length found ",
             "in the file for sequence \"", seqname, "\"")
    contigs0 <- IRanges(start=data$offset+1, width=data$width)
    contigs1 <- toNormalIRanges(contigs0)
    if (length(contigs1) != length(contigs0))
        warning("some contigs are adjacent or overlapping")
    contigs <- Mask(seqlen0, start=start(contigs1), width=width(contigs1))
    ans <- gaps(contigs)
    names(ans) <- "inter-contig gaps"
    ans
}

read.rmMask <- function(file, width, use.IDs=FALSE)
{
    if (!isSingleNumber(width))
        stop("'width' must be a single integer")
    if (!is.integer(width))
        width <- as.integer(width)
    if (!isTRUEorFALSE(use.IDs))
        stop("'use.IDs' must be 'TRUE' or 'FALSE'")
    ## For a description of RepeatMasker output format:
    ##   http://www.repeatmasker.org/webrepeatmaskerhelp.html
    ALLCOLS <- c(
        `SW_score`="integer",
        `perc_div`="numeric",
        `perc_del`="numeric",
        `perc_ins`="numeric",
        `query_sequence`="character",
        `begin_in_query`="integer",
        `end_in_query`="integer",
        `left_in_query`="character",
        `C`="character", 
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
    if (use.IDs) {
        names(ranges) <- data$ID
        if (isNotStrictlySorted(start(ranges)))
            ranges <- ranges[order(start(ranges))]
        if (!isNormal(ranges))
            stop("cannot keep the repeat IDs when some repeats overlap")
        nir1 <- asNormalIRanges(ranges)
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

