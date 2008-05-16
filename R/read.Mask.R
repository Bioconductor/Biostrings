### =========================================================================
### Read a mask from a file
### -----------------------
###
### From an NCBI "agp" file (for chrY in hs b36v3):
###   library(BSgenome.Hsapiens.NCBI.b36v3)
###   file1 <- system.file("extdata", "hs_b36v3_chrY.agp", package="Biostrings")
###   mask1 <- read.agpMask(file1, length(Hsapiens$chrY), seqname="chrY")
###
### From an UCSC "gap" file (for chrY in hg18):
###   library(BSgenome.Hsapiens.UCSC.hg18)
###   file2 <- system.file("extdata", "chrY_gap.txt", package="Biostrings")
###   mask2 <- read.gapMask(file2, length(Hsapiens$chrY), seqname="chrY")
###
### From an UCSC "lift" file (for hg18):
###   file3 <- system.file("extdata", "hg18liftAll.lft", package="Biostrings")
###   mask3 <- read.liftMask(file3, seqname="chr1")
###
### From a RepeatMasker .out file (for chrM in ce2):
###   library(BSgenome.Celegans.UCSC.ce2)
###   file4 <- system.file("extdata", "ce2chrM.fa.out", package="Biostrings")
###   mask4 <- read.rmMask(file4, length(Celegans$chrM)) 
###
### From a Tandem Repeats Finder .bed file (for chrM in ce2):
###   file5 <- system.file("extdata", "ce2chrM.bed", package="Biostrings")
###   mask5 <- read.trfMask(file5, length(Celegans$chrM)) 
###
### -------------------------------------------------------------------------


.emptyMaskWithWarning <- function(mask_name, seqname, width)
{
    warning("no ", mask_name, " for sequence \"", seqname, "\" in this file, ",
            "returning empty mask")
    ans <- Mask(width, start=integer(0), width=integer(0))
    names(ans) <- paste(mask_name, "(empty)")
    ans
}

.read.MaskFromAgpOrGap <- function(agp_or_gap, file, width, seqname,
                                   gap.types, use.gap.types)
{
    if (!isSingleNumber(width))
        stop("'width' must be a single integer")
    if (!is.integer(width))
        width <- as.integer(width)
    if (!isSingleString(seqname))
        stop("'seqname' must be a single string")
    if (!is.null(gap.types) && (!is.character(gap.types)
                                || any(is.na(gap.types))
                                || any(duplicated(gap.types))))
        stop("'gap.types' must be 'NULL' or a character vector ",
             "with no NAs and no duplicated")
    if (!isTRUEorFALSE(use.gap.types))
        stop("'use.gap.types' must be 'TRUE' or 'FALSE'")
    if (agp_or_gap == "agp") {
        ALL_COLS <- c(
            `chrom`="character",
            `chr_start`="integer",
            `chr_stop`="integer",
            `part_no`="integer",
            `part_type`="character",
            `gap_len`="character",
            `gap_type`="character",
            `linkage`="character",
            `empty`="character"
        )
    } else if (agp_or_gap == "gap") {
        ALL_COLS <- c(
            `bin`="integer",
            `chrom`="character",
            `chr_start`="integer",
            `chr_stop`="integer",
            `part_no`="integer",
            `part_type`="character",
            `gap_len`="integer",
            `gap_type`="character",
            `bridge`="character"
        )
    } else {
        stop("Biostrings internal error: please report")
    }
    COLS <- c(
        "chrom",
        "chr_start",
        "chr_stop",
        "part_type",
        "gap_len",
        "gap_type"
    )
    ALL_COLS[!(names(ALL_COLS) %in% COLS)] <- "NULL"
    data <- read.table(file,
                       sep="\t",
                       col.names=names(ALL_COLS),
                       colClasses=ALL_COLS,
                       check.names=FALSE,
                       fill=TRUE)
    if (seqname == "?") {
        found_seqnames <- paste("\"", unique(data$chrom), "\"", sep="")
        found_seqnames <- paste(found_seqnames, collapse=", ")
        stop("seqnames found in this file: ", found_seqnames)
    }
    data <- data[data$chrom == seqname, ]
    ii <- data$part_type == "N"
    if (agp_or_gap == "agp") {
        data <- data[ii, ]
    } else if (!all(ii)) {
        warning("gap file contains gaps with a part_type that is not N")
    }
    if (length(gap.types) == 1 && gap.types == "?") {
        found_types <- paste("\"", unique(data$gap_type), "\"", sep="")
        found_types <- paste(found_types, collapse=", ")
        stop("gap types found in this file for sequence \"", seqname, "\": ", found_types)
    }
    mask_name <- "assembly gaps"
    if (!is.null(gap.types)) {
        data <- data[data$gap_type %in% gap.types, ]
        mask_name <- paste(mask_name, " [type=", paste(gap.types, collapse="|"), "]", sep="")
    }
    if (nrow(data) == 0)
        return(.emptyMaskWithWarning(mask_name, seqname, width))
    if (agp_or_gap == "agp")
        ranges_start <- data$chr_start
    else
        ranges_start <- data$chr_start + 1L
    ranges <- IRanges(start=ranges_start, width=as.integer(data$gap_len))
    ## Sanity check
    if (!identical(end(ranges), data$chr_stop))
        stop("broken \"", agp_or_gap, "\" file: contains inconsistent ",
             "chr_start/chr_stop/gap_len values ",
             "for assembly gaps in sequence \"", seqname, "\"")
    if (use.gap.types) {
        names(ranges) <- data$gap_type
        if (isNotStrictlySorted(start(ranges)))
            ranges <- ranges[order(start(ranges))]
        if (!isNormal(ranges))
            stop("cannot use the gap types when some gaps are adjacent or overlap")
        nir1 <- asNormalIRanges(ranges)
    } else {
        nir1 <- toNormalIRanges(ranges)
    }
    new("MaskCollection", nir_list=list(nir1), width=width, active=TRUE, NAMES=mask_name)
}

read.agpMask <- function(file, width, seqname="?", gap.types=NULL, use.gap.types=FALSE)
{
    .read.MaskFromAgpOrGap("agp", file, width, seqname, gap.types, use.gap.types)
}

read.gapMask <- function(file, width, seqname="?", gap.types=NULL, use.gap.types=FALSE)
{
    .read.MaskFromAgpOrGap("gap", file, width, seqname, gap.types, use.gap.types)
}

read.liftMask <- function(file, seqname="?", width=NA)
{
    if (!isSingleString(seqname))
        stop("'seqname' must be a single string")
    if (!isSingleNumberOrNA(width))
        stop("'width' must be a single integer or 'NA'")
    if (!is.integer(width))
        width <- as.integer(width)
    ALL_COLS <- c(
        `offset`="integer",
        `xxxx`="NULL",  # not sure how to call this
        `width`="integer",
        `seqname`="character",
        `seqlen`="integer"
    )
    data <- read.table(file,
                       col.names=names(ALL_COLS),
                       colClasses=ALL_COLS,
                       check.names=FALSE)
    if (seqname == "?") {
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
    ALL_COLS <- c(
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
    ALL_COLS[!(names(ALL_COLS) %in% COLS)] <- "NULL"
    data <- read.table(file,
                       col.names=names(ALL_COLS),
                       colClasses=ALL_COLS,
                       skip=3,
                       check.names=FALSE)
    ranges <- IRanges(start=data$begin_in_query, end=data$end_in_query)
    if (use.IDs) {
        names(ranges) <- data$ID
        if (isNotStrictlySorted(start(ranges)))
            ranges <- ranges[order(start(ranges))]
        if (!isNormal(ranges))
            stop("cannot use the repeat IDs when some repeats are adjacent or overlap")
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
    ALL_COLS <- c(
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
    ALL_COLS[!(names(ALL_COLS) %in% COLS)] <- "NULL"
    data <- read.table(file,
                       col.names=names(ALL_COLS),
                       colClasses=ALL_COLS,
                       check.names=FALSE)
    ranges <- IRanges(start=data$chromStart+1, end=data$chromEnd)
    nir1 <- toNormalIRanges(ranges)
    #name1 <- "Tandem Repeats Finder [period<=12]"
    new("MaskCollection", nir_list=list(nir1), width=width, active=TRUE, NAMES="trf")
}

