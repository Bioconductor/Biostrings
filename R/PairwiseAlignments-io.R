### =========================================================================
### Input/output of PairwiseAlignments objects
### -------------------------------------------------------------------------
###
### Only output is supported at the moment.
###

.pre2postaligned <- function(pos, axset)
{
    if (!is(axset, "AlignedXStringSet0") || length(axset) != 1L)
        stop("'axset' must be an AlignedXStringSet0 object of length 1")
    ## .postaligned_gap_ranges() below sometimes needs to call
    ## .pre2postaligned() with a 'pos' that is 'end(axset@range) + 1L'.
    ## This happens when there is a gap at the end of the alignment like in
    ##   x <- DNAString("TCAACTTAACTT")
    ##   y <- DNAString("GGGCAACAACGGG")
    ##   pa <- pairwiseAlignment(x, y, type="global-local",
    ##                           gapOpening=-2, gapExtension=-1)
    ##   writePairwiseAlignments(pa)
    stopifnot(all(pos >= start(axset@range)),
              all(pos <= end(axset@range) + 1L))
    lkup <- integer(width(axset@range) + 1L)
    gap_ranges <- indel(axset)[[1L]]
    lkup[start(gap_ranges)] <- width(gap_ranges)
    lkup <- cumsum(lkup + 1L)
    lkup[pos - start(axset@range) + 1L]
}

.test.pre2postaligned <- function(axset)
{
    if (!is(axset, "AlignedXStringSet0") || length(axset) != 1L)
        stop("'axset' must be an AlignedXStringSet0 object of length 1")
    target <- subseq(unaligned(axset)[[1L]],
                     start(axset@range), end(axset@range))
    pos <- as.integer(axset@range)
    current <- aligned(axset)[[1L]][.pre2postaligned(pos, axset)]
    identical(as.character(target), as.character(current))
}

.postaligned_gap_ranges <- function(axset)
{
    if (!is(axset, "AlignedXStringSet0") || length(axset) != 1L)
        stop("'axset' must be an AlignedXStringSet0 object of length 1")
    gap_ranges <- indel(axset)[[1L]]
    prealigned_gap_start <- start(gap_ranges) + start(axset@range) - 1L
    postaligned_gap_start <- .pre2postaligned(prealigned_gap_start, axset) -
                             width(gap_ranges)
    #S4Vectors:::fancy_mseq(width(gap_ranges), postaligned_gap_start - 1L) 
    IRanges(postaligned_gap_start, width=width(gap_ranges))
}

.postaligned_match_ranges <- function(axset)
{
    if (!is(axset, "AlignedXStringSet0") || length(axset) != 1L)
        stop("'axset' must be an AlignedXStringSet0 object of length 1")
    aligned_axset <- aligned(axset) 
    postaligned_width <- width(aligned_axset)
    prealigned_mismatches <- mismatch(axset)[[1L]]
    is_in_range <- prealigned_mismatches >= start(axset@range) &
                   prealigned_mismatches <= end(axset@range)
    prealigned_mismatches <- prealigned_mismatches[is_in_range]
    postaligned_mismatches <- .pre2postaligned(prealigned_mismatches, axset)
    postaligned_mismatch_ranges <- as(postaligned_mismatches, "IRanges")
    postaligned_gap_ranges <- .postaligned_gap_ranges(axset)
    postaligned_mismatches_or_indels <- union(postaligned_mismatch_ranges,
                                              postaligned_gap_ranges)
    setdiff(IRanges(1L, postaligned_width), postaligned_mismatches_or_indels)
}

.makePipes <- function(x)
{
    if (!is(x, "PairwiseAlignments") || length(x) != 1L)
        stop("'x' must be a PairwiseAlignments object of length 1")
    x_pattern <- pattern(x)  # QualityAlignedXStringSet object
    x_subject <- subject(x)  # QualityAlignedXStringSet object
    postaligned_pattern <- aligned(x_pattern)[[1L]]
    postaligned_subject <- aligned(x_subject)[[1L]]
    postaligned_pattern_match_ranges <- .postaligned_match_ranges(x_pattern)
    postaligned_subject_match_ranges <- .postaligned_match_ranges(x_subject)
    postaligned_match_ranges <- intersect(postaligned_pattern_match_ranges,
                                          postaligned_subject_match_ranges)
    ## Sanity check:
    ii <- unlist(postaligned_match_ranges)
    if (!identical(as.character(postaligned_pattern[ii]),
                   as.character(postaligned_subject[ii])))
        stop("Biostrings internal error: mismatches and/or indels ",
             "reported in 'pattern(x)' and 'subject(x)' are inconsistent")
    tmp <- rep.int(" ", length(postaligned_pattern))
    tmp[ii] <- "|"
    BString(paste(tmp, collapse=""))
}

.makeXStringSetFrom2XStrings <- function(x1=NULL, x2=NULL)
{
    #if (is.null(x1) && is.null(x2))
    #    return(BStringSet(c("", "")))
    if (is.null(x1)) {
        x1 <- XString(seqtype(x2), "-")
        x1 <- rep.int(x1, length(x2))
    }
    if (is.null(x2)) {
        x2 <- XString(seqtype(x1), "-")
        x2 <- rep.int(x2, length(x1))
    }
    c(as(x1, "XStringSet"), as(x2, "XStringSet"))
}

.makeXStringFromSpaces <- function(nspace)
{
    rep.int(BString(" "), nspace)
}

.makePostalignedSeqs <- function(x, trim.global=FALSE)
{
    if (!is(x, "PairwiseAlignments") || length(x) != 1L)
        stop("'x' must be a PairwiseAlignments object of length 1")
    x_pattern <- pattern(x)  # QualityAlignedXStringSet object
    x_subject <- subject(x)  # QualityAlignedXStringSet object
    postaligned_pattern <- aligned(x_pattern)[[1L]]
    postaligned_subject <- aligned(x_subject)[[1L]]
    ans_seqs <- c(as(postaligned_pattern, "XStringSet"),
                  as(postaligned_subject, "XStringSet"))
    ## Sanity check:
    if (!isConstant(width(ans_seqs)))
        stop("Biostrings internal error: the 2 post-aligned sequences ",
             "must have the same length")
    unaligned_pattern <- unaligned(x_pattern)
    pattern_name <- names(unaligned(x_pattern))
    if (is.null(pattern_name))
        pattern_name <- ""
    unaligned_subject <- unaligned(x_subject)
    subject_name <- names(unaligned(x_subject))
    if (is.null(subject_name))
        subject_name <- ""
    ans_ranges <- c(x_pattern@range, x_subject@range)
    ans_pipes <- .makePipes(x)
    if (trim.global) {
        names(ans_seqs) <- c(pattern_name, subject_name)
        return(list(ans_seqs, ans_ranges, ans_pipes))
    }
    x_type <- type(x)
    is_pattern_global <- x_type %in% c("global", "global-local")
    is_subject_global <- x_type %in% c("global", "local-global")
    if (is_pattern_global) {
        unaligned_pattern <- unaligned_pattern[[1L]]
        start1 <- start(ans_ranges)[1L]
        if (start1 > 1L) {
            prefix1 <- subseq(unaligned_pattern, end=start1 - 1L)
            prefixes <- .makeXStringSetFrom2XStrings(x1=prefix1)
            ans_seqs <- xscat(prefixes, ans_seqs)
            prefix <- .makeXStringFromSpaces(length(prefix1))
            ans_pipes <- xscat(prefix, ans_pipes)
            start(ans_ranges)[1L] <- 1L
        }
        end1 <- end(ans_ranges)[1L]
        if (end1 < length(unaligned_pattern)) {
            suffix1 <- subseq(unaligned_pattern, start=end1 + 1L)
            suffixes <- .makeXStringSetFrom2XStrings(x1=suffix1)
            ans_seqs <- xscat(ans_seqs, suffixes)
            suffix <- .makeXStringFromSpaces(length(suffix1))
            ans_pipes <- xscat(ans_pipes, suffix)
            end(ans_ranges)[1L] <- length(unaligned_pattern)
        }
    }
    if (is_subject_global) {
        unaligned_subject <- unaligned_subject[[1L]]
        start2 <- start(ans_ranges)[2L]
        if (start2 > 1L) {
            prefix2 <- subseq(unaligned_subject, end=start2 - 1L)
            prefixes <- .makeXStringSetFrom2XStrings(x2=prefix2)
            ans_seqs <- xscat(prefixes, ans_seqs)
            prefix <- .makeXStringFromSpaces(length(prefix2))
            ans_pipes <- xscat(prefix, ans_pipes)
            start(ans_ranges)[2L] <- 1L
        }
        end2 <- end(ans_ranges)[2L]
        if (end2 < length(unaligned_subject)) {
            suffix2 <- subseq(unaligned_subject, start=end2 + 1L)
            suffixes <- .makeXStringSetFrom2XStrings(x2=suffix2)
            ans_seqs <- xscat(ans_seqs, suffixes)
            suffix <- .makeXStringFromSpaces(length(suffix2))
            ans_pipes <- xscat(ans_pipes, suffix)
            end(ans_ranges)[2L] <- length(unaligned_subject)
        }
    }
    names(ans_seqs) <- c(pattern_name, subject_name)
    list(ans_seqs, ans_ranges, ans_pipes)
}

.writePairHeader <- function(x, alignment.length, Identity, Similarity, Gaps, 
                             pattern.name="P1", subject.name="S1",
                             Matrix=NA, file="")
{
    if (!is(x, "PairwiseAlignments") || length(x) != 1L)
        stop("'x' must be a PairwiseAlignments object of length 1")
    if (!isSingleNumber(alignment.length))
        stop("'alignment.length' must be a single number")
    if (!is.integer(alignment.length))
        alignment.length <- as.integer(alignment.length)
    if (!isSingleStringOrNA(Matrix))
        stop("'Matrix' must be a single string or NA")

    Gap_penalty <- sprintf("%.1f", (x@gapOpening + x@gapExtension))
    Extend_penalty <- sprintf("%.1f", x@gapExtension)
    prettyPercentage <- function(ratio)
        sprintf("%.1f%%", round(ratio * 100, digits=1L))
    Identity_percentage <- prettyPercentage(Identity / alignment.length)
    Identity <- paste(format(Identity, width=7L), "/", alignment.length,
                      " (", Identity_percentage, ")", sep="")
    Similarity_percentage <- prettyPercentage(Similarity / alignment.length)
    Similarity <- paste(format(Similarity, width=5L), "/", alignment.length,
                        " (", Similarity_percentage, ")", sep="")
    Gaps_percentage <- prettyPercentage(Gaps / alignment.length)
    Gaps <- paste(format(Gaps, width=11L), "/", alignment.length,
                  " (", Gaps_percentage, ")", sep="")
    Score <- x@score

    cat("#=======================================\n", file=file)
    cat("#\n", file=file)
    cat("# Aligned_sequences: 2\n", file=file)
    cat("# 1: ", pattern.name, "\n", sep="", file=file)
    cat("# 2: ", subject.name, "\n", sep="", file=file)
    cat("# Matrix: ", Matrix, "\n", sep="", file=file)
    cat("# Gap_penalty: ", Gap_penalty, "\n", sep="", file=file)
    cat("# Extend_penalty: ", Extend_penalty, "\n", sep="", file=file)
    cat("#\n", file=file)
    cat("# Length: ", alignment.length, "\n", sep="", file=file)
    cat("# Identity: ", Identity, "\n", sep="", file=file)
    cat("# Similarity: ", Similarity, "\n", sep="", file=file)
    cat("# Gaps: ", Gaps, "\n", sep="", file=file)
    cat("# Score: ", Score, "\n", sep="", file=file)
    cat("#\n#\n", file=file)
    cat("#=======================================\n", file=file)
}

.writePairSequences <- function(top.string, bottom.string, middle.string,
                                top.name="P1", bottom.name="S1",
                                top.start=1L, bottom.start=1L,
                                block.width=50, file="")
{
    if (!isSingleNumber(block.width))
        stop("'block.width' must be a single number")
    if (!is.integer(block.width))
        block.width <- as.integer(block.width)

    alignment_length <- length(top.string)
    start_width <- max(nchar(as.character(top.start + alignment_length)),
                       nchar(as.character(bottom.start + alignment_length)))
    name_width <- max(20L - start_width - 1L,
                      nchar(top.name), nchar(bottom.name))
    nblock <- alignment_length %/% block.width
    if (alignment_length %% block.width != 0L)
        nblock <- nblock + 1L
    start1 <- top.start
    start3 <- bottom.start
    for (i in seq_len(nblock)) {
        to <- i * block.width
        from <- to - block.width + 1L
        if (to > alignment_length)
            to <- alignment_length
        string1 <- as.character(subseq(top.string, from, to))
        string2 <- as.character(subseq(middle.string, from, to))
        string3 <- as.character(subseq(bottom.string, from, to))
        if (i != 1L)
            cat("\n", file=file)
        ## 1st line
        cat(format(top.name, width=name_width), " ",
            format(start1, justify="right", width=start_width), " ",
            string1, file=file, sep="")
        end1 <- start1 + to - from - countPattern("-", string1)
        cat(format(end1, justify="right", width=7L),
            "\n", sep="", file=file)
        ## 2nd line
        cat(format("", width=name_width), " ",
            format("", width=start_width), " ",
            string2, "\n", file=file, sep="")
        ## 2rd line
        cat(format(bottom.name, width=name_width), " ",
            format(start3, justify="right", width=start_width), " ",
            string3, file=file, sep="")
        end3 <- start3 + to - from - countPattern("-", string3)
        cat(format(end3, justify="right", width=7L),
            "\n", sep="", file=file)
        start1 <- end1 + 1L
        start3 <- end3 + 1L
    }
}

writePairwiseAlignments <- function(x, file="", Matrix=NA, block.width=50)
{
    if (!is(x, "PairwiseAlignments"))
        stop("'x' must be a PairwiseAlignments object")
    if (isSingleString(file)) {
        if (file == "") {
            file <- stdout()
        } else if (substring(file, 1L, 1L) == "|") {
            file <- pipe(substring(file, 2L), "w")
            on.exit(close(file))
        } else {
            file <- file(file, "w")
            on.exit(close(file))
        }
    } else if (!is(file, "connection")) {
        stop("'file' must be a single string or a connection object")
    }
    pkgversion <- as.character(packageVersion("Biostrings"))
    Program <- paste("Biostrings (version ", pkgversion, "), ",
                     "a Bioconductor package", sep="")
    cat("########################################\n", file=file)
    cat("# Program: ", Program, "\n", sep="", file=file)
    cat("# Rundate: ", date(), "\n", sep="", file=file)
    cat("########################################\n", file=file)
    x_len <- length(x)
    if (x_len == 0L)
        warning("'x' is an empty PairwiseAlignments object ",
                "-> nothing to write")
    #else if (x_len >= 2L)
    #    warning("'x' contains more than 1 pairwise alignment")
    x_type <- type(x)
    is_pattern_global <- x_type %in% c("global", "global-local")
    is_subject_global <- x_type %in% c("global", "local-global")
    if (length(unaligned(subject(x))) != 1L) {
        bottom_name0 <- ""
    } else {
        bottom_name0 <- "S1"
    }
    for (i in seq_len(x_len)) {
        #if (i != 1L)
        #    cat("\n\n", file=file)
        xi <- x[i]
        postaligned_seqs <- .makePostalignedSeqs(xi)
        seqs <- postaligned_seqs[[1L]]
        ranges <- postaligned_seqs[[2L]]
        pipes <- postaligned_seqs[[3L]]
        name1 <- names(seqs)[1L]
        if (name1 == "")
            name1 <- paste("P", i, sep="")
        name2 <- names(seqs)[2L]
        if (name2 == "") {
            if (bottom_name0 == "") {
                name2 <- paste("S", i, sep="")
            } else {
                name2 <- bottom_name0
            }
        }
        Identity <- countPattern("|", pipes)
        Gaps <- sum(width(union(ranges(matchPattern("-", seqs[[1L]])),
                                ranges(matchPattern("-", seqs[[2L]])))))
        .writePairHeader(xi, length(seqs[[1L]]), Identity, NA, Gaps,
                         pattern.name=name1, subject.name=name2,
                         Matrix=Matrix, file=file)
        cat("\n", file=file)
        start1 <- start(ranges)[1L]
        start2 <- start(ranges[2L])
        .writePairSequences(seqs[[1L]], seqs[[2L]], pipes,
                            top.name=name1, bottom.name=name2,
                            top.start=start1, bottom.start=start2,
                            block.width=block.width, file=file)
        cat("\n\n", file=file)
    }
    cat("#---------------------------------------\n", file=file)
    cat("#---------------------------------------\n", file=file)
    invisible(NULL)
}

