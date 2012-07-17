###

if (FALSE) {
  pattern <- DNAStringSet(c(p1="ACCA", p2="ACGCA", p3="ACGGCA", p4="ACGAGCA", p5="ACGTTGCA"))
  subject <- DNAStringSet(c(s1="TTTTACGTGCATTTTACGCA"))
  x <- pairwiseAlignment(pattern, subject, type="global-local")
}

.pre2postaligned <- function(pos, axset)
{
    if (!is(axset, "AlignedXStringSet0") || length(axset) != 1L)
        stop("'axset' must be an AlignedXStringSet0 object of length 1")
    if (any(pos < start(axset@range)) || any(pos > end(axset@range)))
        stop("'pos' must be >= 'start(axset@range)' and <= 'end(axset@range)'")
    lkup <- integer(width(axset@range))
    gap_ranges <- indel(axset)[[1L]]
    lkup[start(gap_ranges)] <- width(gap_ranges)
    lkup <- cumsum(lkup + 1L)
    lkup[pos - start(axset@range) + 1L]
}

.test.pre2postaligned <- function(axset)
{
    if (!is(axset, "AlignedXStringSet0") || length(axset) != 1L)
        stop("'axset' must be an AlignedXStringSet0 object of length 1")
    target <- subseq(unaligned(axset)[[1L]], start(axset@range), end(axset@range))
    current <- aligned(axset)[[1L]][.pre2postaligned(as.integer(axset@range), axset)]
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
    #IRanges:::fancy_mseq(width(gap_ranges), postaligned_gap_start - 1L) 
    IRanges(postaligned_gap_start, width=width(gap_ranges))
}

.postaligned_match_ranges <- function(axset)
{
    if (!is(axset, "AlignedXStringSet0") || length(axset) != 1L)
        stop("'axset' must be an AlignedXStringSet0 object of length 1")
    aligned_axset <- aligned(axset) 
    postaligned_width <- width(aligned_axset)
    postaligned_mismatches <- .pre2postaligned(mismatch(axset)[[1L]], axset)
    postaligned_mismatch_ranges <- as(postaligned_mismatches, "IRanges")
    postaligned_gap_ranges <- .postaligned_gap_ranges(axset)
    postaligned_mismatches_or_indels <- union(postaligned_mismatch_ranges,
                                              postaligned_gap_ranges)
    setdiff(IRanges(1L, postaligned_width), postaligned_mismatches_or_indels)
}

.onePairTo3Strings <- function(x)
{
    if (!is(x, "PairwiseAlignedXStringSet") || length(x) != 1L)
        stop("'x' must be a PairwiseAlignedXStringSet object of length 1")
    aligned_pattern <- pattern(x)  # QualityAlignedXStringSet object
    aligned_subject <- subject(x)  # QualityAlignedXStringSet object
    postaligned_pattern <- aligned(aligned_pattern)[[1L]]
    postaligned_subject <- aligned(aligned_subject)[[1L]]
    postaligned_pattern_match_ranges <-
            .postaligned_match_ranges(aligned_pattern)
    postaligned_subject_match_ranges <-
            .postaligned_match_ranges(aligned_subject)
    postaligned_match_ranges <- intersect(postaligned_pattern_match_ranges,
                                          postaligned_subject_match_ranges)
    ## Sanity check:
    ii <- unlist(postaligned_match_ranges)
    if (!identical(as.character(postaligned_pattern[ii]),
                   as.character(postaligned_subject[ii])))
        stop("Biostrings internal error: mismatches and/or indels ",
             "reported in 'pattern(x)' and 'subject(x)' are inconsistent")
    top_string <- as.character(postaligned_pattern)
    bottom_string <- as.character(postaligned_subject)
    middle_string <- rep.int(" ", nchar(top_string))
    middle_string[ii] <- "|"
    middle_string <- paste(middle_string, collapse="")
    ans <- c(top_string, middle_string, bottom_string)
    top_name <- names(unaligned(aligned_pattern))
    if (is.null(top_name))
        top_name <- ""
    bottom_name <- names(unaligned(aligned_subject))
    if (is.null(bottom_name))
        bottom_name <- ""
    names(ans) <- c(top_name, "", bottom_name)
    ans
}

.writePairHeader <- function(x, alignment.length, Identity, Similarity, Gaps, 
                             pattern.name="P1", subject.name="S1",
                             Matrix=NA, file="")
{
    if (!is(x, "PairwiseAlignedXStringSet") || length(x) != 1L)
        stop("'x' must be a PairwiseAlignedXStringSet object of length 1")

    Gap_penalty <- sprintf("%.1f", - (x@gapOpening + x@gapExtension))
    Extend_penalty <- sprintf("%.1f", - x@gapExtension)
    Identity_percentage <- round(100 * Identity / alignment.length, digits=1L)
    Identity <- paste(format(Identity, width=7L), "/", alignment.length, " ",
                      "(", sprintf("%.1f", Identity_percentage), "%)", sep="")
    Similarity_percentage <- round(100 * Similarity / alignment.length, digits=1L)
    Similarity <- paste(format(Similarity, width=5L), "/", alignment.length, " ",
                        "(", sprintf("%.1f", Similarity_percentage), "%)", sep="")
    Gaps_percentage <- round(100 * Gaps / alignment.length, digits=1L)
    Gaps <- paste(format(Gaps, width=11L), "/", alignment.length, " ",
                  "(", sprintf("%.1f", Gaps_percentage), "%)", sep="")
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
    cat("\n\n", file=file)
    cat("#=======================================\n", file=file)
}

.writePairSequences <- function(top.string, bottom.string, middle.string,
                                top.name="P1", bottom.name="S1",
                                top.start=1L, bottom.start=1L,
                                block.width=50L, file="")
{
    alignment_length <- nchar(top.string)
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
        string1 <- substr(top.string, from, to)
        string2 <- substr(middle.string, from, to)
        string3 <- substr(bottom.string, from, to)
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

writePairwiseAlignedXStringSet <- function(x, file="", block.width=50L)
{
    if (!is(x, "PairwiseAlignedXStringSet"))
        stop("'x' must be a PairwiseAlignedXStringSet object")
    x_len <- length(x)
    if (x_len == 0L)
        warning("'x' is an empty PairwiseAlignedXStringSet object ",
                "-> nothing to write")
    else if (x_len >= 2L)
        warning("'x' contains more than 1 pairwise alignment")
    if (length(unaligned(subject(x))) != 1L) {
        bottom_name0 <- ""
    } else {
        bottom_name0 <- "S1"
    }
    for (i in seq_len(x_len)) {
        if (i != 1L)
            cat("\n\n", file=file)
        xi <- x[i]
        strings <- .onePairTo3Strings(xi)
        string_names <- names(strings)
        top_name <- string_names[1L]
        if (top_name == "")
            top_name <- paste("P", i, sep="")
        bottom_name <- string_names[3L]
        if (bottom_name == "") {
            if (bottom_name0 == "") {
                bottom_name <- paste("S", i, sep="")
            } else {
                bottom_name <- bottom_name0
            }
        }
        Identity <- countPattern("|", strings[2L])
        Gaps <- sum(width(union(ranges(matchPattern("-", strings[1L])),
                                ranges(matchPattern("-", strings[3L])))))
        .writePairHeader(xi, nchar(strings[1L]), Identity, NA, Gaps,
                         pattern.name=top_name, subject.name=bottom_name,
                         file=file)
        aligned_pattern <- pattern(xi)
        aligned_subject <- subject(xi)
        top_start <- start(aligned_pattern@range)
        bottom_start <- start(aligned_subject@range)
        cat("\n", file=file)
        .writePairSequences(strings[1L], strings[3L], strings[2L],
                            top.name=top_name, bottom.name=bottom_name,
                            top.start=top_start, bottom.start=bottom_start,
                            block.width=block.width, file=file)
        cat("\n\n", file=file)
        cat("#---------------------------------------\n", file=file)
        cat("#---------------------------------------\n", file=file)
    }
    invisible(NULL)
}

