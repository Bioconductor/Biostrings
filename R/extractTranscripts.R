### =========================================================================
### extractTranscripts() & related functions
### -------------------------------------------------------------------------


transcriptWidths <- function(exonStarts=list(), exonEnds=list())
{
    if (is(exonStarts, "RangesList")) {
        if (!identical(exonEnds, list()))
            stop("'exonEnds' cannot be specified ",
                 "when 'exonStarts' is a RangesList object")
        exonEnds <- end(exonStarts)
        exonStarts <- start(exonStarts)
    }
    if (!is.list(exonStarts)) {
        if (!is(exonStarts, "IntegerList"))
            stop("'exonStarts' must be a list of integer vectors ",
                 "or an IntegerList object")
        exonStarts <- as.list(exonStarts)
    }
    if (!is.list(exonEnds)) {
        if (!is(exonEnds, "IntegerList"))
            stop("'exonEnds' must be a list of integer vectors ",
                 "or an IntegerList object")
        exonEnds <- as.list(exonEnds)
    }
    if (length(exonStarts) != length(exonEnds))
        stop("'exonStarts', 'exonEnds' must have the same length")
    .Call("transcript_widths", exonStarts, exonEnds, PACKAGE="Biostrings")
}

extractTranscripts <- function(x, exonStarts=list(), exonEnds=list(),
                               strand=character(0),
                               reorder.exons.on.minus.strand=FALSE)
{
    if (!is(x, "DNAString")) {
        if (!is(x, "MaskedDNAString"))
            stop("'x' must be a DNAString object")
        masks(x) <- NULL
    }
    if (is(exonStarts, "RangesList")) {
        if (!identical(exonEnds, list()))
            stop("'exonEnds' cannot be specified ",
                 "when 'exonStarts' is a RangesList object")
        exonEnds <- end(exonStarts)
        exonStarts <- start(exonStarts)
    }
    if (!is.list(exonStarts)) {
        if (!is(exonStarts, "IntegerList"))
            stop("'exonStarts' must be a list of integer vectors ",
                 "or an IntegerList object")
        exonStarts <- as.list(exonStarts)
    }
    if (!is.list(exonEnds)) {
        if (!is(exonEnds, "IntegerList"))
            stop("'exonEnds' must be a list of integer vectors ",
                 "or an IntegerList object")
        exonEnds <- as.list(exonEnds)
    }
    if (is.factor(strand))
        strand <- as.vector(strand)
    if (!is.character(strand))
        stop("'strand' must be a character vector")
    if (length(exonStarts) != length(strand)
     || length(exonEnds) != length(strand))
        stop("'exonStarts', 'exonEnds' and 'strand' must have the same length")
    if (!isTRUEorFALSE(reorder.exons.on.minus.strand))
        stop("'reorder.exons.on.minus.strand' must be TRUE or FALSE")
    lkup <- getDNAComplementLookup()
    .Call("extract_transcripts",
          x, exonStarts, exonEnds, strand, reorder.exons.on.minus.strand, lkup,
          PACKAGE="Biostrings")
}

transcriptLocs2refLocs <- function(tlocs, exonStarts=list(), exonEnds=list(),
                                   strand=character(0),
                                   reorder.exons.on.minus.strand=FALSE)
{
    if (!is.list(tlocs)) {
        if (!is(tlocs, "IntegerList"))
            stop("'tlocs' must be a list of integer vectors ",
                 "or an IntegerList object")
        tlocs <- as.list(tlocs)
    }
    if (is(exonStarts, "RangesList")) {
        if (!identical(exonEnds, list()))
            stop("'exonEnds' cannot be specified ",
                 "when 'exonStarts' is a RangesList object")
        exonEnds <- end(exonStarts)
        exonStarts <- start(exonStarts)
    }
    if (!is.list(exonStarts)) {
        if (!is(exonStarts, "IntegerList"))
            stop("'exonStarts' must be a list of integer vectors ",
                 "or an IntegerList object")
        exonStarts <- as.list(exonStarts)
    }
    if (!is.list(exonEnds)) {
        if (!is(exonEnds, "IntegerList"))
            stop("'exonEnds' must be a list of integer vectors ",
                 "or an IntegerList object")
        exonEnds <- as.list(exonEnds)
    }
    if (is.factor(strand))
        strand <- as.vector(strand)
    if (!is.character(strand))
        stop("'strand' must be a character vector")
    if (length(tlocs) != length(strand)
     || length(exonStarts) != length(strand)
     || length(exonEnds) != length(strand))
        stop("'tlocs', 'exonStarts', 'exonEnds' and 'strand' ",
             "must have the same length")
    if (!isTRUEorFALSE(reorder.exons.on.minus.strand))
        stop("'reorder.exons.on.minus.strand' must be TRUE or FALSE")
    .Call("tlocs2rlocs",
          tlocs, exonStarts, exonEnds, strand, reorder.exons.on.minus.strand,
          PACKAGE="Biostrings")
}

