### =========================================================================
### MultipleAlignment objects
### -------------------------------------------------------------------------
###

## virtual class for metadata
setClass("AlignmentMetadata",
    representation=representation("VIRTUAL")
)

## TODO: ... classes for each aligner

## virtual class for multiple alignments
setClass("MultipleAlignment",
    representation("VIRTUAL",
#                   metadata="AlignmentMetadata",
                   unmasked="XStringSet",
                   rowmask="NormalIRanges",
                   colmask="NormalIRanges")
)

## concrete classes for multiple alignments
setClass("DNAMultipleAlignment",
    contains="MultipleAlignment",
    representation=representation(unmasked="DNAStringSet")
)
setClass("RNAMultipleAlignment",
    contains="MultipleAlignment",
    representation=representation(unmasked="RNAStringSet")
)
setClass("AAMultipleAlignment",
    contains="MultipleAlignment",
    representation=representation(unmasked="AAStringSet")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.MultipleAlignment.unmasked <- function(object)
{
    if (length(unique(nchar(object@unmasked))) > 1)
        "alignments have an unequal number of characters"
    else
        NULL
}

.valid.MultipleAlignment.rowmask <- function(object)
{
    rng <- range(object@rowmask)
    if (length(rng) > 0 &&
        (start(rng) < 1 || end(rng) > length(object@unmasked)))
        "'rowmask' contains values outside of [1, nrow]"
    else
        NULL
}

.valid.MultipleAlignment.colmask <- function(object)
{
    rng <- range(object@colmask)
    if (length(rng) > 0 &&
        (start(rng) < 1 || end(rng) > nchar(object@unmasked)[1L]))
        "'colmask' contains values outside of [1, ncol]"
    else
        NULL
}

setValidity("MultipleAlignment",
    function(object)
    {
        problems <-
          c(.valid.MultipleAlignment.unmasked(object),
            .valid.MultipleAlignment.rowmask(object),
            .valid.MultipleAlignment.colmask(object))
        if (is.null(problems)) TRUE else problems
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor-like methods.
###

setMethod("unmasked", "MultipleAlignment", function(x) x@unmasked)

setMethod("rownames","MultipleAlignment", function(x) names(x@unmasked))
setReplaceMethod("rownames", "MultipleAlignment",
    function(x,value)
    {
        names(x@unmasked) <- value
        x
    }
)

setGeneric("rowmask", signature="x", function(x) standardGeneric("rowmask"))
setMethod("rowmask", "MultipleAlignment", function(x) x@rowmask)

setGeneric("rowmask<-", function(x,value) standardGeneric("rowmask<-"))
setReplaceMethod("rowmask", signature(x="MultipleAlignment", value="NULL"),
    function(x, value) initialize(x, rowmask = new("NormalIRanges"))
)
setReplaceMethod("rowmask",
    signature(x="MultipleAlignment", value="NormalIRanges"),
    function(x, value) initialize(x, rowmask = value)
)
setReplaceMethod("rowmask",
    signature(x="MultipleAlignment", value="IRanges"),
    function(x, value) initialize(x, rowmask = asNormalIRanges(value))
)

setGeneric("colmask", signature="x", function(x) standardGeneric("colmask"))
setMethod("colmask", "MultipleAlignment", function(x) x@colmask)

setGeneric("colmask<-", function(x,value) standardGeneric("colmask<-"))
setReplaceMethod("colmask", signature(x="MultipleAlignment", value="NULL"),
    function(x, value) initialize(x, colmask = new("NormalIRanges"))
)
setReplaceMethod("colmask",
    signature(x="MultipleAlignment", value="NormalIRanges"),
    function(x, value) initialize(x, colmask = value)
)
setReplaceMethod("colmask",
    signature(x="MultipleAlignment", value="IRanges"),
    function(x, value) initialize(x, colmask = asNormalIRanges(value))
)

setMethod("nrow","MultipleAlignment", function(x) length(x@unmasked))
setMethod("ncol","MultipleAlignment",
    function(x) ifelse(nrow(x) == 0, 0L, nchar(x@unmasked[[1L]]))
)
setMethod("dim","MultipleAlignment", function(x) c(nrow(x), ncol(x)))

setGeneric("maskednrow", signature="x",
    function(x) standardGeneric("maskednrow")
)
setMethod("maskednrow", "MultipleAlignment",
    function(x) sum(width(rowmask(x)))
)

setGeneric("maskedncol", signature="x",
    function(x) standardGeneric("maskedncol")
)
setMethod("maskedncol", "MultipleAlignment",
    function(x) sum(width(colmask(x)))
)

setGeneric("maskeddim", signature="x",
    function(x) standardGeneric("maskeddim")
)
setMethod("maskeddim", "MultipleAlignment",
    function(x) c(maskednrow(x), c(maskedncol(x)))
)

setMethod("maskedratio", "MultipleAlignment",
    function(x) maskeddim(x) / dim(x)
)

setMethod("nchar", "MultipleAlignment",
    function(x) ncol(x) - maskedncol(x)
)

setMethod("xsbasetype", "MultipleAlignment",
    function(x) xsbasetype(unmasked(x))
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructors.
###

DNAMultipleAlignment <-
function(x=character(), start=NA, end=NA, width=NA, use.names=TRUE)
{
    new("DNAMultipleAlignment",
        unmasked=
        DNAStringSet(x=x, start=start, end=end, width=width,
                     use.names=use.names))
}

RNAMultipleAlignment <-
function(x=character(), start=NA, end=NA, width=NA, use.names=TRUE)
{
    new("RNAMultipleAlignment",
        unmasked=
        RNAStringSet(x=x, start=start, end=end, width=width,
                     use.names=use.names))
}

AAMultipleAlignment <-
function(x=character(), start=NA, end=NA, width=NA, use.names=TRUE)
{
    new("AAMultipleAlignment",
        unmasked=
        AAStringSet(x=x, start=start, end=end, width=width,
                    use.names=use.names))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Read function.
###

.read.MultipleAlignment.splitRows <-
function(rows, markupPattern)
{
    markupLines <- grep(markupPattern, rows)
    alnLines <- gaps(as(markupLines, "IRanges"), start=1, end=length(rows))
    nseq <- unique(width(alnLines))
    if (length(nseq) != 1)
        stop("missing alignment rows")
    rows <- seqselect(rows, alnLines)
    spaces <- regexpr("\\s+", rows)
    ids <- substr(rows, 1L, spaces - 1L)
    nsplits <- length(rows) %/% nseq
    if (!identical(ids, rep.int(head(ids, nseq), nsplits)))
        stop("alignment rows out of order")
    alns <- substr(rows, spaces + attr(spaces, "match.length"), nchar(rows))
    structure(do.call(paste,
                      c(split(alns, rep(seq_len(nsplits), each=nseq)), sep="")),
              names = head(ids, nseq))
}

.read.Stockholm <-
function(filepath)
{
    rows <- scan(filepath, what = "", sep = "\n", strip.white = TRUE,
                 quiet = TRUE, blank.lines.skip = FALSE)
    if (length(rows) < 3 ||
        !identical(grep("^# STOCKHOLM", rows[1L]), 1L))
        stop("invalid Stockholm file")
    chartr(".", "-",
           .read.MultipleAlignment.splitRows(rows, "(^\\s*|^#.*|^//\\s*)$"))
}

.read.ClustalAln <-
function(filepath)
{
    rows <- scan(filepath, what = "", sep = "\n", strip.white = TRUE,
                 quiet = TRUE, blank.lines.skip = FALSE)
    if (length(rows) < 3 ||
        !identical(grep("^CLUSTAL", rows[1L]), 1L) ||
        !identical(rows[2:3], c("","")))
        stop("invalid Clustal aln file")
    rows <- tail(rows, -3)
    rows <- sub("^(\\S+\\s+\\S+)\\s*\\d*$", "\\1", rows)
    .read.MultipleAlignment.splitRows(rows, "^(\\s|\\*|:|\\.)*$")
}

.read.MultipleAlignment <-
function(filepath, format)
{
    if (missing(format)) {
        ext <- tolower(sub(".*\\.([^.]*)$", "\\1", filepath))
        format <- switch(ext, "sto" = "stockholm", "aln" = "clustal", "fasta")
    } else {
        format <- match.arg(tolower(format), c("fasta", "stockholm", "clustal"))
    }
    switch(format,
           "stockholm" = .read.Stockholm(filepath),
           "clustal" = .read.ClustalAln(filepath),
           read.DNAStringSet(filepath, format=format))
}

read.DNAMultipleAlignment <-
function(filepath, format)
{
    DNAMultipleAlignment(.read.MultipleAlignment(filepath, format))
}

read.RNAMultipleAlignment <-
function(filepath, format)
{
    RNAMultipleAlignment(.read.MultipleAlignment(filepath, format))
}

read.AAMultipleAlignment <-
function(filepath, format)
{
    AAMultipleAlignment(.read.MultipleAlignment(filepath, format))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

setAs("MultipleAlignment", "DNAStringSet",
    function(from) DNAStringSet(as.character(from))
)
setAs("MultipleAlignment", "RNAStringSet",
    function(from) RNAStringSet(as.character(from))
)
setAs("MultipleAlignment", "AAStringSet",
    function(from) AAStringSet(as.character(from))
)
setAs("MultipleAlignment", "BStringSet",
    function(from) BStringSet(as.character(from))
)

setAs("character", "DNAMultipleAlignment",
      function(from) DNAMultipleAlignment(from)
)
setAs("character", "RNAMultipleAlignment",
      function(from) RNAMultipleAlignment(from)
)
setAs("character", "AAMultipleAlignment",
      function(from) AAMultipleAlignment(from)
)

setMethod("as.character", "MultipleAlignment",
    function(x, use.names=TRUE)
        apply(as.matrix(x, use.names=use.names), 1, paste, collapse="")
)

setMethod("as.matrix", "MultipleAlignment",
    function(x, use.names=TRUE)
    {
        if (nrow(x) == 0) {
            m <- matrix(character(0), nrow=0, ncol=0)
        } else {
            m <- as.matrix(unmasked(x), use.names=use.names)
            if (maskednrow(x) > 0)
                m <- m[- as.integer(rowmask(x)), , drop=FALSE]
            if (maskedncol(x) > 0)
                m <- m[, - as.integer(colmask(x)), drop=FALSE]
        }
        m
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Utilities.
###

setMethod("consensusMatrix","MultipleAlignment",
    function(x, as.prob=FALSE, baseOnly=FALSE)
    {
        strings <- unmasked(x)
        if (maskednrow(x) > 0)
            strings <- strings[- as.integer(rowmask(x))]
        m <- callGeneric(strings, as.prob=as.prob, baseOnly=baseOnly,
                         width=ncol(x))
        if (maskedncol(x) > 0)
            m[, as.integer(colmask(x))] <- NA
        m
    }
)

setMethod("consensusString","MultipleAlignment",
    function(x, ambiguityMap, threshold, codes)
    {
        strings <- unmasked(x)
        if (maskednrow(x) > 0)
            strings <- strings[- as.integer(rowmask(x))]
        cmat <- consensusMatrix(strings, width=ncol(x))[codes, , drop=FALSE]
        col_sums <- colSums(cmat)
        col_sums[col_sums == 0] <- 1  # to avoid division by 0
        cmat <- cmat / rep(col_sums, each=nrow(cmat))
        consensus <-
          consensusString(cmat, ambiguityMap=ambiguityMap, threshold=threshold)
        if (ncol(x) > 0 && length(consensus) == 0) {
            consensus <- paste(rep.int("#", ncol(x)), collapse="")
        } else if (maskedncol(x) > 0) {
            consensus <- safeExplode(consensus)
            consensus[as.integer(colmask(x))] <- "#"
            consensus <- paste(consensus, collapse = "")
        }
        consensus
    }
)

setMethod("consensusString","DNAMultipleAlignment",
    function(x, ambiguityMap=IUPAC_CODE_MAP, threshold=0.25)
    {
        callNextMethod(x, ambiguityMap=ambiguityMap, threshold=threshold,
                       codes=names(IUPAC_CODE_MAP))
    }
)

setMethod("consensusString","RNAMultipleAlignment",
    function(x,
            ambiguityMap=as.character(RNAStringSet(DNAStringSet(IUPAC_CODE_MAP))),
            threshold=0.25)
    {
        callNextMethod(x, ambiguityMap=ambiguityMap, threshold=threshold,
                       codes=
                       as.character(RNAStringSet(DNAStringSet(names(IUPAC_CODE_MAP)))))
    }
)

setMethod("consensusString","AAMultipleAlignment",
    function(x, ambiguityMap="?", threshold=0.5)
    {
        callNextMethod(x, ambiguityMap=ambiguityMap, threshold=threshold,
                       codes=names(AMINO_ACID_CODE))
    }
)

setMethod("alphabetFrequency","MultipleAlignment",
    function(x, as.prob=FALSE, collapse=FALSE)
    {
        if (collapse) {
            strings <- as(x, sprintf("%sStringSet", xsbasetype(x)))
            m <- callGeneric(strings, as.prob=as.prob, collapse=TRUE)
        } else {
            rmasks <- rowmask(x)
            if (length(rmasks) > 0)
                rowmask(x) <- NULL
            strings <- as(x, sprintf("%sStringSet", xsbasetype(x)))
            m <- callGeneric(strings, as.prob=as.prob)
            if (length(rmasks) > 0)
                m[as.integer(rmasks),] <- NA
        }
        m
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method.
###

.MultipleAlignment.show_frame_header <-
function (iW, with.names) 
{
    cat(format("", width = iW + 1), sep="")
    if (with.names) {
        cat(format(" aln", width = getOption("width") - iW - .namesW - 1),
            format("names", width=.namesW, justify="left"), sep="")
    } else {
        cat(" aln")
    }
    cat("\n")
}

.MultipleAlignment.show_frame_line <-
function (x, i, iW)
{
    snippetWidth <- getOption("width") - 2 - iW
    if (!is.null(names(x))) 
        snippetWidth <- snippetWidth - .namesW - 1
    seq_snippet <- toSeqSnippet(x[[i]], snippetWidth)
    if (!is.null(names(x))) 
        seq_snippet <- format(seq_snippet, width=snippetWidth)
    cat(format(paste("[", i, "]", sep=""), width=iW, justify="right"),
        " ", seq_snippet, sep="")
    if (!is.null(names(x))) {
        snippet_name <- names(x)[i]
        if (is.na(snippet_name))
            snippet_name <- "<NA>"
        else if (nchar(snippet_name) > .namesW)
            snippet_name <-
              paste(substr(snippet_name, 1, .namesW - 3), "...", sep="")
        cat(" ", snippet_name, sep="")
    }
    cat("\n")
}

.MultipleAlignment.show_frame <-
function (x, half_nrow=9L)
{
    lx <- length(x)
    iW <- nchar(as.character(lx)) + 2
    .MultipleAlignment.show_frame_header(iW, !is.null(names(x)))
    if (lx <= 2 * half_nrow + 1) {
        for (i in seq_len(lx))
            .MultipleAlignment.show_frame_line(x, i, iW)
    } else {
        for (i in 1:half_nrow)
            .MultipleAlignment.show_frame_line(x, i, iW)
        cat(format("...", width=iW, justify="right"), "...\n")
        for (i in (lx - half_nrow + 1L):lx)
            .MultipleAlignment.show_frame_line(x, i, iW)
    }
}

setMethod("show", "MultipleAlignment",
    function(object)
    {
        nr <- nrow(object)
        nc <- ncol(object)
        cat(class(object), " with ", nr, ifelse(nr == 1, " row and ", 
            " rows and "), nc, ifelse(nc == 1, " column\n", " columns\n"), 
            sep = "")
        if (nr != 0) {
            strings <- unmasked(object)
            mdim <- maskeddim(object)
            if (sum(mdim) > 0) {
                if (mdim[1] > 0) {
                    strings <- BStringSet(strings)
                    maskStrings <-
                      rep(BStringSet(paste(rep.int("#", ncol(object)),
                                           collapse="")), mdim[1])
                    i <- as.integer(rowmask(object))
                    if (!is.null(rownames(object)))
                        names(maskStrings) <- rownames(object)[i]
                    strings[i] <- maskStrings
                }
                if (mdim[2] > 0) {
                    strings <- as.matrix(strings)
                    strings[, as.integer(colmask(object))] <- "#"
                    strings <- BStringSet(apply(strings, 1, paste, collapse=""))
                }
            }
            .MultipleAlignment.show_frame(strings)
        }
    }
)
