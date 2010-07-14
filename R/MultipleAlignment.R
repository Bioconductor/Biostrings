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

.read.ClustalWAln <-
function(filepath)
{
    con <- file(filepath)
    data <- readLines(con)
    close(con)  
    ## drop the header and empty lines
    data <- data[grep("^CLUSTAL", data, invert=TRUE)]
    data <- data[data != ""]
    ## get the index of the 1st line to be a complete blank line
    count <- grep("^(\\s|\\*)+$", data)[1] - 1L
    data <- data[grep("^(\\s|\\*)+$", data, invert=TRUE)]
    ## Therefore we shall gather and then drop the IDs
    ids <- gsub("(^gi\\S+)\\s+.+", "\\1", data)[seq_len(count)]
    data <- gsub("^gi\\S+\\s+", "", data)
    ## And also the positions from the end
    data <- gsub("\\s+\\d+$", "", data)
    ## And our count tells us how to split these up
    data <- split(data, factor(rep(seq_len(count), length(data)/count)))
    data <- unlist(lapply(data, paste, collapse=""))
    names(data) <- ids
    data
}

read.DNAMultipleAlignment <-
function(filepath, format=c("fasta", "clustalw"))
{
    if (!isSingleString(format)) 
        stop("'format' must be a single string")
    format <- match.arg(tolower(format), c("fasta", "clustalw"))
    data <-
      switch(format,
             "clustalw" = .read.ClustalWAln(filepath),
             read.DNAStringSet(filepath, format=format))
    DNAMultipleAlignment(data) 
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
        m <- as.matrix(unmasked(x), use.names=use.names)
        if (maskednrow(x) > 0)
            m <- m[- as.integer(rowmask(x)), , drop=FALSE]
        if (maskedncol(x) > 0)
            m <- m[, - as.integer(colmask(x)), drop=FALSE]
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
        m <- consensusMatrix(strings, as.prob=as.prob, baseOnly=baseOnly)
        if (maskedncol(x) > 0)
            m[, as.integer(colmask(x))] <- NA
        m
    }
)

setMethod("consensusString","MultipleAlignment",
    function(x, ...)
    {
        strings <- unmasked(x)
        if (maskednrow(x) > 0)
            strings <- strings[- as.integer(rowmask(x))]
        consensus <- consensusString(strings, ...)
        if (maskedncol(x) > 0) {
            consensus <- safeExplode(consensus)
            consensus[as.integer(colmask(x))] <- "#"
            consensus <- paste(consensus, collapse = "")
        }
        consensus
    }
)

setMethod("alphabetFrequency","MultipleAlignment",
    function(x, as.prob=FALSE, collapse=FALSE)
    {
        if (collapse) {
            callGeneric(as(x, sprinf("%sStringSet")), as.prob=as.prob,
                        collapse=TRUE)
        } else {
            strings <- unmasked(x)
            if (maskednrow(x) > 0)
                strings <- strings[- as.integer(rowmask(x))]
            m <- alphabetFrequency(strings, as.prob=as.prob)
            if (maskedncol(x) > 0)
                m[, as.integer(colmask(x))] <- NA
            m
        }
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
        cat(format(" seq", width = getOption("width") - iW - .namesW - 1),
            format("names", width=.namesW, justify="left"), sep="")
    } else {
        cat(" seq")
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
