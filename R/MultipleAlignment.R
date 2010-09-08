### =========================================================================
### MultipleAlignment objects
### -------------------------------------------------------------------------
###

## virtual class for multiple alignments
setClass("MultipleAlignment",
    representation("VIRTUAL",
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

.appendMask <- function(mask, append, value){
  append <- match.arg(append, choices=c("union", "replace", "intersect"))
  switch(append,
         "union" = return(union(mask, value)),
         "replace"  = return(value),
         "intersect"  = return(intersect(mask, value)))
}
.setMask <- function(mask, append, invert, length, value){
  if (!isTRUEorFALSE(invert))
    stop("'invert' must be TRUE or FALSE")
  if(invert==TRUE){
    ## 1st invert value using gaps()
    value <- gaps(value, start=1, end=length)     
    value <- .appendMask(mask=mask, append=append, value=value)
  }
  if(invert==FALSE){
    value <- .appendMask(mask=mask, append=append, value=value)
  } 
  value
}
setGeneric("rowmask<-", signature=c("x", "value"),
    function(x, append="union", invert=FALSE, value)
           standardGeneric("rowmask<-")
)
##HOW does the following call .setMask???
setReplaceMethod("rowmask", signature(x="MultipleAlignment", value="NULL"),
    function(x, append="replace", invert=FALSE, value)
        callGeneric(x, append=append, invert=invert,
                    value=new("NormalIRanges"))
)
setReplaceMethod("rowmask",
    signature(x="MultipleAlignment", value="NormalIRanges"),
    function(x, append="union", invert=FALSE, value)
    {
      value <- .setMask(mask=rowmask(x), append=append, invert=invert,
                        length=dim(x)[1], value=value)
      initialize(x, rowmask = value)
    }
)
setReplaceMethod("rowmask",
    signature(x="MultipleAlignment", value="ANY"),
    function(x, append="union", invert=FALSE, value)
        callGeneric(x, append=append, invert=invert,
                    value=as(value, "NormalIRanges"))
)

setGeneric("colmask", signature="x", function(x) standardGeneric("colmask"))
setMethod("colmask", "MultipleAlignment", function(x) x@colmask)

setGeneric("colmask<-", signature=c("x", "value"),
    function(x, append="union", invert=FALSE, value)
           standardGeneric("colmask<-")
)
setReplaceMethod("colmask", signature(x="MultipleAlignment", value="NULL"),
    function(x, append="replace", invert=FALSE, value)
        callGeneric(x, append=append, invert=invert,
                    value=new("NormalIRanges"))
)
setReplaceMethod("colmask",
    signature(x="MultipleAlignment", value="NormalIRanges"),
    function(x, append="union", invert=FALSE, value)
    {
      value <- .setMask(mask=colmask(x), append=append, invert=invert,
                        length=dim(x)[2], value=value)
      initialize(x, colmask = value)
    }
)
setReplaceMethod("colmask",
    signature(x="MultipleAlignment", value="ANY"),
    function(x, append="union",invert=FALSE, value)
        callGeneric(x, append=append, invert=invert,
                    value=as(value, "NormalIRanges"))
)

setMethod("maskMotif", signature(x="MultipleAlignment", motif="ANY"),
    function(x, motif, min.block.width=1, ...)
    {
        cmask <- colmask(x)
        if (length(colmask(x)) > 0)
            colmask(x) <- NULL
        string <- consensusString(x)
        if (length(string) == 1) {
            string <- gsub("#", "+", string)
            string <- as(string, xsbaseclass(unmasked(x)))
            maskedString <-
              callGeneric(string, motif, min.block.width=min.block.width, ...)
            newmask <- nir_list(masks(maskedString))[[1L]]
            colmask(x) <- union(newmask, cmask)
        }
        x
    }
)

setGeneric("maskGaps", signature="x",
    function(x, ...) standardGeneric("maskGaps")
)
setMethod("maskGaps", "MultipleAlignment",
    function(x, min.fraction=0.5, min.block.width=3)
    {
        if (!isSingleNumber(min.fraction) || min.fraction < 0 ||
            min.fraction > 1)
            stop("'min.fraction' must be a number in [0, 1]")
        if (!isSingleNumber(min.block.width) || min.block.width <= 0)
            stop("'min.block.width' must be a positive integer")
        if (!is.integer(min.block.width))
            min.block.width <- as.integer(min.block.width)

        cmask <- colmask(x)
        if (length(colmask(x)) > 0)
            colmask(x) <- NULL
        m <- consensusMatrix(x)
        newmask <- (m["-",] / colSums(m)) >= min.fraction
        newmask[is.na(newmask)] <- FALSE
        newmask <- as(newmask, "NormalIRanges")
        newmask <- newmask[width(newmask) >= min.block.width]
        colmask(x) <- union(newmask, cmask)
        x
    }
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
function(x=character(), start=NA, end=NA, width=NA, use.names=TRUE,
         rowmask=NULL, colmask=NULL)
{
    new("DNAMultipleAlignment",
        unmasked=
        DNAStringSet(x=x, start=start, end=end, width=width,
                     use.names=use.names),
        rowmask=rowmask,
        colmask=colmask)
}

RNAMultipleAlignment <-
function(x=character(), start=NA, end=NA, width=NA, use.names=TRUE,
         rowmask=NULL, colmask=NULL)
{
    new("RNAMultipleAlignment",
        unmasked=
        RNAStringSet(x=x, start=start, end=end, width=width,
                     use.names=use.names),
        rowmask=rowmask,
        colmask=colmask)
}

AAMultipleAlignment <-
function(x=character(), start=NA, end=NA, width=NA, use.names=TRUE,
         rowmask=NULL, colmask=NULL)
{
    new("AAMultipleAlignment",
        unmasked=
        AAStringSet(x=x, start=start, end=end, width=width,
                    use.names=use.names),
        rowmask=rowmask,
        colmask=colmask)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Read function.
###

## markupPattern specifies which lines to skip
.read.MultipleAlignment.splitRows <-
function(rows, markupPattern)
{
    markupLines <- grep(markupPattern, rows, perl=TRUE)
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

## In order to recycle .read.MultipleAlignment.splitRows().
## I need to have the names on each row.
.read.PhylipAln <- 
function(filepath, maskGen=FALSE)
{
    rows <- scan(filepath, what = "", sep = "\n", strip.white = TRUE,
                 quiet = TRUE, blank.lines.skip = FALSE)
    if (length(rows) < 1 ||
        !identical(grep("^\\d+?\\s\\d+?", rows[1L]), 1L))
        stop("invalid Phylip file")
    ##(mask+num rows + blank line) 
    nameLength = as.numeric(sub("(\\d+).*$","\\1", rows[1])) +1 
    rows <- tail(rows, -1)
    names = character()
    names[nameLength] = "" ## an empty string is ALWAYS the last "name"
    offset = 0L
    for(i in seq_len(length(rows))){
      if(i<nameLength){
        rows[i] = sub("(^\\S+)\\s+(\\S+)", "\\1\\|\\2", rows[i])
        rows[i] = gsub("\\s", "", rows[i])
        rows[i] = sub("\\|", " ", rows[i])
        names[i] = sub("(\\S+).*$","\\1",rows[i])
      }else{        
        rows[i] = gsub("\\s", "", rows[i])
        rows[i] = paste(names[i %% nameLength], rows[i])
      }
    }
    rows <- c(" ",rows)
    if(maskGen==FALSE){ ## filter out the Mask values OR blank lines
       .read.MultipleAlignment.splitRows(rows, "^(Mask|\\s)")
    }else{## only retrieve the Mask values
      if(length(grep("^(?!Mask)",rows, perl=TRUE))==length(rows)){
        return(as(IRanges(),"NormalIRanges"))
      }else{
        msk <- .read.MultipleAlignment.splitRows(rows, "^(?!Mask)")
        ## THEN cast them to be a NormalIRanges object.
        splt = strsplit(msk,"") ## split up all chars
        names(splt) <- NULL ## drop the name
        splt = unlist(splt) ## THEN unlist
        lsplt = as.logical(as.numeric(splt)) ## NOW you can get a logical
        return(gaps(as(lsplt,"NormalIRanges"))) ## gaps() inverts the mask
      }
    }
}

.checkFormat <- function(format){
    if (missing(format)) {
        ext <- tolower(sub(".*\\.([^.]*)$", "\\1", filepath))
        format <- switch(ext, "sto" = "stockholm", "aln" = "clustal", "fasta")
    } else {
        format <- match.arg(tolower(format), c("fasta", "stockholm", "clustal",
                                               "phylip"))
    }
    format
}

.read.MultipleAlignment <-
function(filepath, format)
{
    format <- .checkFormat(format)
    switch(format,
           "stockholm" = .read.Stockholm(filepath),
           "clustal" = .read.ClustalAln(filepath),
           "phylip" = .read.PhylipAln(filepath),
           read.DNAStringSet(filepath, format=format))
    ##fasta uses read.DNAStringSet (default)
    ##TODO: BUGs with stockholm??
}

.read.MultipleMask <-
function(filepath, format)
{
    format <- .checkFormat(format)
    switch(format,
           "stockholm" = as(IRanges(),"NormalIRanges"),
           "clustal" = as(IRanges(),"NormalIRanges"),
           "phylip" = .read.PhylipAln(filepath, maskGen=TRUE),
           as(IRanges(),"NormalIRanges"))
}

## TODO: to implement Phylip support:
## 1) modify constructor of each MultipleAlignment object so that they can also take a colmask argument at construction and apply it.
## 2) write .read.MultipleMask(filepath,format) to extract mask or return NULL
## 3) add a boolean maskGen parameter to .read.PhylipAln()

read.DNAMultipleAlignment <-
function(filepath, format)
{
    DNAMultipleAlignment(.read.MultipleAlignment(filepath, format),
                         rowmask=as(IRanges(),"NormalIRanges"),
                         colmask=.read.MultipleMask(filepath,format))
}

read.RNAMultipleAlignment <-
function(filepath, format)
{
    RNAMultipleAlignment(.read.MultipleAlignment(filepath, format),
                         rowmask=as(IRanges(),"NormalIRanges"),
                         colmask=.read.MultipleMask(filepath,format))
}

read.AAMultipleAlignment <-
function(filepath, format)
{
    AAMultipleAlignment(.read.MultipleAlignment(filepath, format),
                        rowmask=as(IRanges(),"NormalIRanges"),
                        colmask=.read.MultipleMask(filepath,format))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Write functions.
###

write.Phylip <- function(x, file){
  if(inherits(origMAlign, "MultipleAlignment")){
    ## 1st, we need to capture the colmask as a vector that can be included.
    msk = colmask(x)
    colmask(x) <- NULL
    ## Then massage this to be a character vector
    ch = as.character(x)
 
    ## Then we have to split it up, appending the names at the beginning.
    writeLines(ch, file)
  }
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
        consensus <- rep.int("#", ncol(x))
        if (ncol(x) > 0) {
            cmat <- consensusMatrix(x, width=ncol(x))
            ngaps <- cmat["-",]
            cmat <- cmat[codes, , drop=FALSE]
            colsums <- colSums(cmat, na.rm=TRUE)
            nzsum <- which(colsums > 1e-6)
            if (length(nzsum) > 0) {
                colsums[- nzsum] <- 1  # to avoid division by 0
                gaplocs <- which(as.logical(ngaps > colsums))
                nzsum <- setdiff(nzsum, gaplocs)
                cmat <- cmat / rep(colsums, each=nrow(cmat))
                consensus[gaplocs] <- "-"
                consensus[nzsum] <-
                  safeExplode(consensusString(cmat[, nzsum, drop=FALSE],
                                              ambiguityMap=ambiguityMap,
                                              threshold=threshold))
            }
            consensus <- paste(consensus, collapse="")
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
            ambiguityMap=
            structure(as.character(RNAStringSet(DNAStringSet(IUPAC_CODE_MAP))),
                      names=
                      as.character(RNAStringSet(DNAStringSet(names(IUPAC_CODE_MAP))))),
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

setGeneric("consensusViews", signature="x",
    function(x, ...) standardGeneric("consensusViews")
)

setMethod("consensusViews","MultipleAlignment",
    function(x, ambiguityMap, threshold)
    {
        cmask <- colmask(x)
        if (length(cmask) > 0)
            colmask(x) <- NULL
        consensus <-
          consensusString(x, ambiguityMap=ambiguityMap, threshold=threshold)
        if (length(consensus) == 0)
            consensus <- ""
        Views(BString(consensus), gaps(cmask, start=1, end=ncol(x)))
    }
)

setMethod("consensusViews","DNAMultipleAlignment",
    function(x, ambiguityMap=IUPAC_CODE_MAP, threshold=0.25)
    {
        callNextMethod(x, ambiguityMap=ambiguityMap, threshold=threshold)
    }
)

setMethod("consensusViews","RNAMultipleAlignment",
    function(x,
             ambiguityMap=as.character(RNAStringSet(DNAStringSet(IUPAC_CODE_MAP))),
             threshold=0.25)
    {
        callNextMethod(x, ambiguityMap=ambiguityMap, threshold=threshold)
    }
)

setMethod("consensusViews","AAMultipleAlignment",
    function(x, ambiguityMap="?", threshold=0.5)
    {
        callNextMethod(x, ambiguityMap=ambiguityMap, threshold=threshold)
    }
)

setMethod("alphabetFrequency","MultipleAlignment",
    function(x, as.prob=FALSE, collapse=FALSE)
    {
        if (collapse) {
            strings <- as(x, sprintf("%sStringSet", xsbasetype(x)))
            m <- callGeneric(strings, as.prob=as.prob, collapse=TRUE)
        } else {
            rmask <- rowmask(x)
            if (length(rmask) > 0)
                rowmask(x) <- NULL
            strings <- as(x, sprintf("%sStringSet", xsbasetype(x)))
            m <- callGeneric(strings, as.prob=as.prob)
            if (length(rmask) > 0)
                m[as.integer(rmask),] <- NA
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
