### =========================================================================
### PDict objects
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "PreprocessedTB" VIRTUAL class.
###

setClass("PreprocessedTB",
    representation(
        "VIRTUAL",
        tb="DNAStringSet",  # always constant width
        dups="Dups"
    )
)

setMethod("length", "PreprocessedTB", function(x) length(x@tb))

setMethod("width", "PreprocessedTB", function(x) width(x@tb))

setGeneric("tb", function(x) standardGeneric("tb"))
setMethod("tb", "PreprocessedTB", function(x) x@tb)

setGeneric("tb.width", function(x) standardGeneric("tb.width"))
setMethod("tb.width", "PreprocessedTB", function(x) width(x@tb)[1])

setGeneric("dups", function(x) standardGeneric("dups"))
setMethod("dups", "PreprocessedTB", function(x) x@dups)

setMethod("initialize", "PreprocessedTB",
    function(.Object, tb, dup2unq)
    {
        .Object@tb <- tb
        .Object@dups <- Dups(dup2unq)
        .Object
    }
)

.PreprocessedTB.showFirstLine <- function(x)
{
    cat("Preprocessed Trusted Band of length ", length(x),
        ", width ", tb.width(x),
        ", and type \"", class(x), "\"\n", sep="")
}

setMethod("duplicated", "PreprocessedTB",
    function(x, incomparables=FALSE, ...) duplicated(dups(x))
)

setMethod("dupFrequency", "PreprocessedTB",
    function(x) dupFrequency(dups(x))
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "Twobit" class.
###
### A low-level container for storing a PreprocessedTB object (preprocessed
### Trusted Band) of type "Twobit".
### With this type of preprocessing, the 2-bit-per-letter signatures of all
### the oligonucleotides in the Trusted Band are computed and the mapping
### from these signatures to the 1-based position of the corresponding
### oligonucleotide in the Trusted Band is stored in a way that allows very
### fast lookup.
###

setClass("Twobit",
    contains="PreprocessedTB",
    representation(
        sign2pos="XInteger",  # length(x@sign2pos) is tb.width(x)^4
        base_codes="integer"
    )
)

setMethod("show", "Twobit",
    function(object)
    {
        .PreprocessedTB.showFirstLine(object)
        cat("  (length of sign2pos lookup table is ",
            length(object@sign2pos), ")\n", sep="")
    }
)

setMethod("initialize", "Twobit",
    function(.Object, tb, dup2unq0)
    {
        base_codes <- codes(tb, baseOnly=TRUE)
        C_ans <- .Call("build_Twobit", tb, dup2unq0, base_codes,
                       PACKAGE="Biostrings")
        .Object <- callNextMethod(.Object, tb, C_ans$dup2unq)
        .Object@sign2pos <- XInteger(0)
        .Object@sign2pos@xdata@xp <- C_ans$twobit_sign2pos_xp
        .Object@base_codes <- base_codes
        .Object
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "ACtree" class.
###
### A low-level container for storing a PreprocessedTB object (preprocessed
### Trusted Band) of type "ACtree".
### With this type of preprocessing, all the oligonucleotides in the Trusted
### Band are stored in a 4-ary Aho-Corasick tree (note that this tree is in
### fact an oriented graph if we consider the failure links or the shortcut
### links).
### The number of integers needed to represent a tree node is the number of
### base letters in the DNA alphabet plus 4.
###

setClass("ACtree",
    contains="PreprocessedTB",
    representation(
        nodes="XInteger",
        base_codes="integer"
    )
)

.ACtree.ints_per_acnode <- function(x) (length(x@base_codes) + 4L)

setMethod("show", "ACtree",
    function(object)
    {
        .PreprocessedTB.showFirstLine(object)
        nnodes <- length(object@nodes) %/% .ACtree.ints_per_acnode(object)
        cat("  (number of nodes in the Aho-Corasick tree is ",
            nnodes, ")\n", sep="")
    }
)

### Implement the 0-based subsetting semantic (like in C).
### Typical use:
###   > pdict <- PDict(c("agg", "aca", "gag", "caa", "agt", "aca"))
###   > pdict@actree[0:2] # look at the 3 first nodes
###   > pdict@actree[] # look at all the nodes
###   > flinks0 <- as.matrix(pdict@actree)[ , "flink"]
###   > flinks0 # no failure link is set yet
###   > end_index <- endIndex(matchPDict(pdict, DNAString("acaagagagt")))
###   > flinks1 <- as.matrix(pdict@actree)[ , "flink"]
###   > flinks1 # some failure links have been set
### As you can see the 'pdict' object "learns" from being used!
###   
setMethod("[", "ACtree",
    function(x, i, j, ..., drop)
    {
        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")
        ints_per_acnode <- .ACtree.ints_per_acnode(x)
        nnodes <- length(x@nodes) %/% ints_per_acnode
        if (missing(i)) {
            i <- 0:(nnodes-1)
        } else {
            if (!is.numeric(i))
                stop("invalid subscript type")
            if (any(is.na(i)))
                stop("subscript contains NAs")
            if (!is.integer(i))
                i <- as.integer(i)
        }
        ii <- rep(i * ints_per_acnode, each=ints_per_acnode) +
                          seq_len(ints_per_acnode)
        ans <- matrix(x@nodes[ii], ncol=ints_per_acnode, byrow=TRUE)
        colnames(ans) <- c("parent_id", "depth", x@base_codes,
                           "flink", "P_id")
        rownames(ans) <- i
        ans
    }
)

setMethod("as.matrix", "ACtree", function(x) x[])

setMethod("initialize", "ACtree",
    function(.Object, tb, dup2unq0)
    {
        base_codes <- codes(tb, baseOnly=TRUE)
        on.exit(.Call("free_actree_nodes_buf", PACKAGE="Biostrings"))
        C_ans <- .Call("build_ACtree", tb, dup2unq0, base_codes,
                       PACKAGE="Biostrings")
        .Object <- callNextMethod(.Object, tb, C_ans$dup2unq)
        .Object@nodes <- XInteger(0)
        .Object@nodes@xdata@xp <- C_ans$actree_nodes_xp
        .Object@base_codes <- base_codes
        .Object
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "PDict3Parts" class.
###

setClass("PDict3Parts",
    representation(
        head="DNAStringSet",
        pptb="PreprocessedTB",
        tail="DNAStringSet"
    )
)

setMethod("length", "PDict3Parts", function(x) length(x@pptb))

setMethod("width", "PDict3Parts",
    function(x) { width(x@head) + width(x@pptb) + width(x@tail) }
)

setMethod("head", "PDict3Parts",
    function(x, ...)
    {
        if (all(width(x@head) == 0L))
            return(NULL)
        x@head
    }
)

setMethod("tb", "PDict3Parts", function(x) tb(x@pptb))

setMethod("tb.width", "PDict3Parts", function(x) tb.width(x@pptb))

setMethod("tail", "PDict3Parts",
    function(x, ...)
    {
        if (all(width(x@tail) == 0L))
            return(NULL)
        x@tail
    }
)

setMethod("dups", "PDict3Parts",
    function(x)
        if (is.null(head(x)) && is.null(tail(x))) dups(x@pptb) else NULL
)

.PDict3Parts <- function(x, tb.start, tb.end, tb.width, type, pptb0)
{
    tb <- DNAStringSet(x, start=tb.start, end=tb.end, width=tb.width,
                          use.names=FALSE)
    head_start <- start(x@ranges)
    head_width <- start(tb@ranges) - start(x@ranges)
    head_ranges <- new2("IRanges", start=head_start, width=head_width, check=FALSE)
    head <- new("DNAStringSet", super=super(x), ranges=head_ranges)
    tail_start <- end(tb@ranges) + 1L
    tail_width <- end(x@ranges) - end(tb@ranges)
    tail_ranges <- new2("IRanges", start=tail_start, width=tail_width, check=FALSE)
    tail <- new("DNAStringSet", super=super(x), ranges=tail_ranges)
    use_pptb0 <- !is.null(pptb0) &&
                     all(head_width == 0L) && all(tail_width == 0L) &&
                     type == "ACtree"
    if (use_pptb0) {
        pptb <- pptb0
    } else {
        if (is.null(pptb0))
            dup2unq0 <- NULL
        else
            dup2unq0 <- dups(pptb0)@dup2unq
        pptb <- new(type, tb, dup2unq0)
    }
    new("PDict3Parts", head=head, pptb=pptb, tail=tail)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "PDict" VIRTUAL class (top level).
###
### A (virtual) container for storing a preprocessed dictionary of DNA
### patterns that can later be passed to the matchPDict() function for fast
### matching.
###

setClass("PDict",
    representation(
        "VIRTUAL",
        dict0="DNAStringSet",
        constant_width="logical",
        dups0="Dups"
    )
)

setMethod("length", "PDict", function(x) length(x@dict0))

setMethod("width", "PDict", function(x) width(x@dict0))

setMethod("names", "PDict", function(x) names(x@dict0))

setMethod("dups", "PDict",
    function(x) if (length(x@dups0) == 0) NULL else x@dups0
)

### Extract the i-th element of a PDict object as DNAString object.
### Note that only the "[[" operator is provided for now. Providing "[" sounds
### like a nice feature too but 'x[i]' would have to return a PDict object of
### the same PDict subtype as 'x' i.e. a preprocessed dictionary where the
### preprocessed data structure reflects the subsetted dictionary.
### For example if the preprocessed data structure is an Aho-Corasick tree,
### then this tree needs to be updated so that it stays in sync with
### 'x@dict0[i]'. This updating operation might be complex and expensive in
### terms of CPU cycles and/or memory usage. It could even be that its cost is
### in fact greater than preprocessing again 'x@dict0[i]' from scratch!
### So in the end, "[" would not have much value (other than providing some
### convenience) over the approach that consists to ask the user to do the
### subsetting upstream i.e. to subset the original dictionary before s/he
### passes it to PDict() again.
setMethod("[[", "PDict",
    function(x, i, j, ...)
    {
        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")
        if (missing(i))
            stop("subscript is missing")
        x@dict0[[i]]
    }
)

setReplaceMethod("[[", "PDict",
    function(x, i, j,..., value)
        stop("attempt to modify the value of a ", class(x), " instance")
)

setMethod("duplicated", "PDict",
    function(x, incomparables=FALSE, ...)
    {
        if (is.null(dups(x)))
            stop("duplicates information not available for this object")
        duplicated(dups(x))
    }
)

setMethod("dupFrequency", "PDict",
    function(x)
    {
        if (is.null(dups(x)))
            stop("duplicates information not available for this object")
        dupFrequency(dups(x))
    }
)

### Just an alias for "dupFrequency".
setGeneric("patternFrequency", function(x) standardGeneric("patternFrequency"))
setMethod("patternFrequency", "PDict", function(x) dupFrequency(x))

.PDict.showFirstLine <- function(x, type)
{
    cat(class(x), " object of length ", length(x), sep="")
    if (x@constant_width) {
        width <- width(x@dict0)[1]
        width_info <- paste("width ", width, sep="")
    } else {
        width_info <- "variable width"
    }
    if (length(type) == 1)
        cat(", ", width_info, ", and type \"", type, "\"", sep="")
    else
        cat(" and ", width_info, sep="")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "TB_PDict" class.
###
### A container for storing a Trusted Band PDict object.
###

setClass("TB_PDict",
    contains="PDict",
    representation(
        threeparts="PDict3Parts"
    )
)

setMethod("head", "TB_PDict", function(x, ...) head(x@threeparts))
setMethod("tb", "TB_PDict", function(x) tb(x@threeparts))
setMethod("tb.width", "TB_PDict", function(x) tb.width(x@threeparts))
setMethod("tail", "TB_PDict", function(x, ...) tail(x@threeparts))

setMethod("dups", "TB_PDict",
    function(x)
    {
        ans <- dups(x@threeparts)
        if (!is.null(ans))
            return(ans)
        callNextMethod()
    }
)

setMethod("show", "TB_PDict",
    function(object)
    {
        pdict_type <- class(object@threeparts@pptb)
        .PDict.showFirstLine(object, pdict_type)
        head <- head(object)
        tail <- tail(object)
        if (is.null(head) && is.null(tail))
             return(cat("\n", sep=""))
        cat(":\n")
        if (is.null(head)) {
            cat("  - with NO head")
        } else {
            cat("  - with a head of ")
            min_width <- min(width(head))
            max_width <- max(width(head))
            if (min_width == max_width)
                cat("width ", min_width, sep="")
            else
                cat("variable width (min=", min_width,
                    " / max=", max_width, ")", sep="")
        }
        cat("\n")
        cat("  - with a Trusted Band of width ", tb.width(object), sep="")
        cat("\n")
        if (is.null(tail)) {
            cat("  - with NO tail")
        } else {
            cat("  - with a tail of ")
            min_width <- min(width(tail))
            max_width <- max(width(tail))
            if (min_width == max_width)
                cat("width ", min_width, sep="")
            else
                cat("variable width (min=", min_width,
                    " / max=", max_width, ")", sep="")
        }
        cat("\n")
    }
)

.TB_PDict <- function(x, tb.start, tb.end, tb.width, type)
{
    constant_width <- min(width(x)) == max(width(x))
    if (constant_width && hasOnlyBaseLetters(x))
        pptb0 <- new("ACtree", x, NULL)
    else
        pptb0 <- NULL
    threeparts <- .PDict3Parts(x, tb.start, tb.end, tb.width, type, pptb0)
    ans <- new("TB_PDict", dict0=x,
                           constant_width=constant_width,
                           threeparts=threeparts)
    if (is.null(dups(threeparts)) && !is.null(pptb0))
        ans@dups0 <- dups(pptb0)
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "MTB_PDict" class.
###
### A container for storing a Multiple Trusted Band PDict object.
###

setClass("MTB_PDict",
    contains="PDict",
    representation(
        threeparts_list="list"
    )
)

setMethod("as.list", "MTB_PDict",
    function(x, ...)
    {
        lapply(x@threeparts_list,
          function(threeparts)
              new("TB_PDict", dict0=x@dict0,
                              constant_width=x@constant_width,
                              dups0=x@dups0,
                              threeparts=threeparts)
        )
    }
)

setMethod("show", "MTB_PDict",
    function(object)
    {
        cat("  ")
        .PDict.showFirstLine(object, NULL)
        cat("\nComponents:\n")
        show(as.list(object))
    }
)

### 'max.mismatch' is assumed to be an integer >= 1
.MTB_PDict <- function(x, max.mismatch, type)
{
    min.TBW <- 3L
    min_width <- min(width(x))
    if (min_width < 2L * min.TBW)
        stop("'max.mismatch >= 1' is supported only if the width ",
             "of dictionary 'x' is >= ", 2L * min.TBW)
    constant_width <- min_width == max(width(x))
    NTB <- max.mismatch + 1L # nb of Trusted Bands
    TBW0 <- min_width %/% NTB
    if (TBW0 < min.TBW) {
        max.max.mismatch <- min_width %/% min.TBW - 1L
        stop("'max.mismatch' must be <= ", max.max.mismatch,
             " given the width of dictionary 'x'")
    }
    all_tbw0 <- rep.int(TBW0, NTB - min_width %% NTB)
    all_tbw1 <- rep.int(TBW0 + 1L, min_width %% NTB)
    all_tbw <- c(all_tbw0, all_tbw1)
    ## R1 is the average number of (perfect) matches that is expected
    ## to be found during stage1 each time the sliding window is moved
    ## to the next position.
    R1 <- length(x) * sum(1/4^all_tbw)
    if (R1 > 10)
        warning("given the characteristics of dictionary 'x', ",
                "this value of 'max.mismatch' will\n",
                "  give poor performance when you call ",
                "matchPDict() on this MTB_PDict object\n",
                "  (it will of course depend ultimately on the ",
                "length of the subject)")
    all_headw <- diffinv(all_tbw)
    if (constant_width)
        pptb0 <- new("ACtree", x, NULL)
    else
        pptb0 <- NULL
    threeparts_list <- lapply(seq_len(NTB),
                         function(i)
                           .PDict3Parts(x, all_headw[i]+1L, all_headw[i+1], NA, type, pptb0)
                       )
    ans <- new("MTB_PDict", dict0=x,
                            constant_width=constant_width,
                            threeparts_list=threeparts_list)
    if (!is.null(pptb0))
        ans@dups0 <- dups(pptb0)
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The PDict() constructor (user-friendly).
###

.PDict <- function(x, max.mismatch, tb.start, tb.end, tb.width,
                      type, skip.invalid.patterns)
{
    if (!is(x, "DNAStringSet"))
        x <- DNAStringSet(x)
    if (length(x) == 0)
        stop("'x' must contain at least one pattern")
    names <- names(x)
    if (!is.null(names)) {
        if (any(names %in% c("", NA)))
            stop("'x' has invalid names")
        if (any(duplicated(names)))
            stop("'x' has duplicated names")
    }
    if (!isSingleNumberOrNA(max.mismatch))
        stop("'max.mismatch' must be a single integer or 'NA'")
    if (!is.integer(max.mismatch))
        max.mismatch <- as.integer(max.mismatch)
    if (!isSingleNumberOrNA(tb.start))
        stop("'tb.start' must be a single integer or 'NA'")
    if (!isSingleNumberOrNA(tb.end))
        stop("'tb.end' must be a single integer or 'NA'")
    if (!isSingleNumberOrNA(tb.width))
        stop("'tb.width' must be a single integer or 'NA'")
    if (!is.character(type))
        stop("'type' must be a character vector")
    if (!identical(skip.invalid.patterns, FALSE))
        stop("'skip.invalid.patterns' must be FALSE for now, sorry")
    is_default_TB <- is.na(tb.start) && is.na(tb.end) && is.na(tb.width)
    if (!is.na(max.mismatch) && !is_default_TB)
            stop("'tb.start', 'tb.end' and 'tb.width' must be NAs ",
                 "when 'max.mismatch' is not NA")
    if (is.na(max.mismatch) || max.mismatch == 0) {
        .TB_PDict(x, tb.start, tb.end, tb.width, type)
    } else {
        if (max.mismatch < 0)
            stop("'max.mismatch' must be 'NA' or >= 0")
        .MTB_PDict(x, max.mismatch, type)
    }
}

setGeneric("PDict", signature="x",
    function(x, max.mismatch=NA, tb.start=NA, tb.end=NA, tb.width=NA,
                type="ACtree", skip.invalid.patterns=FALSE)
        standardGeneric("PDict")
)

setMethod("PDict", "character",
    function(x, max.mismatch=NA, tb.start=NA, tb.end=NA, tb.width=NA,
                type="ACtree", skip.invalid.patterns=FALSE)
        .PDict(x, max.mismatch, tb.start, tb.end, tb.width,
                  type, skip.invalid.patterns)
)

setMethod("PDict", "DNAStringSet",
    function(x, max.mismatch=NA, tb.start=NA, tb.end=NA, tb.width=NA,
                type="ACtree", skip.invalid.patterns=FALSE)
        .PDict(x, max.mismatch, tb.start, tb.end, tb.width,
                  type, skip.invalid.patterns)
)

setMethod("PDict", "XStringViews",
    function(x, max.mismatch=NA, tb.start=NA, tb.end=NA, tb.width=NA,
                type="ACtree", skip.invalid.patterns=FALSE)
    {
        if (!is(subject(x), "DNAString"))
            stop("'subject(x)' must be a DNAString object")
        .PDict(x, max.mismatch, tb.start, tb.end, tb.width,
                  type, skip.invalid.patterns)
    }
)

### Just because of those silly "AsIs" objects found in the probe packages
### (e.g. drosophila2probe$sequence)
setMethod("PDict", "AsIs",
    function(x, max.mismatch=NA, tb.start=NA, tb.end=NA, tb.width=NA,
                type="ACtree", skip.invalid.patterns=FALSE)
        .PDict(x, max.mismatch, tb.start, tb.end, tb.width,
                  type, skip.invalid.patterns)
)

