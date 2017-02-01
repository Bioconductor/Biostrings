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
        exclude_dups0="logical",
        dups="Dups",
        base_codes="integer"
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

setGeneric("nnodes", function(x) standardGeneric("nnodes"))

setGeneric("hasAllFlinks", function(x) standardGeneric("hasAllFlinks"))

setGeneric("computeAllFlinks",
    function(x, ...) standardGeneric("computeAllFlinks"))

setMethod("initialize", "PreprocessedTB",
    function(.Object, tb, pp_exclude, high2low, base_codes)
    {
        .Object@tb <- tb
        .Object@exclude_dups0 <- !is.null(pp_exclude)
        .Object@dups <- Dups(high2low)
        .Object@base_codes <- base_codes  # should be 'xscodes(tb, baseOnly=TRUE)'
        .Object
    }
)

.PreprocessedTB.showFirstLine <- function(x)
{
    cat("Preprocessed Trusted Band\n")
    cat("| length x width = ", length(x), " x ", tb.width(x), "\n", sep="")
    cat("| algorithm = \"", class(x), "\"\n", sep="")
}

setMethod("togrouplength", "PreprocessedTB",
    function(x, j=NULL) togrouplength(dups(x), j=j)
)

setMethod("duplicated", "PreprocessedTB",
    function(x, incomparables=FALSE, ...) duplicated(dups(x))
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "Twobit" class.
###
### A low-level container for storing the PreprocessedTB object (preprocessed
### Trusted Band) obtained with the "Twobit" algo.
### With this algo, the 2-bit-per-letter signatures of all
### the oligonucleotides in the Trusted Band are computed and the mapping
### from these signatures to the 1-based position of the corresponding
### oligonucleotide in the Trusted Band is stored in a way that allows very
### fast lookup.
###

setClass("Twobit",
    contains="PreprocessedTB",
    representation(
        sign2pos="XInteger"  # length(x@sign2pos) is tb.width(x)^4
    )
)

setMethod("show", "Twobit",
    function(object)
    {
        .PreprocessedTB.showFirstLine(object)
        cat("| length of sign2pos lookup table = ",
            length(object@sign2pos), "\n", sep="")
    }
)

setMethod("initialize", "Twobit",
    function(.Object, tb, pp_exclude)
    {
        base_codes <- xscodes(tb, baseOnly=TRUE)
        C_ans <- .Call2("build_Twobit", tb, pp_exclude, base_codes,
                       PACKAGE="Biostrings")
        .Object <- callNextMethod(.Object, tb, pp_exclude, C_ans$high2low, base_codes)
        .Object@sign2pos <- C_ans$sign2pos
        .Object
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "ACtree2" class.
###
### A low-level container for storing the PreprocessedTB object (preprocessed
### Trusted Band) obtained with the "ACtree2" algo.
###

### Big Atomic Buffer of integers.
setClass("IntegerBAB", representation(xp="externalptr"))

setClass("ACtree2",
    contains="PreprocessedTB",
    representation(
        nodebuf_ptr="IntegerBAB",
        nodeextbuf_ptr="IntegerBAB"
    )
)

setMethod("nnodes", "ACtree2",
    function(x) .Call2("ACtree2_nnodes", x, PACKAGE="Biostrings")
)

setMethod("hasAllFlinks", "ACtree2",
    function(x) .Call2("ACtree2_has_all_flinks", x, PACKAGE="Biostrings")
)

setMethod("computeAllFlinks", "ACtree2",
    function(x) .Call2("ACtree2_compute_all_flinks", x, PACKAGE="Biostrings")
)

setMethod("show", "ACtree2",
    function(object)
    {
        .PreprocessedTB.showFirstLine(object)
	invisible(.Call2("ACtree2_summary", object, PACKAGE="Biostrings"))
    }
)

setMethod("initialize", "ACtree2",
    function(.Object, tb, pp_exclude)
    {
        nodebuf_max_nblock <- .Call2("ACtree2_nodebuf_max_nblock",
                                    PACKAGE="Biostrings")
        nodebuf_ptr <- .Call2("IntegerBAB_new", nodebuf_max_nblock,
                             PACKAGE="Biostrings")
        nodeextbuf_max_nblock <- .Call2("ACtree2_nodeextbuf_max_nblock",
                                       PACKAGE="Biostrings")
        nodeextbuf_ptr <- .Call2("IntegerBAB_new", nodeextbuf_max_nblock,
                                PACKAGE="Biostrings")
        base_codes <- xscodes(tb, baseOnly=TRUE)
        C_ans <- .Call2("ACtree2_build",
                       tb, pp_exclude, base_codes,
                       nodebuf_ptr, nodeextbuf_ptr,
                       PACKAGE="Biostrings")
        .Object <- callNextMethod(.Object, tb, pp_exclude, C_ans$high2low, base_codes)
        .Object@nodebuf_ptr <- nodebuf_ptr
        .Object@nodeextbuf_ptr <- nodeextbuf_ptr
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

.PDict3Parts <- function(x, tb.start, tb.end, tb.width, algo, pptb0)
{
    threeparts <- threebands(x, start=tb.start, end=tb.end, width=tb.width)
    head <- threeparts$left
    tb <- threeparts$middle
    tail <- threeparts$right
    if (is.null(pptb0)) {
        pptb <- new(algo, tb, NULL)
    } else {
        use_pptb0 <- algo == class(pptb0) &&
                     all(width(head) == 0L) && all(width(tail) == 0L)
        if (use_pptb0) {
            ## We can avoid doing the expensive preprocessing again by
            ## making the 'pptb' that would be returned by
            ## 'new(algo, tb, high2low(dups(pptb0)))':
            pptb <- pptb0
            pptb@dups <- Dups(rep.int(as.integer(NA), length(pptb)))
            pptb@exclude_dups0 <- TRUE
        } else {
            pptb <- new(algo, tb, high2low(dups(pptb0)))
        }
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

### TODO: Use dups0="Dups_OR_NULL" below like for the ByPos_MIndex class (see
### MIndex-class.R).
setClass("PDict",
    contains="List",
    representation(
        "VIRTUAL",
        dict0="DNAStringSet",
        constant_width="logical",
        dups0="Dups"  # TODO: use Dups_OR_NULL instead!
    ),
    prototype(
        elementType="DNAString"
    )
)

setMethod("length", "PDict", function(x) length(x@dict0))

setMethod("width", "PDict", function(x) width(x@dict0))

setMethod("names", "PDict", function(x) names(x@dict0))

setReplaceMethod("names", "PDict",
    function(x, value)
        stop("attempt to modify the names of a ", class(x), " instance")
)

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
        i <- normalizeDoubleBracketSubscript(i, x)
        x@dict0[[i]]
    }
)

setMethod("togrouplength", "PDict",
    function(x, j=NULL)
    {
        if (is.null(dups(x)))
            stop("duplicates information not available for this object")
        togrouplength(dups(x), j=j)
    }
)

setMethod("duplicated", "PDict",
    function(x, incomparables=FALSE, ...)
    {
        if (is.null(dups(x)))
            stop("duplicates information not available for this object")
        duplicated(dups(x))
    }
)

### Just an alias for "togrouplength".
setGeneric("patternFrequency", function(x) standardGeneric("patternFrequency"))
setMethod("patternFrequency", "PDict", function(x) togrouplength(x))

.PDict.showFirstLine <- function(x, algo)
{
    cat(class(x), " object of length ", length(x), sep="")
    if (x@constant_width) {
        width <- width(x@dict0)[1]
        width_info <- paste("width ", width, sep="")
    } else {
        width_info <- "variable width"
    }
    cat(" and ", width_info, sep="")
    if (!is.null(algo))
        cat(" (preprocessing algo=\"", algo, "\")", sep="")
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

setMethod("show", "TB_PDict",
    function(object)
    {
        algo <- class(object@threeparts@pptb)
        .PDict.showFirstLine(object, algo)
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

.TB_PDict <- function(x, tb.start, tb.end, tb.width, algo)
{
    constant_width <- isConstant(width(x))
    if (constant_width && hasOnlyBaseLetters(x))
        pptb0 <- new("ACtree2", x, NULL)  # because ACtree2 supports big input
    else
        pptb0 <- NULL
    threeparts <- .PDict3Parts(x, tb.start, tb.end, tb.width, algo, pptb0)
    ans <- new("TB_PDict", dict0=x,
                           constant_width=constant_width,
                           threeparts=threeparts)
    ## '!is.null(pptb0)' should be the same as 'threeparts@pptb@exclude_dups0'
    if (!is.null(pptb0))
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
.MTB_PDict <- function(x, max.mismatch, algo)
{
    min.TBW <- 3L
    min_width <- min(width(x))
    if (min_width < 2L * min.TBW)
        stop("'max.mismatch >= 1' is supported only if the width ",
             "of dictionary 'x' is >= ", 2L * min.TBW)
    constant_width <- isConstant(width(x))
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
        pptb0 <- new("ACtree2", x, NULL)  # because ACtree2 supports big input
    else
        pptb0 <- NULL
    threeparts_list <- lapply(seq_len(NTB),
                         function(i)
                           .PDict3Parts(x, all_headw[i]+1L, all_headw[i+1L], NA, algo, pptb0)
                       )
    ans <- new("MTB_PDict", dict0=x,
                            constant_width=constant_width,
                            threeparts_list=threeparts_list)
    if (!is.null(pptb0))
        ans@dups0 <- dups(pptb0)
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "Expanded_TB_PDict" class.
###

setClass("Expanded_TB_PDict",
    contains="TB_PDict",
    representation(
        expanded_dict0="DNAStringSetList"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The PDict() constructor (user-friendly).
###

.PDict <- function(x, max.mismatch, tb.start, tb.end, tb.width,
                      algo, skip.invalid.patterns)
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
    if (!is.character(algo))
        stop("'algorithm' must be a character vector")
    if ("ACtree" %in% algo) {
        warning("support for ACtree preprocessing algo has been ",
                "dropped, using ACtree2 algo")
        algo[!is.na(match(algo, "ACtree"))] <- "ACtree2"
    }
    if (!identical(skip.invalid.patterns, FALSE))
        stop("'skip.invalid.patterns' must be FALSE for now, sorry")
    is_default_TB <- is.na(tb.start) && is.na(tb.end) && is.na(tb.width)
    if (!is.na(max.mismatch) && !is_default_TB)
            stop("'tb.start', 'tb.end' and 'tb.width' must be NAs ",
                 "when 'max.mismatch' is not NA")
    if (is.na(max.mismatch) || max.mismatch == 0) {
        .TB_PDict(x, tb.start, tb.end, tb.width, algo)
    } else {
        if (max.mismatch < 0)
            stop("'max.mismatch' must be 'NA' or >= 0")
        .MTB_PDict(x, max.mismatch, algo)
    }
}

setGeneric("PDict", signature="x",
    function(x, max.mismatch=NA, tb.start=NA, tb.end=NA, tb.width=NA,
                algorithm="ACtree2", skip.invalid.patterns=FALSE)
        standardGeneric("PDict")
)

setMethod("PDict", "character",
    function(x, max.mismatch=NA, tb.start=NA, tb.end=NA, tb.width=NA,
                algorithm="ACtree2", skip.invalid.patterns=FALSE)
        .PDict(x, max.mismatch, tb.start, tb.end, tb.width,
                  algorithm, skip.invalid.patterns)
)

setMethod("PDict", "DNAStringSet",
    function(x, max.mismatch=NA, tb.start=NA, tb.end=NA, tb.width=NA,
                algorithm="ACtree2", skip.invalid.patterns=FALSE)
        .PDict(x, max.mismatch, tb.start, tb.end, tb.width,
                  algorithm, skip.invalid.patterns)
)

setMethod("PDict", "XStringViews",
    function(x, max.mismatch=NA, tb.start=NA, tb.end=NA, tb.width=NA,
                algorithm="ACtree2", skip.invalid.patterns=FALSE)
    {
        if (!is(subject(x), "DNAString"))
            stop("'subject(x)' must be a DNAString object")
        .PDict(x, max.mismatch, tb.start, tb.end, tb.width,
                  algorithm, skip.invalid.patterns)
    }
)

### 2 extra "PDict" methods to deal with the probe sequences stored
### in the *probe annotation packages (e.g. drosophila2probe).
setMethod("PDict", "AsIs",
    function(x, max.mismatch=NA, tb.start=NA, tb.end=NA, tb.width=NA,
                algorithm="ACtree2", skip.invalid.patterns=FALSE)
        .PDict(x, max.mismatch, tb.start, tb.end, tb.width,
                  algorithm, skip.invalid.patterns)
)
setMethod("PDict", "probetable",
    function(x, max.mismatch=NA, tb.start=NA, tb.end=NA, tb.width=NA,
                algorithm="ACtree2", skip.invalid.patterns=FALSE)
        PDict(x$sequence, max.mismatch=max.mismatch,
              tb.start=tb.start, tb.end=tb.end, tb.width=tb.width,
              algorithm=algorithm, skip.invalid.patterns=skip.invalid.patterns)
)

