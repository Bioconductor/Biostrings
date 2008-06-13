### =========================================================================
### PDict objects
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "ACtree" class (NOT EXPORTED).
###
### A low-level container for storing the Aho-Corasick tree built from the
### input dictionary (note that this tree is actually an oriented graph if we
### consider the failure links). When the input is based on an alphabet of 4
### letters (DNA input) then the Aho-Corasick tree is a 4-ary tree and the
### base_codes slot is a vector of 4 integers that will be filled with the
### internal codes (byte values) of the unique letters found in the input.
### The number of integers needed to represent a tree node is the number of
### letters in the alphabet plus 3 (i.e. 7 for DNA input) hence this internal
### representation of the Aho-Corasick tree can only be used for input based
### on a very small alphabet.
###
### Slot description:
###
###   nodes: an external integer vector (XInteger object) for the storage
###       of the nodes.
### 
###   base_codes: an integer vector filled with the internal codes (byte
###       values) of the unique letters found in the input dictionary during
###       its preprocessing.

setClass("ACtree",
    representation(
        nodes="XInteger",
        base_codes="integer"
    )
)

.ACtree.ints_per_acnode <- function(x) (length(x@base_codes) + 4L)

setMethod("length", "ACtree", function(x) length(x@nodes) %/% .ACtree.ints_per_acnode(x))

setMethod("show", "ACtree",
    function(object)
    {
        cat(length(object), "-node ", length(object@base_codes), "-ary Aho-Corasick tree\n", sep="")
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
        lx <- length(x)
        if (missing(i)) {
            i <- 0:(lx-1)
        } else {
            if (!is.numeric(i))
                stop("invalid subscript type")
            if (any(is.na(i)))
                stop("subscript contains NAs")
            if (!is.integer(i))
                i <- as.integer(i)
        }
        ints_per_acnode <- .ACtree.ints_per_acnode(x)
        ii <- rep(i * ints_per_acnode, each=ints_per_acnode) + seq_len(ints_per_acnode)
        ans <- matrix(x@nodes[ii], ncol=ints_per_acnode, byrow=TRUE)
        colnames(ans) <- c("parent_id", "depth", x@base_codes,
                           "flink", "P_id")
        rownames(ans) <- i
        ans
    }
)

setMethod("as.matrix", "ACtree", function(x) x[])

setMethod("initialize", "ACtree",
    function(.Object, xp, base_codes)
    {
        .Object@nodes <- XInteger(1)
        .Object@nodes@xp <- xp
        .Object@base_codes <- base_codes
        .Object
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "Dups" class.
###
### A general container for storing the mapping between duplicated elements
### of an arbitrary vector.

setClass("Dups",
    representation(
        dup2unq="integer",  # many-to-one integer mapping
        unq2dup="list"      # one-to-many integer mapping
    )
)

.reverse.dup2unq <- function(dup2unq)
{
    ans <- vector(mode="list", length=length(dup2unq))
    sparse_ans <- split(seq_along(dup2unq), dup2unq)
    ans[as.integer(names(sparse_ans))] <- sparse_ans
    ans
}

.valid.Dups <- function(object)
{
    if (!is.integer(object@dup2unq))
        return("the 'dup2unq' slot must contain an integer vector")
    if (!all(object@dup2unq >= 1L, na.rm=TRUE))
        return("the 'dup2unq' slot must contain integer values >= 1")
    if (!all(object@dup2unq < seq_along(object@dup2unq), na.rm=TRUE))
        return("when mapped, values in the 'dup2unq' slot must be mapped to lower values")
    if (!all(is.na(object@dup2unq[object@dup2unq])))
        return("when mapped, values in the 'dup2unq' slot must be mapped to unmapped values")
    if (!is.list(object@unq2dup))
        return("the 'unq2dup' slot must contain a list")
    if (length(object@dup2unq) != length(object@unq2dup))
        return("the 'dup2unq' and 'unq2dup' slots must have the same length")
    if (!identical(.reverse.dup2unq(object@dup2unq), object@unq2dup))
        return("the 'unq2dup' slot must contain the reverse mapping of the 'dup2unq' slot")
    NULL
}
setValidity("Dups",
    function(object)
    {
        problems <- .valid.Dups(object)
        if (is.null(problems)) TRUE else problems
    }
)

Dups <- function(dup2unq)
    new("Dups", dup2unq=dup2unq, unq2dup=.reverse.dup2unq(dup2unq))

setMethod("length", "Dups", function(x) length(x@dup2unq))

setMethod("duplicated", "Dups",
    function(x, incomparables=FALSE, ...) !is.na(x@dup2unq)
)

setGeneric("dupFrequency", function(x) standardGeneric("dupFrequency"))

setMethod("dupFrequency", "Dups",
    function(x)
    {
        ans <- rep.int(1L, length(x@unq2dup))
        mapped_unqs <- setdiff(unique(x@dup2unq), NA)
        for (unq in mapped_unqs) {
            ii <- as.integer(c(unq, x@unq2dup[[unq]]))
            ans[ii] <- length(ii)
        }
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "PDict" class.
###
### Very general container for any kind of preprocessed pattern dictionary.
###

setClass("PDict", representation("VIRTUAL"))

setGeneric("patternFrequency", function(x) standardGeneric("patternFrequency"))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "ACtree_PDict" class.
###
### A container for storing a preprocessed constant width dictionary (or set)
### of DNA patterns.
###
### Slot description:
###
###   length: the dictionary length L i.e. the number of patterns
###       (e.g. L=10^6).
###
###   width: the dictionary width W i.e. the number of chars per pattern
###       (e.g. W=25).
###
###   actree: the Aho-Corasick tree built from the input dictionary.
###

setClass("ACtree_PDict",
    contains="PDict",
    representation(
        length="integer",
        width="integer",
        actree="ACtree",
        dups="Dups",
        NAMES="character",  # R doesn't like @names !!
        geom="list"
    )
)

### Not intended to be used directly by the user.
setMethod("initialize", "ACtree_PDict",
    function(.Object, length, pp_Cans, names)
    {
        .Object@length <- length
        .Object@width <- pp_Cans$geom$width
        .Object@actree <- new("ACtree", pp_Cans$actree_nodes_xp, pp_Cans$actree_base_codes)
        .Object@dups <- Dups(pp_Cans$dup2unq)
        .Object@NAMES <- if (is.null(names)) as.character(NA) else names
        .Object@geom <- pp_Cans$geom
        .Object
    }
)

setMethod("length", "ACtree_PDict", function(x) x@length)

setMethod("width", "ACtree_PDict", function(x) x@width)

setMethod("names", "ACtree_PDict",
    function(x)
        if (length(x@NAMES) == 1 && is.na(x@NAMES)) NULL else x@NAMES
)

setMethod("show", "ACtree_PDict",
    function(object)
    {
        cat(length(object), "-pattern constant width PDict object of width ",
            width(object), sep="")
        if (is.null(names(object)))
            cat(" (patterns have no names)")
        cat("\n")
    }
)

setMethod("duplicated", "ACtree_PDict",
    function(x, incomparables=FALSE, ...) duplicated(x@dups)
)

setMethod("dupFrequency", "ACtree_PDict",
    function(x) dupFrequency(x@dups)
)
setMethod("patternFrequency", "ACtree_PDict",
    function(x) dupFrequency(x@dups)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "TBdna_PDict" class.
###
### A container for storing a preprocessed "Trusted Band" dictionary (or
### set) of DNA patterns.
### 2 important particular cases of "Trusted Band" dictionaries are:
###   1) "Trusted Prefix": a "Trusted Prefix" dictionary is a "Trusted Band"
###      dictionary with no head.
###   2) "Trusted Suffix": a "Trusted Suffix" dictionary is a "Trusted Band"
###      dictionary with no tail.
###

### 'head' and 'tail' must be DNAStringSet objects of the same length.
### One of them can be NULL but they can't both be NULL at the same time.
setClass("TBdna_PDict",
    contains="ACtree_PDict",
    representation(
        head="DNAStringSet", # can be NULL
        tail="DNAStringSet"  # can be NULL
    )
)

### Not intended to be used directly by the user.
setMethod("initialize", "TBdna_PDict",
    function(.Object, length, pp_Cans, names)
        callNextMethod(.Object, length, pp_Cans, names)
)

setMethod("head", "TBdna_PDict", function(x, ...) x@head)
setMethod("tail", "TBdna_PDict", function(x, ...) x@tail)

setMethod("show", "TBdna_PDict",
    function(object)
    {
        if (is.null(object@head))
            trusted_part <- "Prefix"
        else if (is.null(object@tail))
            trusted_part <- "Suffix"
        else
            trusted_part <- "Band"
        cat(length(object), "-pattern PDict object with a Trusted ",
            trusted_part, " of width ", width(object), sep="")
        if (is.null(names(object)))
            cat(" (patterns have no names)")
        cat("\n")
    }
)

setMethod("duplicated", "TBdna_PDict",
    function(x, incomparables=FALSE, ...)
    {
        stop("not ready yet, sorry!")
    }
)

setMethod("patternFrequency", "TBdna_PDict",
    function(x)
    {
        stop("not ready yet, sorry!")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The PDict() constructor (user-friendly).
###

.PDict <- function(dict, names, tb.start, tb.end,
                   drop.head, drop.tail, skip.invalid.patterns)
{
    if (!is.null(names)) {
        if (any(names %in% c("", NA)))
            stop("'dict' has invalid names")
        if (any(duplicated(names)))
            stop("'dict' has duplicated names")
    }
    if (!isSingleNumberOrNA(tb.start))
        stop("'tb.start' must be a single integer or NA")
    if (!is.integer(tb.start))
        tb.start <- as.integer(tb.start)
    if (!isSingleNumberOrNA(tb.end))
        stop("'tb.end' must be a single integer or NA")
    if (!is.integer(tb.end))
        tb.end <- as.integer(tb.end)
    if (!isTRUEorFALSE(drop.head))
        stop("'drop.head' must be 'TRUE' or 'FALSE'")
    if (!isTRUEorFALSE(drop.tail))
        stop("'drop.tail' must be 'TRUE' or 'FALSE'")
    on.exit(.Call("CWdna_free_actree_nodes_buf", PACKAGE="Biostrings"))
    if (is.character(dict)) {
        pp_Cans <- .Call("build_ACtree_PDict_from_CHARACTER",
                         dict,
                         tb.start, tb.end,
                         PACKAGE="Biostrings")
    } else if (is(dict, "DNAStringSet")) {
        pp_Cans <- .Call("build_ACtree_PDict_from_XStringSet",
                         dict,
                         tb.start, tb.end,
                         PACKAGE="Biostrings")
    } else {
        stop("unsuported 'dict' type")
    }
    if (is.na(tb.start) || tb.start == 1L)
        drop.head <- TRUE
    if (is.na(tb.end) || tb.end == -1L)
        drop.tail <- TRUE
    if (drop.head && drop.tail)
        return(new("ACtree_PDict", length(dict), pp_Cans, names))
    ans <- new("TBdna_PDict", length(dict), pp_Cans, names)
    if (!drop.head) 
        ans@head <- DNAStringSet(dict, end=tb.start-1L)
    if (!drop.tail)
        ans@tail <- DNAStringSet(dict, start=tb.end+1L)
    ans
}

### The input dictionary 'dict' must be:
###   - of length >= 1
### The supported types for the input dictionary are:
###   - character vector
###   - DNAStringSet object
###   - XStringViews object
### Typical use:
###   > library(Biostrings)
###   > dict <- c("abb", "aca", "bab", "caa", "abd", "aca")
###   > pdict <- PDict(dict)
###
setGeneric("PDict", signature="dict",
    function(dict, tb.start=1, tb.end=NA,
             drop.head=FALSE, drop.tail=FALSE, skip.invalid.patterns=FALSE)
        standardGeneric("PDict")
)

setMethod("PDict", "character",
    function(dict, tb.start=1, tb.end=NA,
             drop.head=FALSE, drop.tail=FALSE, skip.invalid.patterns=FALSE)
    {
        if (length(dict) == 0)
            stop("'dict' must contain at least one pattern")
        .PDict(dict, names(dict), tb.start, tb.end,
               drop.head, drop.tail, skip.invalid.patterns)
    }
)

setMethod("PDict", "DNAStringSet",
    function(dict, tb.start=1, tb.end=NA,
             drop.head=FALSE, drop.tail=FALSE, skip.invalid.patterns=FALSE)
    {
        if (length(dict) == 0)
            stop("'dict' must contain at least one pattern")
        .PDict(dict, names(dict), tb.start, tb.end,
               drop.head, drop.tail, skip.invalid.patterns)
    }
)

setMethod("PDict", "XStringViews",
    function(dict, tb.start=1, tb.end=NA,
             drop.head=FALSE, drop.tail=FALSE, skip.invalid.patterns=FALSE)
    {
        dict <- DNAStringSet(dict)
        PDict(dict, tb.start=tb.start, tb.end=tb.end,
                    drop.head=drop.head, drop.tail=drop.tail,
                    skip.invalid.patterns=skip.invalid.patterns)
    }
)

### Just because of those silly "AsIs" objects found in the probe packages
### (e.g. drosophila2probe$sequence)
setMethod("PDict", "AsIs",
    function(dict, tb.start=1, tb.end=NA,
             drop.head=FALSE, drop.tail=FALSE, skip.invalid.patterns=FALSE)
    {
        if (!is.character(dict))
            stop("unsuported 'dict' type")
        class(dict) <- "character" # keeps the names (unlike as.character())
        PDict(dict, tb.start=tb.start, tb.end=tb.end,
                    drop.head=drop.head, drop.tail=drop.tail,
                    skip.invalid.patterns=skip.invalid.patterns)
    }
)

