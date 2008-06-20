### =========================================================================
### PDict objects
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "Dups" class.
###
### A general container for storing the direct and reverse mappings between
### duplicated and unique elements of an arbitrary vector.
###

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
        return("when mapped, values in the 'dup2unq' slot must be mapped ",
               "to lower values")
    if (!all(is.na(object@dup2unq[object@dup2unq])))
        return("when mapped, values in the 'dup2unq' slot must be mapped ",
               "to unmapped values")
    if (!is.list(object@unq2dup))
        return("the 'unq2dup' slot must contain a list")
    if (length(object@dup2unq) != length(object@unq2dup))
        return("the 'dup2unq' and 'unq2dup' slots must have the same length")
    if (!identical(.reverse.dup2unq(object@dup2unq), object@unq2dup))
        return("the 'unq2dup' slot must contain the reverse mapping ",
               "of the 'dup2unq' slot")
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

setMethod("show", "Dups",
    function(object)
    {
        percentage <- 100 * sum(duplicated(object)) / length(object)
        cat(class(object), " object of length ", length(object),
            " (", percentage, "% of duplicates)\n", sep="")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "PDict" class.
###
### A (virtual) container for storing a preprocessed dictionary of DNA
### patterns that can later be passed to the matchPDict() function for fast
### matching.
### Derive this virtual class for each type of preprocessing that you want
### to implement.
###

setClass("PDict",
    representation(
        "VIRTUAL",
        dict0="DNAStringSet",
        head="DNAStringSet",
        tb="DNAStringSet",    # the Trusted Band (always of constant width)
        tail="DNAStringSet",
        tb.dups="Dups"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### PDict accessor methods.
###

setMethod("length", "PDict", function(x) length(x@dict0))

setMethod("width", "PDict", function(x) width(x@dict0))

setMethod("names", "PDict", function(x) names(x@dict0))

setMethod("head", "PDict",
    function(x, ...)
    {
        if (all(width(x@head) == 0L))
            return(NULL)
        x@head
    }
)

setGeneric("tb", function(x) standardGeneric("tb"))
setMethod("tb", "PDict", function(x) x@tb)

setGeneric("tb.width", function(x) standardGeneric("tb.width"))
setMethod("tb.width", "PDict", function(x) width(x@tb)[1])

setMethod("tail", "PDict",
    function(x, ...)
    {
        if (all(width(x@tail) == 0L))
            return(NULL)
        x@tail
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting a PDict object.
###
### Only the "[[" operator is provided for now.
### Providing "[" sounds nice too but it has to return a PDict object of
### the same type made of the selected patterns. For an ACtree_PDict object
### this means that the @actree slot must be updated and this could cost
### (in terms of amount of CPU and memory used) as much as preprocessing
### again the entire set of selected patterns.
### So in the end, "[" would not have much value (except some convenience)
### over the solution that consists to ask the user to subset upstream i.e.
### to subset the original dictionary before s/he passes it to PDict().
### 

### Extract the i-th element of a PDict object as DNAString object.
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
    {
        stop("attempt to modify the value of a ", class(x), " instance")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### PDict initialization.
###

### Not intended to be used directly by the user.
setMethod("initialize", "PDict",
    function(.Object, dict0, head, tb, tail, dup2unq)
    {
        .Object@dict0 <- dict0
        .Object@head <- head
        .Object@tb <- tb
        .Object@tail <- tail
        .Object@tb.dups <- Dups(dup2unq)
        .Object
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The PDict "show" method.
###

setMethod("show", "PDict",
    function(object)
    {
        cat(length(object), "-pattern PDict object", sep="")
        head <- head(object)
        tail <- tail(object)
        if (is.null(head) && is.null(tail)) {
             cat(" of width ", tb.width(object), "\n", sep="")
             return(invisible(NULL))
        }
        cat(":\n")
        if (is.null(head)) {
            cat("  - with NO head")
        } else {
            cat("  - with a head of ")
            min_width <- min(width(head))
            max_width <- max(width(head))
            if (min_width == max_width)
                cat("constant width ", min_width, sep="")
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
                cat("constant width ", min_width, sep="")
            else
                cat("variable width (min=", min_width,
                    " / max=", max_width, ")", sep="")
        }
        cat("\n")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### PDict other methods.
###

setMethod("duplicated", "PDict",
    function(x, incomparables=FALSE, ...)
    {
        head <- head(x)
        tail <- tail(x)
        if (is.null(head) && is.null(tail))
            return(duplicated(x@tb.dups))
        stop("duplicated() is not yet working on a PDict object ",
             "with a head or a tail, sorry!")
    }
)

setMethod("dupFrequency", "PDict",
    function(x)
    {
        head <- head(x)
        tail <- tail(x)
        if (is.null(head) && is.null(tail))
            return(dupFrequency(x@tb.dups))
        stop("dupFrequency() is not yet working on a PDict object ",
             "with a head or a tail, sorry!")
    }
)

setGeneric("patternFrequency", function(x) standardGeneric("patternFrequency"))

setMethod("patternFrequency", "PDict", function(x) dupFrequency(x))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "Twobit" class.
###
### A low-level container for storing the mapping from the 2bit-per-letter
### signatures of a set of oligonucleotides of constant width to the 1-based
### positions of these oligonucleotides in the set.
###
### Slots:
###
###   width:
###
###   sign2pos: XInteger object of length width^4
### 
###   base_codes:
###

setClass("Twobit",
    representation(
        width="integer",
        sign2pos="XInteger",
        base_codes="integer"
    )
)

setMethod("length", "Twobit", function(x) length(x@sign2pos))

setMethod("width", "Twobit", function(x) x@width)

setMethod("show", "Twobit",
    function(object)
    {
        cat(class(object), " object of width ", width(object),
            " and length ", length(object), "\n", sep="")
    }
)

setMethod("initialize", "Twobit",
    function(.Object, width, sign2pos_xp, base_codes)
    {
        .Object@width <- width
        .Object@sign2pos <- XInteger(0)
        .Object@sign2pos@xp <- sign2pos_xp
        .Object@base_codes <- base_codes
        .Object
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "Twobit_PDict" class.
###
### A specific PDict container where the Trusted Band is stored in a
### Twobit object.
###

setClass("Twobit_PDict",
    contains="PDict",
    representation(
        twobit="Twobit"
    )
)

### Not intended to be used directly by the user.
setMethod("initialize", "Twobit_PDict",
    function(.Object, dict0, head, tb, tail, base_codes)
    {
        C_ans <- .Call("build_Twobit_PDict", tb, base_codes,
                       PACKAGE="Biostrings")
        .Object <- callNextMethod(.Object, dict0, head, tb, tail,
                                           C_ans$dup2unq)
        .Object@twobit <- new("Twobit", width(tb)[1],
                                        C_ans$twobit_sign2pos_xp,
                                        base_codes)
        .Object
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "ACtree" class.
###
### A low-level container for storing the Aho-Corasick tree built from the
### input dictionary (note that this tree is actually an oriented graph if we
### consider the failure links). When the input is based on an alphabet of 4
### letters (DNA input) then the Aho-Corasick tree is a 4-ary tree and the
### base_codes slot is a vector of 4 integers that will be filled with the
### internal codes (byte values) of the unique letters found in the input.
### The number of integers needed to represent a tree node is the number of
### letters in the alphabet plus 4 (i.e. 4 for DNA input) hence this internal
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
###

setClass("ACtree",
    representation(
        nodes="XInteger",
        base_codes="integer"
    )
)

.ACtree.ints_per_acnode <- function(x) (length(x@base_codes) + 4L)

setMethod("length", "ACtree",
    function(x) length(x@nodes) %/% .ACtree.ints_per_acnode(x)
)

setMethod("show", "ACtree",
    function(object)
    {
        cat(length(object), "-node ",
            length(object@base_codes), "-ary Aho-Corasick tree\n",
            sep="")
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
    function(.Object, nodes_xp, base_codes)
    {
        .Object@nodes <- XInteger(0)
        .Object@nodes@xp <- nodes_xp
        .Object@base_codes <- base_codes
        .Object
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "ACtree_PDict" class.
###
### A specific PDict container where the Trusted Band is stored in an
### ACtree object.
###

setClass("ACtree_PDict",
    contains="PDict",
    representation(
        actree="ACtree"
    )
)

### Not intended to be used directly by the user.
setMethod("initialize", "ACtree_PDict",
    function(.Object, dict0, head, tb, tail, base_codes)
    {
        on.exit(.Call("free_actree_nodes_buf", PACKAGE="Biostrings"))
        C_ans <- .Call("build_ACtree_PDict", tb, base_codes,
                       PACKAGE="Biostrings")
        .Object <- callNextMethod(.Object, dict0, head, tb, tail,
                                           C_ans$dup2unq)
        .Object@actree <- new("ACtree", C_ans$actree_nodes_xp, base_codes)
        .Object
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The PDict() constructor (user-friendly).
###

.PDict <- function(x, tb.start, tb.end, tb.width, type, skip.invalid.patterns)
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
    tb <- DNAStringSet(x, start=tb.start, end=tb.end, width=tb.width,
                          use.names=FALSE)
    head_start <- start(x)
    head_width <- start(tb) - start(x)
    head <- new("DNAStringSet", super(x), head_start, head_width, names=NULL)
    tail_start <- end(tb) + 1L
    tail_width <- end(x) - end(tb)
    tail <- new("DNAStringSet", super(x), tail_start, tail_width, names=NULL)
    base_codes <- codes(tb, baseOnly=TRUE)
    if (type == "ACtree") {
        ans <- new("ACtree_PDict", x, head, tb, tail, base_codes)
    } else if (type == "Twobit") {
        ans <- new("Twobit_PDict", x, head, tb, tail, base_codes)
    } else {
        stop("unknown type")
    }
    ans
}

setGeneric("PDict", signature="x",
    function(x, tb.start=NA, tb.end=NA, tb.width=NA,
                type="ACtree", skip.invalid.patterns=FALSE)
        standardGeneric("PDict")
)

setMethod("PDict", "character",
    function(x, tb.start=NA, tb.end=NA, tb.width=NA,
                type="ACtree", skip.invalid.patterns=FALSE)
        .PDict(x, tb.start, tb.end, tb.width, type, skip.invalid.patterns)
)

setMethod("PDict", "DNAStringSet",
    function(x, tb.start=NA, tb.end=NA, tb.width=NA,
                type="ACtree", skip.invalid.patterns=FALSE)
        .PDict(x, tb.start, tb.end, tb.width, type, skip.invalid.patterns)
)

setMethod("PDict", "XStringViews",
    function(x, tb.start=NA, tb.end=NA, tb.width=NA,
                type="ACtree", skip.invalid.patterns=FALSE)
    {
        if (!is(subject(x), "DNAString"))
            stop("'subject(x)' must be a DNAString object")
        .PDict(x, tb.start, tb.end, tb.width, type, skip.invalid.patterns)
    }
)

### Just because of those silly "AsIs" objects found in the probe packages
### (e.g. drosophila2probe$sequence)
setMethod("PDict", "AsIs",
    function(x, tb.start=NA, tb.end=NA, tb.width=NA,
                type="ACtree", skip.invalid.patterns=FALSE)
        .PDict(x, tb.start, tb.end, tb.width, type, skip.invalid.patterns)
)

