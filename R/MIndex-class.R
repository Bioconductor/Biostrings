### =========================================================================
### MIndex objects
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "MIndex" API.
### 
### An API for manipulating the match container returned by matchPDict().
###

setClass("MIndex", representation("VIRTUAL"))


setGeneric("startIndex", signature="x",
    function(x, all.names=FALSE) standardGeneric("startIndex"))

setGeneric("endIndex", signature="x",
    function(x, all.names=FALSE) standardGeneric("endIndex"))

setGeneric("countIndex", signature="x",
    function(x, all.names=FALSE) standardGeneric("countIndex"))

setMethod("countIndex", "MIndex",
    function(x, all.names=FALSE)
    {
        if (!is.null(names(x))) {
            end_index <- endIndex(x, all.names=all.names)
            if (length(end_index) == 0)
                return(integer(0))
            return(sapply(end_index, length))
        }
        if (!missing(all.names))
            warning("'all.names' is ignored when patterns have no names")
        end_index <- endIndex(x)
        if (length(end_index) == 0)
            return(integer(0))
        sapply(end_index, length)
    }
)

### Return a single integer or string (not NA).
setMethod("[[", "MIndex",
    function(x, i, j, ...)
    {
        # 'x' is guaranteed to be an MIndex object (if it's not, then
        # the method dispatch algo would not have called this method in the
        # first place), so nargs() is guaranteed to be >= 1
        if (nargs() >= 3)
            stop("too many subscripts")
        subscripts <- list(...)
        if (!missing(i))
            subscripts$i <- i
        if (!missing(j))
            subscripts$j <- j
        # At this point, 'subscripts' should be guaranteed
        # to be of length <= 1
        if (length(subscripts) == 0)
            stop("no index specified")
        key <- subscripts[[1]]
        if (!is.character(key) && !is.numeric(key))
            stop("wrong argument for subsetting an object of class \"MIndex\"")
        if (length(key) < 1)
            stop("attempt to select less than one element")
        if (length(key) > 1)
            stop("attempt to select more than one element")
        if (is.na(key))
            stop("subsetting an object of class \"MIndex\" with NA is not supported")
        if (is.numeric(key)) {
            if (!is.integer(key))
                key <- as.integer(key)
            if (key < 1 || length(x) < key)
                stop("subscript out of bounds")
        }
        key
    }
)

setMethod("$", "MIndex", function(x, name) x[[name]])


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "ByPos_MIndex" class.
### 
### When 'pdict' is a PDict object with no names then matchPDict(pdict, ...)
### returns the matches in a ByPos_MIndex object.
###
### Note that in normal operations the user NEVER needs to create a
### ByPos_MIndex object explicitely or to modify an existing one:
### ByPos_MIndex objects are created by the matchPDict() function
### and have a read-only semantic.
###
### Slot description:
###
###   ends: a list of integer vectors.
###
###   width: temporary hack. In the future we will probably want to store the
###       starts of the matches when the original dictionary has a variable
###       width. Another solution would be to keep the width slot and to make
###       it the same length as the ends slot (it's currently of length 1
###       only).
###

setClass("ByPos_MIndex",
    contains="MIndex",
    representation(
        dups0="Dups",
        ends="list",
        width="integer"
    )
)

setMethod("length", "ByPos_MIndex", function(x) length(x@ends))

setMethod("names", "ByPos_MIndex", function(x) NULL)

setMethod("show", "ByPos_MIndex",
    function(object)
    {
        cat(length(object), "-pattern MIndex object\n", sep="")
    }
)

setMethod("[[", "ByPos_MIndex",
    function(x, i, j, ...)
    {
        key <- callNextMethod()
        if (is.character(key))
            stop("MIndex object has no names")
        if (length(x@dups0) != 0 && !is.na(key2 <- x@dups0@dup2unq[key]))
            key <- key2
        ans_end <- x@ends[[key]]
        ans_width <- rep.int(x@width, length(ans_end))
        ans_start <- ans_end - ans_width + 1L
        new2("IRanges", start=ans_start, width=ans_width, check=FALSE)
    }
)

### An example of a ByPos_MIndex object of length 5 where only the
### 2nd and 5th pattern have matches:
###   > dups0 <- Biostrings:::Dups(c(NA, NA, NA, NA, 2L))
###   > ends <- vector(mode="list", length=5)
###   > ends[[2]] <- c(199L, 402L)
###   > mindex <- new("ByPos_MIndex", dups0=dups0, ends=ends, width=10L)
###   > mindex[[1]]
###   > mindex[[2]]
###   > mindex[[6]] # Error in mindex[[6]] : subscript out of bounds
###   > startIndex(mindex)
###   > endIndex(mindex)
###   > countIndex(mindex)
###
setMethod("startIndex", "ByPos_MIndex",
    function(x, all.names=FALSE)
    {
        if (!missing(all.names))
            warning("'all.names' is ignored when MIndex object has no names")
        .Call("ByPos_MIndex_endIndex", x@dups0@dup2unq, x@ends, 1L - x@width, PACKAGE="Biostrings")
    }
)
setMethod("endIndex", "ByPos_MIndex",
    function(x, all.names=FALSE)
    {
        if (!missing(all.names))
            warning("'all.names' is ignored when MIndex object has no names")
        .Call("ByPos_MIndex_endIndex", x@dups0@dup2unq, x@ends, 0L, PACKAGE="Biostrings")
    }
)

ByPos_MIndex.combine <- function(mi_list, ans_width)
{
    ends_listlist <- lapply(mi_list, function(mi) mi@ends)
    #mergeEnds <- function(...)
    #{
    #    ans <- unlist(list(...), recursive=FALSE, use.names=FALSE)
    #    if (is.null(ans))
    #        return(NULL)
    #    sort(unique(ans))
    #}
    #args <- c(list(FUN=mergeEnds), ends_listlist, list(SIMPLIFY=FALSE))
    #ans_ends <- do.call(mapply, args)
    ans_ends <- .Call("ByPos_MIndex_combine",
                      ends_listlist,
                      PACKAGE="Biostrings")
    new("ByPos_MIndex", ends=ans_ends, width=ans_width)
}

setMethod("combine", signature(x="ByPos_MIndex", y="ByPos_MIndex"),
    function(x, y, ...)
    {
        stop("not ready yet, sorry!")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "ByName_MIndex" class.
### 
### When 'pdict' is a PDict object with no names then matchPDict(pdict, ...)
### returns the matches in a ByName_MIndex object.
###
### Note that in normal operations the user NEVER needs to create a
### ByName_MIndex object explicitely or to modify an existing one:
### ByName_MIndex objects are created by the matchPDict() function
### and have a read-only semantic.
###
### Slot description:
###
###   length: the length of the original dictionary.
###
###   ends_envir: a key-value list (environment) where the values are integer
###       vectors containing the ending positions of the pattern whose
###       position in the original dictionary is given by the key (the keys are
###       strings representing positive integers).
###       
###   width: temporary hack. In the future we will probably want to store the
###       starts of the matches when the original dictionary has a variable
###       width. Another solution would be to keep the width slot and to make
###       it the same length as the ends_envir slot (it's currently of length
###       1 only).
###
###   names: a character vector containing the _unique_ pattern names.
###

setClass("ByName_MIndex",
    contains="MIndex",
    representation(
        length="integer",
        ends_envir="environment",
        width="integer",
        NAMES="character" # R doesn't like @names !!
    )
)

setMethod("length", "ByName_MIndex", function(x) x@length)

setMethod("names", "ByName_MIndex", function(x) x@NAMES)

setMethod("show", "ByName_MIndex",
    function(object)
    {
        cat(length(object), "-pattern sparse MIndex object\n", sep="")
    }
)

setMethod("[[", "ByName_MIndex",
    function(x, i, j, ...)
    {
        key <- callNextMethod()
        if (is.character(key)) {
            pos <- match(key, names(x))
            if (is.na(pos))
                stop("pattern name \"", key, "\" not found")
            key <- pos
        } 
        ans_end <- x@ends_envir[[formatC(key, width=10, format="d", flag="0")]]
        if (is.null(ans_end))
            ans_end <- integer(0)
        ans_width <- rep.int(x@width, length(ans_end))
        ans_start <- ans_end - ans_width + 1L
        new2("IRanges", start=ans_start, width=ans_width, check=FALSE)
    }
)

### An example of a ByName_MIndex object of length 5 where only the
### 2nd pattern has matches:
###   > ends_envir <- new.env(hash=TRUE, parent=emptyenv())
###   > ends_envir[['0000000002']] <- c(199L, 402L)
###   > mindex <- new("ByName_MIndex", length=5L, ends_envir=ends_envir, width=10L, NAMES=letters[1:5])
###   > mindex[[1]]
###   > mindex[[2]]
###   > mindex[[6]] # Error in mindex[[6]] : subscript out of bounds
###   > names(mindex)
###   > mindex[["a"]]
###   > mindex[["b"]]
###   > mindex[["aa"]] # Error in mindex[["aa"]] : pattern name ‘aa’ not found
###   > startIndex(mindex)
###   > startIndex(mindex, all.names=TRUE)
###   > endIndex(mindex)
###   > endIndex(mindex, all.names=TRUE)
###   > countIndex(mindex)
###   > countIndex(mindex, all.names=TRUE)
###
setMethod("startIndex", "ByName_MIndex",
    function(x, all.names=FALSE)
    {
        if (!isTRUEorFALSE(all.names))
            stop("'all.names' must be 'TRUE' or 'FALSE'")
        .Call("ByName_MIndex_endIndex", x@ends_envir, 1L - x@width, x@NAMES, all.names, PACKAGE="Biostrings")
    }
)
setMethod("endIndex", "ByName_MIndex",
    function(x, all.names=FALSE)
    {
        if (!isTRUEorFALSE(all.names))
            stop("'all.names' must be 'TRUE' or 'FALSE'")
        .Call("ByName_MIndex_endIndex", x@ends_envir, 0L, x@NAMES, all.names, PACKAGE="Biostrings")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "unlist" method and the "extractAllMatches" function.
###

setMethod("unlist", "MIndex",
    function(x, recursive=TRUE, use.names=TRUE)
    {
        use.names <- normargUseNames(use.names)
        start_index <- startIndex(x)
        if (length(start_index) == 0)
            ans_start <- integer(0)
        else
            ans_start <- unlist(start_index, recursive=FALSE, use.names=FALSE)
        end_index <- endIndex(x)
        if (length(end_index) == 0)
            ans_end <- integer(0)
        else
            ans_end <- unlist(end_index, recursive=FALSE, use.names=FALSE)
        if (use.names) {
            ans_names <- names(end_index)
            if (!is.null(ans_names))
                ans_names <- rep.int(ans_names, times=sapply(end_index, length))
        } else {
            ans_names <- NULL
        }
        ans_width <- ans_end - ans_start + 1L
        ans <- new2("IRanges", start=ans_start, width=ans_width, check=FALSE)
        names(ans) <- ans_names
        ans
    }
)

extractAllMatches <- function(subject, mindex)
{
    if (is(subject, "MaskedXString"))
        subject <- unmasked(subject)
    else if (!is(subject, "XString"))
        stop("'subject' must be an XString or MaskedXString object")
    if (!is(mindex, "MIndex"))
        stop("'mindex' must be an MIndex object")
    if (is.null(names(mindex)))
        stop("extractAllMatches() works only with a \"MIndex\" object with names")
    allviews <- unlist(mindex)
    ans <- unsafe.newXStringViews(subject, start(allviews), width(allviews))
    names(ans) <- names(allviews)
    ans
}

