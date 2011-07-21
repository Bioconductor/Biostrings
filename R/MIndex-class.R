### =========================================================================
### MIndex objects
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "MIndex" VIRTUAL class.
###
### This class serves as a base class for a family of match containers that
### can be manipulated thru a common API. This common API is the MIndex API
### defined in this section.
### Note that, in normal operations, the user should never need to create
### MIndex objects directly or to modify existing ones. Those objects are
### typically returned by a sequence matching/alignment function like
### vmatchPattern() or matchPDict().
###

setClass("MIndex",
    contains="IRangesList",
    representation(
        "VIRTUAL",
        width0="integer",
        NAMES="characterORNULL"
    )
)

setMethod("length", "MIndex", function(x) length(x@width0))

setGeneric("width0", function(x) standardGeneric("width0"))

setMethod("width0", "MIndex", function(x) x@width0)

setMethod("names", "MIndex", function(x) x@NAMES)

setReplaceMethod("names", "MIndex",
    function(x, value)
        stop("attempt to modify the names of a ", class(x), " instance")
)

setGeneric("startIndex", function(x) standardGeneric("startIndex"))

setGeneric("endIndex", function(x) standardGeneric("endIndex"))

setMethod("elementLengths", "MIndex",
    function(x) elementLengths(endIndex(x))
)

setGeneric("countIndex", function(x) standardGeneric("countIndex"))

setMethod("countIndex", "MIndex", function(x) elementLengths(x))

setMethod("unlist", "MIndex",
    function(x, recursive=TRUE, use.names=TRUE)
    {
        use.names <- normargUseNames(use.names)
        start_index <- startIndex(x)
        ans_start <- unlist(start_index, recursive=FALSE, use.names=FALSE)
        if (is.null(ans_start))
            ans_start <- integer(0)
        end_index <- endIndex(x)
        ans_end <- unlist(end_index, recursive=FALSE, use.names=FALSE)
        if (is.null(ans_end))
            ans_end <- integer(0)
        if (use.names) {
            ans_names <- names(x)
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

setAs("MIndex", "CompressedIRangesList",
    function(from)
    {
        IRanges:::newCompressedList("CompressedIRangesList",
                                    unlistData = unlist(from),
                                    end = cumsum(countIndex(from)),
                                    NAMES = names(from))
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
    allviews <- unlist(mindex)
    ans <- unsafe.newXStringViews(subject, start(allviews), width(allviews))
    names(ans) <- names(allviews)
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "ByPos_MIndex" class.
### 

setClass("ByPos_MIndex",
    contains="MIndex",
    representation(
        dups0="Dups",
        ends="list"  # same length as the "width0" slot
    )
)

setMethod("show", "ByPos_MIndex",
    function(object)
    {
        cat("MIndex object of length ", length(object), "\n", sep="")
    }
)

setMethod("[[", "ByPos_MIndex",
    function(x, i, j, ...)
    {
        i <- IRanges:::checkAndTranslateDbleBracketSubscript(x, i)
        if (length(x@dups0) != 0 && !is.na(i2 <- high2low(x@dups0)[i]))
            i <- i2
        ans_end <- x@ends[[i]]
        if (is.null(ans_end))
            ans_end <- integer(0)
        ans_width <- rep.int(x@width0[i], length(ans_end))
        ans_start <- ans_end - x@width0[i] + 1L
        new2("IRanges", start=ans_start, width=ans_width, check=FALSE)
    }
)

### An example of a ByPos_MIndex object of length 5 where only the
### 2nd and 5th pattern have matches:
###   > width0 <- c(9L, 10L, 8L, 4L, 10L)
###   > dups0 <- Dups(c(NA, NA, NA, NA, 2))
###   > ends <- vector(mode="list", length=5)
###   > ends[[2]] <- c(199L, 402L)
###   > mindex <- new("ByPos_MIndex", width0=width0, NAMES=letters[1:5], dups0=dups0, ends=ends)
###   > mindex[[1]]
###   > mindex[[2]]
###   > mindex[[6]] # Error in mindex[[6]] : subscript out of bounds
###   > startIndex(mindex)
###   > endIndex(mindex)
###   > countIndex(mindex)
###
setMethod("startIndex", "ByPos_MIndex",
    function(x)
    {
        .Call2("ByPos_MIndex_endIndex",
              high2low(x@dups0), x@ends, x@width0,
              PACKAGE="Biostrings")
    }
)
setMethod("endIndex", "ByPos_MIndex",
    function(x)
    {
        .Call2("ByPos_MIndex_endIndex",
              high2low(x@dups0), x@ends, NULL,
              PACKAGE="Biostrings")
    }
)

ByPos_MIndex.combine <- function(mi_list)
{
    if (length(mi_list) == 0)
        stop("cannot combine an empty list of MIndex objects")
    ans_width0 <- mi_list[[1]]@width0  # all 'mi_list[[i]]@width0' should be the same!
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
    ans_ends <- .Call2("ByPos_MIndex_combine",
                      ends_listlist,
                      PACKAGE="Biostrings")
    new("ByPos_MIndex", width0=ans_width0, ends=ans_ends)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "SparseMIndex" class (DISABLED FOR NOW).
### 
### Slot description:
###
###   ends_envir: a key-value list (environment) where the values are integer
###       vectors containing the ending positions of the pattern whose
###       position in the original dictionary is given by the key (the keys are
###       strings representing positive integers).
### 

if (FALSE) {

  setClass("SparseMIndex",
    contains="MIndex",
    representation(
        ends_envir="environment"
    )
  )

  setMethod("show", "SparseMIndex",
    function(object)
    {
        cat("Sparse MIndex object of length ", length(object), "\n", sep="")
    }
  )

  setMethod("[[", "SparseMIndex",
    function(x, i, j, ...)
    {
        i <- IRanges:::checkAndTranslateDbleBracketSubscript(x, i)
        ans_end <- x@ends_envir[[formatC(i, width=10, format="d", flag="0")]]
        if (is.null(ans_end))
            ans_end <- integer(0)
        ans_width <- rep.int(x@width0[i], length(ans_end))
        ans_start <- ans_end - x@width0[i] + 1L
        new2("IRanges", start=ans_start, width=ans_width, check=FALSE)
    }
  )

  ### An example of a SparseMIndex object of length 5 where only the
  ### 2nd pattern has matches:
  ###   > width0 <- c(9L, 10L, 8L, 4L, 10L)
  ###   > ends_envir <- new.env(hash=TRUE, parent=emptyenv())
  ###   > ends_envir[['0000000002']] <- c(199L, 402L)
  ###   > mindex <- new("SparseMIndex", width0=width0, NAMES=letters[1:5], ends_envir=ends_envir)
  ###   > mindex[[1]]
  ###   > mindex[[2]]
  ###   > mindex[[6]] # Error in mindex[[6]] : subscript out of bounds
  ###   > names(mindex)
  ###   > mindex[["a"]]
  ###   > mindex[["b"]]
  ###   > mindex[["aa"]] # Error in mindex[["aa"]] : pattern name ‘aa’ not found
  ###   > startIndex(mindex)
  ###   > endIndex(mindex)
  ###   > countIndex(mindex)
  ###
  setMethod("startIndex", "SparseMIndex",
    function(x)
    {
        all.names <- TRUE
        .Call2("SparseMIndex_endIndex",
              x@ends_envir, x@width0, x@NAMES, all.names,
              PACKAGE="Biostrings")
    }
  )
  setMethod("endIndex", "SparseMIndex",
    function(x)
    {
        all.names <- TRUE
        .Call2("SparseMIndex_endIndex",
              x@ends_envir, NULL, x@NAMES, all.names,
              PACKAGE="Biostrings")
    }
  )

}
