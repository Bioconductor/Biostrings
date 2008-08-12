### =========================================================================
### Dups objects
### -------------------------------------------------------------------------
###
### The Dups class is a general container for storing the direct and reverse
### mappings between duplicated and unique elements of an arbitrary vector.
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

### Returns the unq2dup map, not the full Dups object!
Dups.diff <- function(x, y)
{
    .Call("Dups_diff", x@unq2dup, y@dup2unq, PACKAGE="Biostrings")
}

