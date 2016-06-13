### =========================================================================
### Comparing and ordering the elements in one or more XStringSet objects
### -------------------------------------------------------------------------
###


### Returns a character vector of length 2 containing the 2 XStringSet direct
### concrete subclasses that 'x' and 'y' need to be coerced to before they can
### actually be compared.
.coerce_to <- function(x, y)
{
    seqtype1 <- try(seqtype(x), silent=TRUE)
    if (is(seqtype1, "try-error"))
        seqtype1 <- "B"
    seqtype2 <- try(seqtype(y), silent=TRUE)
    if (is(seqtype2, "try-error"))
        seqtype2 <- "B"
    if (seqtype1 != seqtype2) {
        if ((seqtype1 != "B" && seqtype2 == "AA")
         || (seqtype2 != "B" && seqtype1 == "AA"))
            stop("comparison between a \"", class(x), "\" instance ",
                 "and a \"", class(y), "\" instance\n",
                 "  is not supported")
        if (seqtype1 == "B" && seqtype2 != "AA")
            seqtype1 <- seqtype2
        if (seqtype2 == "B" && seqtype1 != "AA")
            seqtype2 <- seqtype1
    }
    class1 <- paste0(seqtype1, "StringSet")
    class2 <- paste0(seqtype2, "StringSet")
    c(class1, class2)
}

.coerce_and_call_next_method <- function(f, x, y, ...)
{
    classes <- .coerce_to(x, y)
    class1 <- classes[[1L]]
    class2 <- classes[[2L]]
    if (!is(x, class1))
        x <- as(x, class1)
    if (!is(y, class2))
        y <- as(y, class2)
    ## We cannot use callNextMethod() in this context (only from within the
    ## body of a method definition), so we use getMethod() instead.
    XRawList_method <- getMethod(f, c("XRawList", "XRawList"))
    XRawList_method(x, y, ...)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### pcompare().
###

### Method signatures for binary comparison operators.
.OP2_SIGNATURES <- list(
    c("XStringSet", "XStringSet"),
    c("XStringSet", "Vector"),
    c("XStringSet", "vector"),
    c("Vector", "XStringSet"),
    c("vector", "XStringSet")
)

.pcompare_XStringSet <- function(x, y)
    .coerce_and_call_next_method("pcompare", x, y)

setMethods("pcompare", .OP2_SIGNATURES, .pcompare_XStringSet)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### match().
###

.XStringSet.match <- function(x, table,
                              nomatch=NA_integer_, incomparables=NULL)
{
    .coerce_and_call_next_method("match", x, table,
                                 nomatch=nomatch, incomparables=incomparables)
}

setMethods("match", .OP2_SIGNATURES, .XStringSet.match)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### is.na() and related methods
###

setMethod("is.na", "XStringSet", function(x) rep(FALSE, length(x)))

setMethod("anyNA", "XStringSet", function(x, recursive=FALSE) FALSE)
