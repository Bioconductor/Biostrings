### =========================================================================
### Element-wise comparison of XStringSet objects
### -------------------------------------------------------------------------
###


### Method signatures for binary comparison operators.
.OP2_SIGNATURES <- list(
    c("XStringSet", "XStringSet"),
    c("XStringSet", "Vector"),
    c("XStringSet", "vector"),
    c("Vector", "XStringSet"),
    c("vector", "XStringSet")
)

.coerce_to_XStringSet_and_call_next_method <- function(f, x, y, ...,
                                                       what_op=NULL)
{
    if (is.null(what_op))
        what_op <- paste0(f, "()")
    seqtypes <- get_seqtype_switches_before_binary_op(x, y, what_op=what_op)
    class1 <- paste0(seqtypes[[1L]], "StringSet")
    class2 <- paste0(seqtypes[[2L]], "StringSet")
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
### pcompare()
###

.pcompare_XStringSet <- function(x, y)
    .coerce_to_XStringSet_and_call_next_method("pcompare", x, y,
                                               what_op="comparison")

setMethods("pcompare", .OP2_SIGNATURES, .pcompare_XStringSet)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### match()
###

.XStringSet.match <- function(x, table,
                              nomatch=NA_integer_, incomparables=NULL)
{
    .coerce_to_XStringSet_and_call_next_method("match", x, table,
                              nomatch=nomatch, incomparables=incomparables)
}

setMethods("match", .OP2_SIGNATURES, .XStringSet.match)

