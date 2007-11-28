### =========================================================================
### The findPalindromes() generic & related functions
### -------------------------------------------------------------------------


debug_find_palindromes <- function()
{
    invisible(.Call("find_palindromes_debug", PACKAGE="Biostrings"))
}

### Return a list with the "start" and the "end" components.
.find.palindromes <- function(subject, min.armlength, max.ngaps, L2R_lkup)
{
    
    views <- .Call("find_palindromes",
                   subject@data@xp, subject@offset, subject@length,
                   min.armlength, max.ngaps, L2R_lkup,
                   PACKAGE="Biostrings")
    new("BStringViews", subject=subject, views=data.frame(start=views$start, end=views$end))
}

setGeneric(
    "findPalindromes", signature="subject",
    function(subject, min.armlength=4, max.ngaps=1, anti)
        standardGeneric("findPalindromes")
)

setMethod("findPalindromes", "BString",
    function(subject, min.armlength, max.ngaps, anti)
    {
        if (!is.numeric(min.armlength) || length(min.armlength) != 1 || is.na(min.armlength))
            stop("'min.armlength' must be a single integer")
        min.armlength <- as.integer(min.armlength)
        if (min.armlength < 2)
            stop("'min.armlength' must be >= 2")
        if (!is.numeric(max.ngaps) || length(max.ngaps) != 1 || is.na(max.ngaps))
            stop("'max.ngaps' must be a single integer")
        max.ngaps <- as.integer(max.ngaps)
        if (max.ngaps < 0)
            stop("'max.ngaps' must be a non-negative integer")
        L2R_lkup <- NULL
        if (class(subject) %in% c("DNAString", "RNAString")) {
            if (missing(anti) || anti)
                L2R_lkup <- getDNAComplementLookup()
        } else if (!missing(anti)) {
            stop("you can specify 'anti' only for a DNAString or RNAString subject")
        }
        .find.palindromes(subject, min.armlength, max.ngaps, L2R_lkup)
    }
)

