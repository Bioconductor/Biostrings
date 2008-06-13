### =========================================================================
### AlignedXStringSet objects
### -------------------------------------------------------------------------
### An AlignedXStringSet object contains an alignment.


setClass("AlignedXStringSet",
    representation(
        unaligned="XStringSet",
        quality="XStringSet",
        range="IRanges",
        indels="list"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Initialization.
###

setMethod("initialize", "AlignedXStringSet",
    function(.Object, unaligned, quality, range, indels, check = TRUE)
    {
        slot(.Object, "unaligned", check = check) <- unaligned
        slot(.Object, "quality", check = check) <- quality
        slot(.Object, "range", check = check) <- range
        slot(.Object, "indels", check = check) <- indels
        .Object
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.AlignedXStringSet <- function(object)
{
    message <- character(0)
    if (length(object@range) != length(indels(object)))
        message <- c(message, "length(range) != length(indels)")
    if (!(length(object@unaligned) %in% c(1, length(object@range))))
        message <- c(message, "length(unaligned) != 1 or length(range)")
    if (!(length(object@quality) %in% c(1, length(object@range))))
        message <- c(message, "length(quality) != 1 or length(range)")
    if (length(message) == 0)
        message <- NULL
    message
}

setValidity("AlignedXStringSet",
    function(object)
    {
        problems <- .valid.AlignedXStringSet(object)
        if (is.null(problems)) TRUE else problems
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods.
###

setGeneric("unaligned", function(x) standardGeneric("unaligned"))
setMethod("unaligned", "AlignedXStringSet", function(x) x@unaligned)

setGeneric("aligned", function(x) standardGeneric("aligned"))
setMethod("aligned", "AlignedXStringSet",
          function(x) {
              codecX <- codec(x)
              if (is.null(codecX)) {
                  gapCode <- charToRaw("-")
              } else {
                  letters2codes <- codecX@codes
                  names(letters2codes) <- codecX@letters
                  gapCode <- as.raw(letters2codes[["-"]])
              }
              .Call("AlignedXStringSet_align_aligned", x, gapCode, PACKAGE="Biostrings")
          })

setMethod("start", "AlignedXStringSet", function(x) start(x@range))
setMethod("end", "AlignedXStringSet", function(x) end(x@range))
setMethod("width", "AlignedXStringSet", function(x) width(x@range))
setGeneric("indels", function(x) standardGeneric("indels"))
setMethod("indels", "AlignedXStringSet", function(x) x@indels)
setMethod("length", "AlignedXStringSet", function(x) length(x@range))
setMethod("nchar", "AlignedXStringSet",
          function(x, type="chars", allowNA=FALSE) .Call("AlignedXStringSet_nchar", x, PACKAGE="Biostrings"))
setMethod("alphabet", "AlignedXStringSet", function(x) alphabet(unaligned(x)))
setMethod("codec", "AlignedXStringSet", function(x) codec(unaligned(x)))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method.
###

### TODO: Make the "show" method to format the alignment in a SGD fashion
### i.e. split in 60-letter blocks and use the "|" character to highlight
### exact matches.
setMethod("show", "AlignedXStringSet",
    function(object)
    {
        if (width(object)[1] == 0)
          cat("[1] \"\"\n")
        else
          cat(paste("[", start(object)[1], "]", sep = ""),
              toSeqSnippet(aligned(object)[[1]], getOption("width") - 8), "\n")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "as.character" method.
###

setMethod("as.character", "AlignedXStringSet",
    function(x)
    {
        as.character(aligned(x))
    }
)

setMethod("toString", "AlignedXStringSet", function(x, ...) as.character(x))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

.safe.subset.XStringSet <- function(x, i)
{
    if (length(x) == 1) x else x[i]
}

setMethod("[", "AlignedXStringSet",
    function(x, i, j, ..., drop)
    {
        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")
        if (missing(i))
            return(x)
        if (!is.numeric(i) || any(is.na(i)))
            stop("invalid subsetting")
        if (any(i < 1) || any(i > length(x)))
            stop("subscript out of bounds")
        new("AlignedXStringSet",
            unaligned = .safe.subset.XStringSet(x@unaligned, i),
            quality = .safe.subset.XStringSet(x@quality, i),
            range = x@range[i,,drop = FALSE], indels = x@indels[i])
    }
)

setReplaceMethod("[", "AlignedXStringSet",
    function(x, i, j,..., value)
    {
        stop("attempt to modify the value of a ", class(x), " instance")
    }
)
