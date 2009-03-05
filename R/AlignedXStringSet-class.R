### =========================================================================
### AlignedXStringSet objects
### -------------------------------------------------------------------------
### An AlignedXStringSet object contains an alignment.


setClass("AlignedXStringSet",
    representation(
        unaligned="XStringSet",
        range="IRanges",
        mismatch = "IntegerList",
        indel="IRangesList"
    )
)

setClass("QualityAlignedXStringSet", contains="AlignedXStringSet")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Initialization.
###

setMethod("initialize", "AlignedXStringSet",
    function(.Object, unaligned, range, mismatch, indel, check = TRUE)
    {
        if (is(unaligned, "QualityScaledXStringSet"))
            stop("'unaligned' must not be of class 'QualityScaledXStringSet'")
        slot(.Object, "unaligned", check = check) <- unaligned
        slot(.Object, "range", check = check) <- range
        slot(.Object, "mismatch", check = check) <- mismatch
        slot(.Object, "indel", check = check) <- indel
        .Object
    }
)

setMethod("initialize", "QualityAlignedXStringSet",
    function(.Object, unaligned, range, mismatch, indel, check = TRUE)
    {
        if (!is(unaligned, "QualityScaledXStringSet"))
            stop("'unaligned' must be of class 'QualityScaledXStringSet'")
        slot(.Object, "unaligned", check = check) <- unaligned
        slot(.Object, "range", check = check) <- range
        slot(.Object, "mismatch", check = check) <- mismatch
        slot(.Object, "indel", check = check) <- indel
        .Object
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.AlignedXStringSet <- function(object)
{
    message <- character(0)
    if (length(object@range) != length(mismatch(object)))
        message <- c(message, "length(range) != length(mismatch)")
    if (length(mismatch(object)) != length(indel(object)))
        message <- c(message, "length(mismatch) != length(indel)")
    if (!(length(object@unaligned) %in% c(1, length(object@range))))
        message <- c(message, "length(unaligned) != 1 or length(range)")
    if (is(object@unaligned, "QualityScaledXStringSet"))
        message <- c(message, "'unaligned' must not be of class 'QualityScaledXStringSet'")
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

.valid.QualityAlignedXStringSet <- function(object)
{
    message <- character(0)
    if (length(object@range) != length(mismatch(object)))
        message <- c(message, "length(range) != length(mismatch)")
    if (length(mismatch(object)) != length(indel(object)))
        message <- c(message, "length(mismatch) != length(indel)")
    if (!(length(object@unaligned) %in% c(1, length(object@range))))
        message <- c(message, "length(unaligned) != 1 or length(range)")
    if (!(length(object@quality) %in% c(1, length(object@range))))
        message <- c(message, "length(quality) != 1 or length(range)")
    if (!is(object@unaligned, "QualityScaledXStringSet"))
        message <- c(message, "'unaligned' must be of class 'QualityScaledXStringSet'")
    if (length(message) == 0)
        message <- NULL
    message
}

setValidity("QualityAlignedXStringSet",
    function(object)
    {
        problems <- .valid.QualityAlignedXStringSet(object)
        if (is.null(problems)) TRUE else problems
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods.
###

setGeneric("unaligned", function(x) standardGeneric("unaligned"))
setMethod("unaligned", "AlignedXStringSet", function(x) x@unaligned)

setGeneric("aligned", function(x, ...) standardGeneric("aligned"))
setMethod("aligned", "AlignedXStringSet",
          function(x, degap = FALSE) {
              if (degap) {
                  if (length(unaligned(x)) == 1) {
                      value <-
                        Views(unaligned(x)[[1]], start=start(x@range), end=end(x@range))
                  } else {
                      value <-
                        narrow(as(unaligned(x), "XStringSet"), start=start(x@range), end=end(x@range))
                  }
              } else {
                  codecX <- xscodec(x)
                  if (is.null(codecX)) {
                      gapCode <- charToRaw("-")
                  } else {
                      letters2codes <- codecX@codes
                      names(letters2codes) <- codecX@letters
                      gapCode <- as.raw(letters2codes[["-"]])
                  }
                  value <- 
                    .Call("AlignedXStringSet_align_aligned", x, gapCode, PACKAGE="Biostrings")
              }
              value
          })

setMethod("start", "AlignedXStringSet", function(x) start(x@range))
setMethod("end", "AlignedXStringSet", function(x) end(x@range))
setMethod("width", "AlignedXStringSet", function(x) width(x@range))
setGeneric("indel", function(x) standardGeneric("indel"))
setMethod("indel", "AlignedXStringSet", function(x) x@indel)
setGeneric("nindel", function(x) standardGeneric("nindel"))
setMethod("nindel", "AlignedXStringSet", function(x) summary(indel(x)))
setMethod("length", "AlignedXStringSet", function(x) length(x@range))
setMethod("nchar", "AlignedXStringSet",
          function(x, type="chars", allowNA=FALSE) .Call("AlignedXStringSet_nchar", x, PACKAGE="Biostrings"))
setMethod("xsbasetype", "AlignedXStringSet", function(x) xsbasetype(unaligned(x)))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method.
###

### TODO: Make the "show" method to format the alignment in a SGD fashion
### i.e. split in 60-letter blocks and use the "|" character to highlight
### exact matches.
setMethod("show", "AlignedXStringSet",
    function(object)
    {
        if (length(object) == 0)
            cat("Empty ", class(object), "\n", sep = "")
        else {
            if (length(object) > 1)
                cat(class(object), " (1 of ", length(object), ")\n", sep = "")
            if (width(object)[1] == 0)
                cat("[1] \"\"\n")
            else
                cat(paste("[", start(object)[1], "]", sep = ""),
                    toSeqSnippet(aligned(object)[[1]], getOption("width") - 8), "\n")
        }
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

setMethod("toString", "AlignedXStringSet", function(x, ...) toString(as.character(x), ...))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

setMethod("[", "AlignedXStringSet",
    function(x, i, j, ..., drop)
    {
        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")
        if (missing(i) || (is.logical(i) && all(i)))
            return(x)
        if (is.logical(i))
            i <- which(i)
        if (!is.numeric(i) || any(is.na(i)))
            stop("invalid subsetting")
        if (any(i < 1) || any(i > length(x)))
            stop("subscript out of bounds")
        new("AlignedXStringSet",
            unaligned = .safe.subset.XStringSet(x@unaligned, i),
            range = x@range[i,,drop = FALSE],
            mismatch = x@mismatch[i], indel = x@indel[i])
    }
)

setMethod("[", "QualityAlignedXStringSet",
    function(x, i, j, ..., drop)
    {
        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")
        if (missing(i) || (is.logical(i) && all(i)))
            return(x)
        if (is.logical(i))
            i <- which(i)
        if (!is.numeric(i) || any(is.na(i)))
            stop("invalid subsetting")
        if (any(i < 1) || any(i > length(x)))
            stop("subscript out of bounds")
        new("QualityAlignedXStringSet",
            unaligned = .safe.subset.XStringSet(x@unaligned, i),
            range = x@range[i,,drop = FALSE],
            mismatch = x@mismatch[i], indel = x@indel[i])
    }
)

setReplaceMethod("[", "AlignedXStringSet",
    function(x, i, j,..., value)
    {
        stop("attempt to modify the value of a ", class(x), " instance")
    }
)

setMethod("rep", "AlignedXStringSet",
    function(x, times)
		x[rep.int(seq_len(length(x)), times)]
)
