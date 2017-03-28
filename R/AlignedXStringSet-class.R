### =========================================================================
### AlignedXStringSet objects
### -------------------------------------------------------------------------
### An AlignedXStringSet object contains an alignment.
###
### NOTE: Because the 'unaligned' slot of an AlignedXStringSet object
### must not be a QualityScaledXStringSet object (see Validity below), then
### the QualityAlignedXStringSet class cannot contain the AlignedXStringSet
### class. Otherwise, any QualityAlignedXStringSet object would be invalid!

setClass("AlignedXStringSet0",
    contains="Vector",
    representation(
        "VIRTUAL",
        unaligned="XStringSet",
        range="IRanges",
        mismatch="CompressedIntegerList",
        indel="CompressedIRangesList"
    )
)

### Combine the new parallel slots with those of the parent class. Make sure
### to put the new parallel slots *first*.
setMethod("parallelSlotNames", "AlignedXStringSet0",
    function(x) c("unaligned", "range", "mismatch", "indel", callNextMethod())
)

setClass("AlignedXStringSet", contains="AlignedXStringSet0")

setClass("QualityAlignedXStringSet",
    contains="AlignedXStringSet0",
    representation(
        unaligned="QualityScaledXStringSet"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.AlignedXStringSet0 <- function(object)
{
    message <- NULL
    if (length(object@range) != length(mismatch(object)))
        message <- c(message, "length(range) != length(mismatch)")
    if (length(mismatch(object)) != length(indel(object)))
        message <- c(message, "length(mismatch) != length(indel)")
    if (!(length(object@unaligned) %in% c(1, length(object@range))))
        message <- c(message, "length(unaligned) != 1 or length(range)")
    message
}

setValidity("AlignedXStringSet0",
    function(object)
    {
        problems <- .valid.AlignedXStringSet0(object)
        if (is.null(problems)) TRUE else problems
    }
)

.valid.AlignedXStringSet <- function(object)
{
    message <- NULL
    if (is(object@unaligned, "QualityScaledXStringSet"))
        message <- c(message, "'unaligned' must not be a QualityScaledXStringSet object")
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
    message <- NULL
    ## FIXME: surely something different is meant here because the
    ## QualityAlignedXStringSet class has no 'quality' slot!
    #if (!(length(object@quality) %in% c(1, length(object@range))))
    #    message <- c(message, "length(quality) != 1 or length(range)")
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
setMethod("unaligned", "AlignedXStringSet0", function(x) x@unaligned)

setGeneric("aligned", function(x, ...) standardGeneric("aligned"))
setMethod("aligned", "AlignedXStringSet0",
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
                    .Call2("AlignedXStringSet_align_aligned", x, gapCode, PACKAGE="Biostrings")
              }
              value
          })

setMethod("start", "AlignedXStringSet0", function(x) start(x@range))
setMethod("end", "AlignedXStringSet0", function(x) end(x@range))
setMethod("width", "AlignedXStringSet0", function(x) width(x@range))
setMethod("ranges", "AlignedXStringSet0", function(x) x@range)
setGeneric("indel", function(x) standardGeneric("indel"))
setMethod("indel", "AlignedXStringSet0", function(x) x@indel)
setGeneric("nindel", function(x) standardGeneric("nindel"))
setMethod("nindel", "AlignedXStringSet0", function(x) summary(indel(x)))
setMethod("length", "AlignedXStringSet0", function(x) length(x@range))
setMethod("nchar", "AlignedXStringSet0",
          function(x, type="chars", allowNA=FALSE) .Call2("AlignedXStringSet_nchar", x, PACKAGE="Biostrings"))
setMethod("seqtype", "AlignedXStringSet0", function(x) seqtype(unaligned(x)))

setMethod("parallelVectorNames", "AlignedXStringSet0",
          function(x) c("unaligned", "range", "mismatch", "indel",
                        "start", "end", "width", "nindel"))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method.
###

### TODO: Make the "show" method to format the alignment in a SGD fashion
### i.e. split in 60-letter blocks and use the "|" character to highlight
### exact matches.
setMethod("show", "AlignedXStringSet0",
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
### Coercion.
###

setMethod("as.character", "AlignedXStringSet0",
    function(x)
    {
        as.character(aligned(x))
    }
)

setMethod("toString", "AlignedXStringSet0",
    function(x, ...) toString(as.character(x), ...)
)

