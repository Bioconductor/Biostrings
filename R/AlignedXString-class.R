### =========================================================================
### AlignedXString objects
### -------------------------------------------------------------------------
### An AlignedXString object contains an alignment.


setClass("AlignedXString",
    representation(
        unaligned="XString",
        quality="XString",
        range="IRanges",
        inserts="IRanges",
        profile="XNumeric"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Initialization.
###

setMethod("initialize", "AlignedXString",
    function(.Object, unaligned, quality, range, inserts, profile, check = TRUE)
    {
        slot(.Object, "unaligned", check = check) <- unaligned
        slot(.Object, "quality", check = check) <- quality
        slot(.Object, "range", check = check) <- range
        slot(.Object, "inserts", check = check) <- inserts
        slot(.Object, "profile", check = check) <- profile
        .Object
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.AlignedXString <- function(object)
{
    message <- character(0)
    if (length(message) == 0)
        message <- NULL
    message
}

setValidity("AlignedXString",
    function(object)
    {
        problems <- .valid.AlignedXString(object)
        if (is.null(problems)) TRUE else problems
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods.
###

setGeneric("unaligned", function(x) standardGeneric("unaligned"))
setMethod("unaligned", "AlignedXString", function(x) x@unaligned)

setGeneric("aligned", function(x) standardGeneric("aligned"))
setMethod("aligned", "AlignedXString",
          function(x) {
            if (width(x) == 0) {
              value <- XString("")
            } else {
              string <- subXString(unaligned(x), start(x), end(x))
              startInserts <- start(x@inserts)
              endInserts <- end(x@inserts)
              widthInserts <- width(x@inserts)
              ncharString <- nchar(string)
              if (length(startInserts) == 0) {
                value <- string
              } else {
                gapStrings <- c("", sapply(widthInserts, function(x) paste(rep("-", x), collapse = "")))
                subdividedString <-
                  substring(as.character(string), c(1, startInserts), c(startInserts - 1, ncharString))
                value <-
                  XString(class(unaligned(x)), paste(gapStrings, subdividedString, sep = "", collapse = ""))
              }
            }
            value
          })

setMethod("start", "AlignedXString", function(x) if (width(x) == 0) NA else start(x@range))
setMethod("end", "AlignedXString", function(x) if (width(x) == 0) NA else end(x@range))
setMethod("width", "AlignedXString", function(x) width(x@range))
setMethod("profile", "AlignedXString", function(fitted, ...) as.numeric(fitted@profile))
setMethod("length", "AlignedXString", function(x) width(x@range) + sum(width(x@inserts)))
setMethod("nchar", "AlignedXString", function(x, type="chars", allowNA=FALSE) length(x))
setMethod("alphabet", "AlignedXString", function(x) alphabet(unaligned(x)))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method.
###

### TODO: Make the "show" method to format the alignment in a SGD fashion
### i.e. split in 60-letter blocks and use the "|" character to highlight
### exact matches.
setMethod("show", "AlignedXString",
    function(object)
    {
        if (width(object) == 0)
          cat("[1] \"\"\n")
        else
          cat(paste("[", start(object), "]", sep = ""),
              toSeqSnippet(aligned(object), getOption("width") - 8), "\n")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "as.character" method.
###

setMethod("as.character", "AlignedXString",
    function(x)
    {
        as.character(aligned(x))
    }
)

setMethod("toString", "AlignedXString", function(x, ...) as.character(x))
