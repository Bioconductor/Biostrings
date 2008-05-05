### =========================================================================
### XStringAlign objects
### -------------------------------------------------------------------------
### An XStringAlign object contains the result of the alignment of 2 XString
### objects of the same subtype.


setClass("XStringAlign",
    representation(
        string1="XString",
        string2="XString",
        quality1="XString",
        quality2="XString",
        match1="IRanges",
        match2="IRanges",
        inserts1="IRanges",
        inserts2="IRanges",
        type="character",
        score="numeric",
        constantMatrix="matrix",
        gapOpening="numeric",
        gapExtension="numeric"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Initialization.
###

setMethod("initialize", "XStringAlign",
    function(.Object, string1, string2, quality1, quality2, match1, match2,
             inserts1, inserts2, type, score, constantMatrix, gapOpening,
             gapExtension, check = TRUE)
    {
        if (!identical(class(string1), class(string2)))
            stop("'string1' and 'string2' must be XString objects of the same subtype")
        if (length(type) != 1 || !(type %in% c("global", "local", "overlap")))
            stop("'type' must be one of 'global', 'local', or 'overlap'")
        gapOpening <- as.double(- abs(gapOpening))
        if (is.na(gapOpening) || length(gapOpening) != 1)
            stop("'gapOpening' must be a non-positive numeric vector of length 1")
        gapExtension <- as.double(- abs(gapExtension))
        if (is.na(gapExtension) || length(gapExtension) != 1)
            stop("'gapExtension' must be a non-positive numeric vector of length 1")
        slot(.Object, "string1", check = check) <- string1
        slot(.Object, "string2", check = check) <- string2
        slot(.Object, "quality1", check = check) <- quality1
        slot(.Object, "quality2", check = check) <- quality2
        slot(.Object, "match1", check = check) <- match1
        slot(.Object, "match2", check = check) <- match2
        slot(.Object, "inserts1", check = check) <- inserts1
        slot(.Object, "inserts2", check = check) <- inserts2
        slot(.Object, "type", check = check) <- type
        slot(.Object, "score", check = check) <- score
        slot(.Object, "constantMatrix", check = check) <- constantMatrix
        slot(.Object, "gapOpening", check = check) <- gapOpening
        slot(.Object, "gapExtension", check = check) <- gapExtension
        .Object
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.XStringAlign <- function(object)
{
    message <- character(0)
    if (!identical(class(object@string1), class(object@string2)))
        message <- c(message, "'string1' and 'string2' must be XString objects of the same subtype")
    if (length(object@type) != 1 || !(object@type %in% c("global", "local", "overlap")))
        message <- c(message, "'type' must be one of 'global', 'local', or 'overlap'")
    if (is.na(object@gapOpening) || length(object@gapOpening) != 1)
        message <- c(message, "'gapOpening' must be a non-positive numeric vector of length 1")
    if (is.na(object@gapExtension) || length(object@gapExtension) != 1)
        message <- c(message, "'gapExtension' must be a non-positive numeric vector of length 1")
    if (length(message) == 0)
        message <- NULL
    message
}

setValidity("XStringAlign",
    function(object)
    {
        problems <- .valid.XStringAlign(object)
        if (is.null(problems)) TRUE else problems
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods.
###

.addInserts <-
function(string, inserts)
{
    startInserts <- start(inserts)
    endInserts <- end(inserts)
    widthInserts <- width(inserts)
    ncharString <- nchar(string)
    if (length(startInserts) == 0) {
        value <- string
    } else {
        gapStrings <- c("", sapply(widthInserts, function(x) paste(rep("-", x), collapse = "")))
        subdividedString <-
          substring(as.character(string), c(1, endInserts + 1), c(endInserts, ncharString))
		value <-
          XString(class(string), paste(gapStrings, subdividedString, sep = "", collapse = ""))
    }
    value
}

setGeneric("align1", function(x) standardGeneric("align1"))
setMethod("align1", "XStringAlign",
          function(x) .addInserts(subXString(x@string1, start(x@match1), end(x@match1)),
                                  x@inserts1))

setGeneric("align2", function(x) standardGeneric("align2"))
setMethod("align2", "XStringAlign",
          function(x) .addInserts(subXString(x@string2, start(x@match2), end(x@match2)),
                                  x@inserts2))

setGeneric("type", function(x) standardGeneric("type"))
setMethod("type", "XStringAlign", function(x) x@type)

setGeneric("score", function(x) standardGeneric("score"))
setMethod("score", "XStringAlign", function(x) x@score)

setMethod("length", "XStringAlign", function(x) width(x@match1) + sum(width(x@inserts1)))
setMethod("nchar", "XStringAlign", function(x, type="chars", allowNA=FALSE) length(x))

setMethod("alphabet", "XStringAlign", function(x) alphabet(string1(x)))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method.
###

### TODO: Make the "show" method to format the alignment in a SGD fashion
### i.e. split in 60-letter blocks and use the "|" character to highlight
### exact matches.
setMethod("show", "XStringAlign",
    function(object)
    {
        alignStarts <-
          format(paste("[", c(start(object@match1), start(object@match2)), "]", sep = ""),
                 justify = "right")
        cat(switch(type(object), "global" = "Global", "overlap" = "Overlap",
                   "local" = "Local"), "Pairwise Alignment\n")
        cat("align1:", alignStarts[1],
            toSeqSnippet(align1(object), getOption("width") - 8), "\n")
        cat("align2:", alignStarts[2],
            toSeqSnippet(align2(object), getOption("width") - 8), "\n")
        cat("score:", score(object), "\n")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "as.character" method.
###

setMethod("as.character", "XStringAlign",
    function(x)
    {
        c(align1 = as.character(align1(x)), align2 = as.character(align2(x)))
    }
)

