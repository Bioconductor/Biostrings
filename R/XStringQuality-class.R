### =========================================================================
### XStringQuality objects
### -------------------------------------------------------------------------
### An XStringQuality object contains quality information for an XString.

setClass("XStringQuality",
    contains="BStringSet",
    representation("VIRTUAL")
)

setClass("PhredQuality", contains="XStringQuality")

setClass("SolexaQuality", contains="XStringQuality")

setClass("IlluminaQuality", contains="XStringQuality")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Create a methodology for managing quality conversions.
###

setGeneric("offset", function(x) standardGeneric("offset"))
setMethod("offset", "PhredQuality", function(x) 33L)
setMethod("offset", "SolexaQuality", function(x) 64L)
setMethod("offset", "IlluminaQuality", function(x) 64L)

setGeneric("minQuality", function(x) standardGeneric("minQuality"))
setMethod("minQuality", "PhredQuality", function(x) 0L)
setMethod("minQuality", "SolexaQuality", function(x) -5L)
setMethod("minQuality", "IlluminaQuality", function(x) 0L)

setGeneric("maxQuality", function(x) standardGeneric("maxQuality"))
setMethod("maxQuality", "PhredQuality", function(x) 93L) # 126 - 33; valid printable ASCII
setMethod("maxQuality", "SolexaQuality", function(x) 99L)
setMethod("maxQuality", "IlluminaQuality", function(x) 99L)

setGeneric("q2p", function(x) standardGeneric("q2p"))
setMethod("q2p", "PhredQuality", function(x) function(q) 10^(-q/10))
setMethod("q2p", "SolexaQuality",
          function(x) function(q) 1 - 1/(1 + 10^(-q/10)))
setMethod("q2p", "IlluminaQuality", function(x) function(q) 10^(-q/10))

setGeneric("p2q", function(x) standardGeneric("p2q"))
setMethod("p2q", "PhredQuality", function(x) function(p) -10 * log10(p))
setMethod("p2q", "SolexaQuality",
          function(x) function(p) -10 * (log10(p) - log10(1 - p)))
setMethod("p2q", "IlluminaQuality", function(x) function(p) -10 * log10(p))

qualityConverter <- function(x, qualityClass, outputType) {
    .BStringSet2integer <- function(x, scale) {
        if (length(x) == 0)
            value <- integer(0)
        else
            value <- as.integer(charToRaw(as.character(x))) - offset(scale)
        value
    }
    .integer2BStringSet <- function(x, scale) {
        if (length(x) == 0)
            value <- BStringSet()
        else
            value <-
              BStringSet(rawToChar(as.raw(pmax.int(minQuality(scale),
                                          pmin.int(maxQuality(scale), x)) +
                                          offset(scale))))
        value
    }
    .numeric2BStringSet <- function(x, scale) {
        if (length(x) == 0) {
            value <- BStringSet()
        } else if (any(is.na(x)) || any(x < 0) || any(x > 1)) {
            stop("'x' must be numbers between 0 and 1 inclusive")
        } else {
            x <- pmax.int(x, 1e-16)
            x <- pmin.int(x, 1 - 1e-16)
            value <- .integer2BStringSet(as.integer(round(p2q(scale)(x))), scale)
        }
        value
    }
    scale <- new(qualityClass)
    outputType <- match.arg(outputType, c("BStringSet", "integer", "numeric"))
    transform <- paste(class(x), "2", outputType, sep = "")
    switch(transform,
           "BStringSet2BStringSet" =, "character2BStringSet" =,
           "integer2integer" =, "numeric2numeric" = x,
           "BStringSet2integer" =, "character2integer" =
           .BStringSet2integer(x, scale),
           "BStringSet2numeric" =, "character2numeric" =
           q2p(scale)(.BStringSet2integer(x, scale)),
           "integer2BStringSet" = .integer2BStringSet(x, scale),
           "numeric2BStringSet" = .numeric2BStringSet(x, scale),
           "integer2numeric" = q2p(scale)(x),
           "numeric2integer" = as.integer(p2q(scale)(x)))
}



### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

.characterToXStringQuality <- function(from, qualityClass)
    as(BStringSet(from), qualityClass)
.BStringToXStringQuality <- function(from, qualityClass)
    as(BStringSet(from), qualityClass)
.BStringSetToXStringQuality <-
function(from, qualityClass)
{
    ans <- new2(qualityClass, pool=from@pool, ranges=from@ranges, check=FALSE)
    names(ans) <- names(from)
    ans
}
.integerToXStringQuality <- function(from, qualityClass)
    as(qualityConverter(from, qualityClass, "BStringSet"), qualityClass)
.numericToXStringQuality <- function(from, qualityClass)
    as(qualityConverter(from, qualityClass, "BStringSet"), qualityClass)

.IntegerOrNumericListToXStringQuality <- function(from, qualityClass)
    as(relist(as(unlist(from, use.names=FALSE), qualityClass)[[1L]], from),
       qualityClass)

.XStringQualityToInteger <- function(from, qualityClass)
	qualityConverter(BStringSet(from), qualityClass, "integer")
.XStringQualityToNumeric <- function(from, qualityClass)
	qualityConverter(BStringSet(from), qualityClass, "numeric")

.XStringQualityToIntegerMatrix <- function(x)
{
    if (!isConstant(width(x)))
        stop("'x' must be rectangular (i.e. have a constant width)")
    ans <- matrix(as.integer(unlist(x)) - offset(x),
                  nrow=length(x), byrow=TRUE)
    rownames(ans) <- names(x)
    ans
}

setAs("character", "PhredQuality",
      function(from) .characterToXStringQuality(from, "PhredQuality"))
setAs("BString", "PhredQuality",
      function(from) .BStringToXStringQuality(from, "PhredQuality"))
setAs("BStringSet", "PhredQuality",
      function(from) .BStringSetToXStringQuality(from, "PhredQuality"))
setAs("integer", "PhredQuality",
      function(from) .integerToXStringQuality(from, "PhredQuality"))
setAs("numeric", "PhredQuality",
      function(from) .numericToXStringQuality(from, "PhredQuality"))
setAs("IntegerList", "PhredQuality",
    function(from) .IntegerOrNumericListToXStringQuality(from, "PhredQuality"))
setAs("NumericList", "PhredQuality",
    function(from) .IntegerOrNumericListToXStringQuality(from, "PhredQuality"))

setAs("character", "SolexaQuality",
      function(from) .characterToXStringQuality(from, "SolexaQuality"))
setAs("BString", "SolexaQuality",
      function(from) .BStringToXStringQuality(from, "SolexaQuality"))
setAs("BStringSet", "SolexaQuality",
      function(from) .BStringSetToXStringQuality(from, "SolexaQuality"))
setAs("integer", "SolexaQuality",
      function(from) .integerToXStringQuality(from, "SolexaQuality"))
setAs("numeric", "SolexaQuality",
      function(from) .numericToXStringQuality(from, "SolexaQuality"))
setAs("IntegerList", "SolexaQuality",
    function(from) .IntegerOrNumericListToXStringQuality(from, "SolexaQuality"))
setAs("NumericList", "SolexaQuality",
    function(from) .IntegerOrNumericListToXStringQuality(from, "SolexaQuality"))

setAs("character", "IlluminaQuality",
      function(from) .characterToXStringQuality(from, "IlluminaQuality"))
setAs("BString", "IlluminaQuality",
      function(from) .BStringToXStringQuality(from, "IlluminaQuality"))
setAs("BStringSet", "IlluminaQuality",
      function(from) .BStringSetToXStringQuality(from, "IlluminaQuality"))
setAs("integer", "IlluminaQuality",
      function(from) .integerToXStringQuality(from, "IlluminaQuality"))
setAs("numeric", "IlluminaQuality",
      function(from) .numericToXStringQuality(from, "IlluminaQuality"))
setAs("IntegerList", "IlluminaQuality",
    function(from) .IntegerOrNumericListToXStringQuality(from, "IlluminaQuality"))
setAs("NumericList", "IlluminaQuality",
    function(from) .IntegerOrNumericListToXStringQuality(from, "IlluminaQuality"))

setMethod("as.vector", "XStringQuality",
    function(x, mode="any")
    {
        if (!isSingleString(mode)) 
            stop("'mode' must be a single string")
        if (mode %in% "integer")  # return the quality scores
            return(.XStringQualityToInteger(x, class(x)))
        if (mode %in% "numeric")  # return the probabilities
            return(.XStringQualityToNumeric(x, class(x)))
        if (mode %in% c("any", "character"))
            return(as.character(x))
        stop("'mode' can be \"integer\", \"numeric\", \"character\" ",
             "or \"any\" when 'x' is an XStringQuality object")
    }
)
### Return the quality scores.
setAs("XStringQuality", "IntegerList",
    function(from)
        relist(as.integer(as(unlist(from, use.names=FALSE), class(from))),
               from)
)
### Return the probabilities.
setAs("XStringQuality", "NumericList",
    function(from)
        relist(as.numeric(as(unlist(from, use.names=FALSE), class(from))),
               from)
)

### Return the quality scores.
setMethod("as.matrix", "XStringQuality",
    function(x, ...) .XStringQualityToIntegerMatrix(x)
)
setAs("XStringQuality", "matrix", function(from) as.matrix(from))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The user-friendly versatile constructors.

PhredQuality <- function(x) as(x, "PhredQuality")
SolexaQuality <- function(x) as(x, "SolexaQuality")
IlluminaQuality <- function(x) as(x, "IlluminaQuality")

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### alphabet and encoding
###

setMethod("alphabet", "XStringQuality",
    function(x) 
{
    alf <- strsplit(rawToChar(as.raw(33:126)), "")[[1]]
    len <- maxQuality(x) - minQuality(x) + 1L
    alf[seq(offset(x) + minQuality(x) - 32L, length.out=len)]
})

setGeneric("encoding", function(x) standardGeneric("encoding"))

setMethod("encoding", "XStringQuality",
    function(x) 
{
    alf <- alphabet(x)
    setNames(seq(minQuality(x), length.out=length(alf)), alf)
})

