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
setMethod("maxQuality", "PhredQuality", function(x) 99L)
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
setAs("PhredQuality", "integer",
      function(from) .XStringQualityToInteger(from, "PhredQuality"))
setMethod("as.integer", "PhredQuality", function(x) as(x, "integer"))
setAs("PhredQuality", "numeric",
      function(from) .XStringQualityToNumeric(from, "PhredQuality"))
setMethod("as.numeric", "PhredQuality", function(x) as(x, "numeric"))

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
setAs("SolexaQuality", "integer",
      function(from) .XStringQualityToInteger(from, "SolexaQuality"))
setMethod("as.integer", "SolexaQuality", function(x) as(x, "integer"))
setAs("SolexaQuality", "numeric",
      function(from) .XStringQualityToNumeric(from, "SolexaQuality"))
setMethod("as.numeric", "SolexaQuality", function(x) as(x, "numeric"))

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
setAs("IlluminaQuality", "integer",
      function(from) .XStringQualityToInteger(from, "IlluminaQuality"))
setMethod("as.integer", "IlluminaQuality", function(x) as(x, "integer"))
setAs("IlluminaQuality", "numeric",
      function(from) .XStringQualityToNumeric(from, "IlluminaQuality"))
setMethod("as.numeric", "IlluminaQuality", function(x) as(x, "numeric"))

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

setMethod("encoding", "XStringQuality",
    function(x) 
{
    alf <- alphabet(x)
    setNames(seq(minQuality(x), length.out=length(alf)), alf)
})
