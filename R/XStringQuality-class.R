### =========================================================================
### XStringQuality objects
### -------------------------------------------------------------------------
### An XStringQuality object contains quality information for an XString.


setClass("XStringQuality", representation("VIRTUAL"))

setClass("PhredQuality",
    contains="XStringSet",
    representation("XStringQuality", super = "BString")
)

setClass("SolexaQuality",
    contains="XStringSet",
    representation("XStringQuality", super = "BString")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Create a methodology for managing quality converstions.
###

setGeneric("offset", function(x) standardGeneric("offset"))
setMethod("offset", "PhredQuality", function(x) 33L)
setMethod("offset", "SolexaQuality", function(x) 64L)

setGeneric("minQuality", function(x) standardGeneric("minQuality"))
setMethod("minQuality", "PhredQuality", function(x) 0L)
setMethod("minQuality", "SolexaQuality", function(x) -5L)

setGeneric("maxQuality", function(x) standardGeneric("maxQuality"))
setMethod("maxQuality", "PhredQuality", function(x) 99L)
setMethod("maxQuality", "SolexaQuality", function(x) 99L)

setGeneric("q2p", function(x) standardGeneric("q2p"))
setMethod("q2p", "PhredQuality", function(x) function(q) 10^(-q/10))
setMethod("q2p", "SolexaQuality", function(x) function(q) 1 - 1/(1 + 10^(-q/10)))

setGeneric("p2q", function(x) standardGeneric("p2q"))
setMethod("p2q", "PhredQuality", function(x) function(p) -10 * log10(p))
setMethod("p2q", "SolexaQuality", function(x) function(p) -10 * log10(p/(1 - p)))

qualityConverter <- function(x, qualityClass, outputType) {
    .BStringSet2integer <- function(x, scale)
        as.integer(charToRaw(as.character(x))) - offset(scale)
    .integer2BStringSet <- function(x, scale)
        BStringSet(rawToChar(as.raw(pmax.int(minQuality(scale), pmin.int(maxQuality(scale), x)) + offset(scale))))
    .numeric2BStringSet <- function(x, scale) {
        if (any(is.na(x)) || any(x < 0) || any(x > 1))
            stop("'x' must be numbers between 0 and 1 inclusive")
        x[x < 1e-12] <- 1e-12
        .integer2BStringSet(as.integer(round(p2q(scale)(x))), scale)
    }
    scale <- new(qualityClass)
    outputType <- match.arg(outputType, c("BStringSet", "integer", "numeric"))
    transform <- paste(class(x), "2", outputType, sep = "")
    switch(transform,
           "BStringSet2BStringSet" =, "character2BStringSet" =, "integer2integer" =, "numeric2numeric" = x,
           "BStringSet2integer" =, "character2integer" = .BStringSet2integer(x, scale),
           "BStringSet2numeric" =, "character2numeric" = q2p(scale)(.BStringSet2integer(x, scale)),
           "integer2BStringSet" = .integer2BStringSet(x, scale),
           "numeric2BStringSet" = .numeric2BStringSet(x, scale),
           "integer2numeric" = q2p(scale)(x),
           "numeric2integer" = as.integer(p2q(scale)(x)))
}



### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

.characterToXStringQuality <- function(from, qualityClass) as(BStringSet(from), qualityClass)
.BStringToXStringQuality <- function(from, qualityClass) as(BStringSet(from), qualityClass)
.BStringSetToXStringQuality <- function(from, qualityClass) {
    ans <- new(qualityClass, super=super(from), start=start(from), width=width(from))
    names(ans) <- names(from)
    ans
}
.integerToXStringQuality <- function(from, qualityClass)
    as(qualityConverter(from, qualityClass, "BStringSet"), qualityClass)
.numericToXStringQuality <- function(from, qualityClass)
    as(qualityConverter(from, qualityClass, "BStringSet"), qualityClass)
.XStringQualityToInteger <- function(from, qualityClass)
	as(qualityConverter(from, qualityClass, "integer"), qualityClass)
.XStringQualityToNumeric <- function(from, qualityClass)
	as(qualityConverter(from, qualityClass, "numeric"), qualityClass)

setAs("character", "PhredQuality", function(from) .characterToXStringQuality(from, "PhredQuality"))
setAs("BString", "PhredQuality", function(from) .BStringToXStringQuality(from, "PhredQuality"))
setAs("BStringSet", "PhredQuality", function(from) .BStringSetToXStringQuality(from, "PhredQuality"))
setAs("integer", "PhredQuality", function(from) .integerToXStringQuality(from, "PhredQuality"))
setAs("numeric", "PhredQuality", function(from) .numericToXStringQuality(from, "PhredQuality"))
setAs("PhredQuality", "integer", function(from) .XStringQualityToInteger(from, "PhredQuality"))
setMethod("as.integer", "PhredQuality", function(x) as(x, "integer"))
setAs("PhredQuality", "numeric", function(from) .XStringQualityToNumeric(from, "PhredQuality"))
setMethod("as.numeric", "PhredQuality", function(x) as(x, "numeric"))

setAs("character", "SolexaQuality", function(from) .characterToXStringQuality(from, "SolexaQuality"))
setAs("BString", "SolexaQuality", function(from) .BStringToXStringQuality(from, "SolexaQuality"))
setAs("BStringSet", "SolexaQuality", function(from) .BStringSetToXStringQuality(from, "SolexaQuality"))
setAs("integer", "SolexaQuality", function(from) .integerToXStringQuality(from, "SolexaQuality"))
setAs("numeric", "SolexaQuality", function(from) .numericToXStringQuality(from, "SolexaQuality"))
setAs("SolexaQuality", "integer", function(from) .XStringQualityToInteger(from, "SolexaQuality"))
setMethod("as.integer", "SolexaQuality", function(x) as(x, "integer"))
setAs("SolexaQuality", "numeric", function(from) .XStringQualityToNumeric(from, "SolexaQuality"))
setMethod("as.numeric", "SolexaQuality", function(x) as(x, "numeric"))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The user-friendly versatile constructors.

PhredQuality <- function(x) as(x, "PhredQuality")
SolexaQuality <- function(x) as(x, "SolexaQuality")
