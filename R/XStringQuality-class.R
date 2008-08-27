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
### Coercion.
###

setAs("character", "PhredQuality",
      function(from) as(BStringSet(from), "PhredQuality")
)

setAs("BString", "PhredQuality",
      function(from) as(BStringSet(from), "PhredQuality")
)

setAs("BStringSet", "PhredQuality",
      function(from)
      new("PhredQuality", super(from), start(from), width(from), names(from))
)

setAs("integer", "PhredQuality",
      function(from) {
          if (any(is.na(from)) || any(from < 0L) || any(from > 99L))
              stop("'from' must be integers between 0 and 99 inclusive")
          new("PhredQuality",
              super = BString(rawToChar(as.raw(33L + from))),
              start = 1,
              width = length(from),
              names = NULL)
      }
)

setAs("numeric", "PhredQuality",
      function(from) {
          if (any(is.na(from)) || any(from < 0) || any(from > 1))
              stop("'from' must be numbers between 0 and 1 inclusive")
          from[from < 1e-12] <- 1e-12
          as(pmin(99L, as.integer(round(-10 * log10(from)))), "PhredQuality")
      }
)

setAs("PhredQuality", "integer",
      function(from) {
          as.integer(charToRaw(as.character(from))) - 33L
      }
)
setMethod("as.integer", "PhredQuality", function(x) as(x, "integer"))

setAs("PhredQuality", "numeric",
      function(from) {
          Q <- as.integer(charToRaw(as.character(from))) - 33L
          10^(-Q/10)
      }
)
setMethod("as.numeric", "PhredQuality", function(x) as(x, "numeric"))

setAs("character", "SolexaQuality",
      function(from) as(BStringSet(from), "SolexaQuality")
)

setAs("BString", "SolexaQuality",
      function(from) as(BStringSet(from), "SolexaQuality")
)

setAs("BStringSet", "SolexaQuality",
      function(from)
      new("SolexaQuality", super(from), start(from), width(from), names(from))
)

setAs("integer", "SolexaQuality",
      function(from) {
          if (any(is.na(from)) || any(from < -5L) || any(from > 99L))
              stop("'from' must be integers between -5 and 99 inclusive")
          new("SolexaQuality",
              super = BString(rawToChar(as.raw(64L + pmax(-5L, from)))),
              start = 1,
              width = length(from),
              names = NULL)
      }
)

setAs("numeric", "SolexaQuality",
      function(from) {
          if (any(is.na(from)) || any(from < 0) || any(from > 1))
              stop("'from' must be numbers between 0 and 1 inclusive")
          from[from < 1e-12] <- 1e-12
          as(pmax(-5L, pmin(99L, as.integer(round(-10 * log10(from/(1 - from)))))), "SolexaQuality")
      }
)

setAs("SolexaQuality", "integer",
      function(from) {
          as.integer(charToRaw(as.character(from))) - 64L
      }
)
setMethod("as.integer", "SolexaQuality", function(x) as(x, "integer"))

setAs("SolexaQuality", "numeric",
      function(from) {
          Q <- as.integer(charToRaw(as.character(from))) - 64L
          1 - 1/(1 + 10^(-Q/10))
      }
)
setMethod("as.numeric", "SolexaQuality", function(x) as(x, "numeric"))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The user-friendly versatile constructors.

PhredQuality <- function(x) as(x, "PhredQuality")
SolexaQuality <- function(x) as(x, "SolexaQuality")
