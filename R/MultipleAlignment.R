### =========================================================================
### MultipleAlignment objects
### -------------------------------------------------------------------------
###
### Everything is now in the MultipleAlignment package!

.call_fun_in_MultipleAlignment <- function(fun, ...)
{
     S4Vectors:::load_package_gracefully("MultipleAlignment",
              "starting with BioC 3.24, calling ", fun, "()")
     msg <- c(fun, "() has moved from Biostrings to the MultipleAlignment ",
              "package, and is formally deprecated in Biostrings >= 2.81.3. ",
              "Please call MultipleAlignment::", fun, "() to get rid of ",
              "this warning.")
     .Deprecated(msg=wmsg(msg))
     FUN <- base::get(fun, envir=asNamespace("MultipleAlignment"),
                      inherits=FALSE)
     do.call(FUN, list(...))
}

rowmask <- function(x)
    .call_fun_in_MultipleAlignment("rowmask", x)

`rowmask<-` <- function(x, append="union", invert=FALSE, value)
    .call_fun_in_MultipleAlignment("rowmask<-",
                                   x, append="union", invert=FALSE, value)

colmask <- function(x)
    .call_fun_in_MultipleAlignment("colmask", x)

`colmask<-` <- function(x, append="union", invert=FALSE, value)
    .call_fun_in_MultipleAlignment("colmask<-",
                                   x, append="union", invert=FALSE, value)

maskGaps <- function(x, ...)
    .call_fun_in_MultipleAlignment("maskGaps", x, ...)

maskednrow <- function(x)
    .call_fun_in_MultipleAlignment("maskednrow", x)

maskedncol <- function(x)
    .call_fun_in_MultipleAlignment("maskedncol", x)

maskeddim <- function(x)
    .call_fun_in_MultipleAlignment("maskeddim", x)

DNAMultipleAlignment <- function(...)
    .call_fun_in_MultipleAlignment("DNAMultipleAlignment", ...)

RNAMultipleAlignment <- function(...)
    .call_fun_in_MultipleAlignment("RNAMultipleAlignment", ...)

AAMultipleAlignment <- function(...)
    .call_fun_in_MultipleAlignment("AAMultipleAlignment", ...)

readDNAMultipleAlignment <- function(...)
    .call_fun_in_MultipleAlignment("readDNAMultipleAlignment", ...)

readRNAMultipleAlignment <- function(...)
    .call_fun_in_MultipleAlignment("readRNAMultipleAlignment", ...)

readAAMultipleAlignment <- function(...)
    .call_fun_in_MultipleAlignment("readAAMultipleAlignment", ...)

write.phylip <- function(...)
    .call_fun_in_MultipleAlignment("write.phylip", ...)

consensusViews <- function(x, ...)
    .call_fun_in_MultipleAlignment("consensusViews", x, ...)

