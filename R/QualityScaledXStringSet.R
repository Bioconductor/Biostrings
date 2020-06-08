### =========================================================================
### QualityScaledXStringSet objects
### -------------------------------------------------------------------------
###


setClass("QualityScaledXStringSet",
    contains="XStringSet",
    representation(
        "VIRTUAL",
        quality="XStringQuality"
    )
)

### Combine the new "parallel slots" with those of the parent class. Make
### sure to put the new parallel slots **first**. See R/Vector-class.R file
### in the S4Vectors package for what slots should or should not be considered
### "parallel".
setMethod("parallel_slot_names", "QualityScaledXStringSet",
    function(x) c("quality", callNextMethod())
)

### QualityScaledXStringSet subclasses
setClass("QualityScaledBStringSet",
    contains=c("BStringSet", "QualityScaledXStringSet")
)
setClass("QualityScaledDNAStringSet",
    contains=c("DNAStringSet", "QualityScaledXStringSet")
)
setClass("QualityScaledRNAStringSet",
    contains=c("RNAStringSet", "QualityScaledXStringSet")
)
setClass("QualityScaledAAStringSet",
    contains=c("AAStringSet", "QualityScaledXStringSet")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

.valid.QualityScaledXStringSet <- function(object)
{
    message <- NULL
    if (!all(nchar(object@quality) == 1 | nchar(object@quality) == nchar(object)))
        message <- c(message, "'nchar(quality)' must equal 1 or nchar of 'XStringSet'")
    message
}

setValidity("QualityScaledXStringSet",
    function(object)
    {
        problems <- .valid.QualityScaledXStringSet(object)
        if (is.null(problems)) TRUE else problems
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods
###

setGeneric("quality", function(x) standardGeneric("quality"), useAsDefault = function(x) x@quality)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The user-friendly versatile constructors
###

### The returned 'quality' is guaranteed to have the shape of 'x' (i.e. same
### length() and width()).
.normarg_quality <- function(quality, x)
{
    if (!is(quality, "XStringQuality"))
        stop("'quality' must be of class 'XStringQuality'")
    quality_width <- width(quality)
    x_width <- width(x)
    if (length(quality) == length(x)) {
        recycle_me <- quality_width != x_width
        if (any(recycle_me & quality_width != 1L))
            stop(wmsg("the quality strings must be of length 1 or have the ",
                      "same length as their corresponding string in 'x'"))
        recycle_idx <- which(recycle_me)
        width2 <- x_width[recycle_idx]
        idx <- relist(rep.int(1L, sum(width2)), PartitioningByWidth(width2))
        quality[recycle_idx] <- quality[recycle_idx][idx]
        return(quality)
    }
    if (length(quality) == 1L) {
        if (all(x_width == quality_width))
            return(rep.int(quality, length(x)))
        if (quality_width != 1L)
            stop(wmsg("when 'quality' is a single string it must be ",
                      "a single letter or have the same width as all ",
                      "the strings in 'x'"))
        quality <- PhredQuality(BStringSet(rep.int(quality[[1L]], max(x_width)),
                                           start=1L, end=x_width))
        return(quality)
    }
    stop("'length(quality)' must equal 'length(x)' or 1")
}

QualityScaledXStringSet <- function(x, quality) {
    if (!is(x, "XStringSet"))
        stop("'x' must be of class 'XStringSet'")
    quality <- .normarg_quality(quality, x)
    output <- as(x, paste0("QualityScaled", class(x)))
    slot(output, "quality", check=FALSE) <- quality
    output
}

QualityScaledBStringSet <- function(x, quality) QualityScaledXStringSet(BStringSet(x), quality)
QualityScaledDNAStringSet <- function(x, quality) QualityScaledXStringSet(DNAStringSet(x), quality)
QualityScaledRNAStringSet <- function(x, quality) QualityScaledXStringSet(RNAStringSet(x), quality)
QualityScaledAAStringSet <- function(x, quality) QualityScaledXStringSet(AAStringSet(x), quality)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Overwrite some endomorphic methods for XStringSet objects
###

### The default "narrow" method calls windows() so we only need to implement
### a "windows" method for QualityScaledXStringSet objects to make narrow()
### work on these objects.
setMethod("windows", "QualityScaledXStringSet",
    function(x, start=NA, end=NA, width=NA)
    {
        x@quality <- windows(x@quality, start=start, end=end, width=width)
        callNextMethod()
    }
)

setMethod("reverse", "QualityScaledXStringSet",
    function(x)
    {
        x@quality <- reverse(x@quality)
        callNextMethod()
    }
)

setMethod("reverseComplement", "QualityScaledDNAStringSet",
    function(x)
    {
        x@quality <- reverse(x@quality)
        callNextMethod()
    }
)

setMethod("reverseComplement", "QualityScaledRNAStringSet",
    function(x)
    {
        x@quality <- reverse(x@quality)
        callNextMethod()
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show()
###

setMethod("show", "QualityScaledXStringSet",
    function(object)
    {
        cat("  A ", class(object), " instance containing:\n", sep="")
        cat("\n")
        selectMethod("show", "XStringSet")(as(object, "XStringSet"))
        cat("\n")
        show(quality(object))
        cat("\n")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### readQualityScaledDNAStringSet() / writeQualityScaledXStringSet()
###

readQualityScaledDNAStringSet <- function(filepath,
                       quality.scoring=c("phred", "solexa", "illumina"),
                       nrec=-1L, skip=0L, seek.first.rec=FALSE,
                       use.names=TRUE)
{
    quality.scoring <- match.arg(quality.scoring)
    x <- readDNAStringSet(filepath, format="fastq",
                          nrec, skip, seek.first.rec,
                          use.names, with.qualities=TRUE)
    qualities <- mcols(x)[ , "qualities"]
    quals <- switch(quality.scoring,
                    phred=PhredQuality(qualities),
                    solexa=SolexaQuality(qualities),
                    illumina=IlluminaQuality(qualities))
    QualityScaledDNAStringSet(x, quals)
}

writeQualityScaledXStringSet <- function(x, filepath,
                       append=FALSE, compress=FALSE, compression_level=NA)
{
    if (!is(x, "QualityScaledXStringSet"))
        stop(wmsg("'x' must be a QualityScaledXStringSet object"))
    writeXStringSet(x, filepath, append, compress, compression_level,
                       format="fastq", qualities=quality(x))
}

