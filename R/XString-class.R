### =========================================================================
### XString objects
### -------------------------------------------------------------------------
###
### The XString virtual class is a general container for storing an "external
### string".
###

setClass("XString", contains="XRaw", representation("VIRTUAL"))

### XString subclasses (no additional slots)
setClass("BString", contains="XString")
setClass("DNAString", contains="XString")
setClass("RNAString", contains="XString")
setClass("AAString", contains="XString")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor-like methods.
###

setMethod("nchar", "XString", function(x, type="chars", allowNA=FALSE) length(x))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "seqtype" and "seqtype<-" methods.
###

setMethod("seqtype", "BString", function(x) "B")
setMethod("seqtype", "DNAString", function(x) "DNA")
setMethod("seqtype", "RNAString", function(x) "RNA")
setMethod("seqtype", "AAString", function(x) "AA")

.copySubSharedRaw <- function(x, start=1, nchar=NA, lkup=NULL)
{
    ans <- SharedRaw(nchar)
    SharedVector.copy(ans, start, start + nchar - 1L, src=x, lkup=lkup)
}

### Downgrades 'x' to a B/DNA/RNA/AAString instance!
setReplaceMethod("seqtype", "XString",
    function(x, value)
    {
        from_seqtype <- seqtype(x)
        to_seqtype <- value
        ans_class <- paste(to_seqtype, "String", sep="")
        lkup <- get_seqtype_conversion_lookup(from_seqtype, to_seqtype)
        if (is.null(lkup))
            return(new(ans_class, shared=x@shared, offset=x@offset, length=x@length))
        shared <- .copySubSharedRaw(x@shared, start=x@offset+1L, nchar=x@length, lkup=lkup)
        new(ans_class, shared=shared, length=length(shared))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### XString.readCodes()
###

XString.readCodes <- function(x, i, imax=integer(0))
{
    SharedRaw.readInts(x@shared, x@offset + i, x@offset + imax)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extract_character_from_XString_by_positions() and
### extract_character_from_XString_by_ranges()
###
### Low-level generics called by the as.character(), show(), and letter()
### methods for XString and XStringViews objects. Not intended to be called
### directly by the end user.
### Purpose is to facilitate support for XString derivatives defined in
### other packages. For example, defining the following methods in the
### Modstrings package will make as.character(), show(), and letter()
### work as expected on ModString and ModStringViews objects (granted
### that seqtype() works properly on ModString derivatives via appropriate
### methods):
###
###   setMethod("extract_character_from_XString_by_positions", "ModString",
###       function(x, pos, collapse=FALSE)
###       {
###           ans <- callNextMethod()
###           codec <- modscodec(seqtype(x))
###           .convert_one_byte_codes_to_letters(ans, codec)
###       }
###   )
###   setMethod("extract_character_from_XString_by_ranges", "ModString",
###       function(x, start, width, collapse=FALSE)
###       {
###           ans <- callNextMethod()
###           codec <- modscodec(seqtype(x))
###           .convert_one_byte_codes_to_letters(ans, codec)
###       }
###   )
###

setGeneric("extract_character_from_XString_by_positions", signature="x",
    function(x, pos, collapse=FALSE)
    {
        ## Only light checking of 'pos' (i.e. we don't check that it contains
        ## valid positions on 'x').
        stopifnot(is(x, "XString"), is.integer(pos))
        ans <- standardGeneric("extract_character_from_XString_by_positions")
        stopifnot(is.character(ans))
        ans
    }
)

### Default method.
setMethod("extract_character_from_XString_by_positions", "XString",
    function(x, pos, collapse=FALSE)
    {
        XVector:::extract_character_from_XRaw_by_positions(x, pos,
                                                           collapse=collapse,
                                                           lkup=xs_dec_lkup(x))
    }
)

setGeneric("extract_character_from_XString_by_ranges", signature="x",
    function(x, start, width, collapse=FALSE)
    {
        ## Only light checking of 'start' and 'width' (i.e. we don't check
        ## that they have the same length and define valid ranges on 'x').
        stopifnot(is(x, "XString"), is.integer(start), is.integer(width))
        ans <- standardGeneric("extract_character_from_XString_by_ranges")
        stopifnot(is.character(ans))
        ans
    }
)

### Default method.
setMethod("extract_character_from_XString_by_ranges", "XString",
    function(x, start, width, collapse=FALSE)
    {
        XVector:::extract_character_from_XRaw_by_ranges(x, start, width,
                                                        collapse=collapse,
                                                        lkup=xs_dec_lkup(x))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### make_XString_from_string()
###
### Low-level generic called by XString() constructor. Not intended to be
### called directly by the end user.
### Purpose is to make it easy to extend the XString() constructor to
### support XString derivatives defined in other packages. For example,
### defining the following method in the Modstrings package will make calls
### of the form 'XString("ModDNA", ...)' work (granted that seqtype() works
### properly on ModDNAString objects via appropriate methods):
###
###   setMethod("make_XString_from_string", "ModString",
###       function(x0, string, start, width)
###       {
###           codec <- modscodec(seqtype(x0))
###           string <- .convert_letters_to_one_byte_codes(string, codec)
###           callNextMethod()
###       }
###   )
###

setGeneric("make_XString_from_string", signature="x0",
    function(x0, string, start, width)
    {
        ## Only light checking of 'start' and 'width' (i.e. we don't check
        ## that they define a valid range on 'string').
        stopifnot(is(x0, "XString"),
                  isSingleInteger(start),
                  isSingleInteger(width))
        if (!isSingleString(string))
            stop(wmsg("input must be a single non-NA string"))
        ans <- standardGeneric("make_XString_from_string")
        stopifnot(class(ans) == class(x0))
        ans
    }
)

### Default method.
setMethod("make_XString_from_string", "XString",
    function(x0, string, start, width)
    {
        lkup <- get_seqtype_conversion_lookup("B", seqtype(x0))
        .Call2("new_XString_from_CHARACTER",
               class(x0), string, start, width, lkup,
               PACKAGE="Biostrings")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The XString() constructor. NOT exported.
###
### This constructor and its helper functions use the uSEW (user-specified
### Start/End/Width) interface.
###

setGeneric("XString", signature="x",
    function(seqtype, x, start=NA, end=NA, width=NA)
        standardGeneric("XString")
)

.charToXString <- function(seqtype, string, start, end, width)
{
    if (!isSingleString(string))
        stop(wmsg("input must be a single non-NA string"))
    x0 <- new2(paste0(seqtype, "String"), check=FALSE)
    solved_SEW <- solveUserSEW(width(string),
                               start=start, end=end, width=width)
    make_XString_from_string(x0, string, start(solved_SEW), width(solved_SEW))
}

setMethod("XString", "character",
    function(seqtype, x, start=NA, end=NA, width=NA)
    {
        if (is.null(seqtype))
            seqtype <- "B"
        .charToXString(seqtype, x, start, end, width)
    }
)

setMethod("XString", "factor",
    function(seqtype, x, start=NA, end=NA, width=NA)
    {
        if (is.null(seqtype))
            seqtype <- "B"
        .charToXString(seqtype, as.character(x), start, end, width)
    }
)

setMethod("XString", "XString",
    function(seqtype, x, start=NA, end=NA, width=NA)
    {
        ans <- subseq(x, start=start, end=end, width=width)
        ## `seqtype<-` must be called even when user supplied 'seqtype' is
        ## NULL because we want to enforce downgrade to a B/DNA/RNA/AAString
        ## instance
        if (is.null(seqtype))
            seqtype <- seqtype(x)
        seqtype(ans) <- seqtype
        ans
    }
)

### Just because of the silly "AsIs" objects found in the probe packages
### (e.g. drosophila2probe$sequence)
setMethod("XString", "AsIs",
    function(seqtype, x, start=NA, end=NA, width=NA)
    {
        if (!is.character(x))
            stop("unsupported input type")
        class(x) <- "character" # keeps the names (unlike as.character())
        XString(seqtype, x, start=start, end=end, width=width)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The user interfaces to the XString() constructor.
###

BString <- function(x="", start=1, nchar=NA)
    XString("B", x, start=start, width=nchar)

DNAString <- function(x="", start=1, nchar=NA)
    XString("DNA", x, start=start, width=nchar)

RNAString <- function(x="", start=1, nchar=NA)
    XString("RNA", x, start=start, width=nchar)

AAString <- function(x="", start=1, nchar=NA)
    XString("AA", x, start=start, width=nchar)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

setAs("XString", "BString",
    function(from) {seqtype(from) <- "B"; from}
)
setAs("XString", "DNAString",
    function(from) {seqtype(from) <- "DNA"; from}
)
setAs("XString", "RNAString",
    function(from) {seqtype(from) <- "RNA"; from}
)
setAs("XString", "AAString",
    function(from) {seqtype(from) <- "AA"; from}
)

setAs("character", "BString", function(from) BString(from))
setAs("character", "DNAString", function(from) DNAString(from))
setAs("character", "RNAString", function(from) RNAString(from))
setAs("character", "AAString", function(from) AAString(from))
setAs("character", "XString", function(from) BString(from))

setMethod("as.character", "XString",
    function(x)
        extract_character_from_XString_by_ranges(x, 1L, length(x))
)

setMethod("toString", "XString", function(x, ...) as.character(x))

### FIXME: Sometimes returns a vector sometimes a factor. This needs to be
### sorted out. The use case is that as.data.frame() relies on this.
setMethod("as.vector", "XString",
    function(x)
    {
        codes <- xscodes(x)
        x_alphabet <- names(codes)
        if (is.null(x_alphabet)) {
            ans <- rawToChar(as.raw(x), multiple=TRUE)
            x_alphabet <- alphabet(x)
            if (!is.null(x_alphabet))
                ans <- factor(ans, levels=x_alphabet)
            return(ans)
        }
        code2pos <- integer(length(codes))
        code2pos[codes] <- seq_along(codes)
        ans <- code2pos[as.integer(x)]
        attributes(ans) <- list(levels=x_alphabet, class="factor")
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method.
###

compact_ellipsis <- rawToChar(as.raw(c(0xe2, 0x80, 0xa6)))

### NOT exported but used in the BSgenome package.
### 'x' must be a single character string, or an XString or
### MaskedXString object.
### Return a character vector possibly with a class attribute on it for
### later S3 dispatch in add_colors().
toSeqSnippet <- function(x, width)
{
    if (width < 7L)
        width <- 7L
    ## Do NOT use nchar() here as it wouldn't do the right thing on a
    ## MaskedXString object!
    x_len <- length(x)
    if (x_len <= width) {
        ans <- as.character(x)
    } else {
        w1 <- (width - 2L) %/% 2L
        w2 <- (width - 3L) %/% 2L
        ans <- paste0(as.character(subseq(x, start=1, width=w1)),
                      #compact_ellipsis,
                      "...",
                      as.character(subseq(x, end=x_len, width=w2)))
    }
    if (is(x, "XString") || is(x, "MaskedXString"))
        class(ans) <- c(seqtype(x), class(ans))  # for S3 dispatch
                                                 # in add_colors()
    ans
}

setMethod("show", "XString",
    function(object)
    {
        object_len <- object@length
        cat(object_len, "-letter ", class(object), " object\n", sep="")
        snippet <- toSeqSnippet(object, getOption("width") - 5L)
        cat("seq: ", add_colors(snippet), "\n", sep="")
    }
)

setMethod("showAsCell", "XString",
    function(object)
    {
        ans <- safeExplode(as.character(object))
        class(ans) <- c(seqtype(object), class(ans))
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Equality.
###
### We want:
###   BString("ab") == "ab" # TRUE
###   DNAString("TG") == RNAString("UG") # TRUE!!!
###   library(BSgenome.Hsapiens.UCSC.hg18)
###   dna <- Hsapiens$chr1
###   dna != Hsapiens$chr1 # FALSE
###   dnav <- Views(dna, start=1:7, end=101:107)
###   dnav[[1]] == dnav[[7]] # TRUE
###   dnav <- Views(dna, start=1:7, end=(length(dna)-6):length(dna))
### This is fast:
###   dnav[[1]] == dnav[[7]] # FALSE
### But this would have killed your machine:
###   s1 <- toString(dnav[[1]])
###   s7 <- toString(dnav[[7]])
###   s1 == s7

### 'x' and 'y' must be XString objects
.XString.equal <- function(x, y)
{
    if (x@length != y@length)
        return(FALSE)
    ans <- !SharedVector.compare(x@shared, x@offset + 1L, y@shared, y@offset + 1L, x@length)
    as.logical(ans)
}

setMethod("==", signature(e1="XString", e2="XString"),
    function(e1, e2)
    {
        if (!comparable_seqtypes(seqtype(e1), seqtype(e2))) {
            class1 <- class(e1)
            class2 <- class(e2)
            stop("comparison between a \"", class1, "\" instance ",
                 "and a \"", class2, "\" instance ",
                 "is not supported")
        }
        .XString.equal(e1, e2)
    }
)
setMethod("==", signature(e1="BString", e2="character"),
    function(e1, e2)
    {
        if (length(e2) != 1 || e2 %in% c("", NA))
            stop("comparison between a \"BString\" object and a character vector ",
                 "of length != 1 or an empty string or an NA ",
                 "is not supported")
        .XString.equal(e1, BString(e2))
    }
)
setMethod("==", signature(e1="character", e2="BString"),
    function(e1, e2) e2 == e1
)

setMethod("!=", signature(e1="XString", e2="XString"),
    function(e1, e2) !(e1 == e2)
)
setMethod("!=", signature(e1="BString", e2="character"),
    function(e1, e2) !(e1 == e2)
)
setMethod("!=", signature(e1="character", e2="BString"),
    function(e1, e2) !(e1 == e2)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "substr" and "substring" methods.
###

setMethod("substr", "XString",
    function(x, start, stop) subseq(x, start=start, end=stop)
)

setMethod("substring", "XString",
    function(text, first, last=1000000L) subseq(text, start=first, end=last)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### updateObject()
###

### Update XString objects created before the big internal renaming I made
### in IRanges 1.3.76.
setMethod("updateObject", "XString",
    function(object, ..., verbose=FALSE)
    {
        if (!is(try(object@shared, silent=TRUE), "try-error"))
            return(object)
        xdata <- object@xdata
        ans_shared <- new("SharedRaw")
        ans_shared@xp <- xdata@xp
        ans_shared@.link_to_cached_object=xdata@.link_to_cached_object
        new2(class(object),
             shared=ans_shared,
             offset=object@offset,
             length=object@length,
             check=FALSE)
    }
)

### Update AAString objects created before AA_ALPHABET was enforced
setMethod("updateObject", "AAString",
    function(object, ..., verbose=FALSE)
    {
        ## Start by calling the updateObject() method for XString objects.
        object <- callNextMethod()

        codec <- xscodec(AAString())
        class(object) <- "BString"
        mapping <- vapply(uniqueLetters(object), utf8ToInt, integer(1L))
        missingVals <- is.na(codec@enc_lkup[mapping+1L])
        if(any(missingVals)){
            errorChars <- paste(names(mapping)[which(missingVals)],
                                collapse=', ')
            stop("Cannot decode, AAString contains invalid character(s): ",
                  errorChars)
        }
        AAString(object)
    }
)
