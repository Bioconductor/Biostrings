# ===========================================================================
# Constructor-like functions and generics for BStringViews objects
# ---------------------------------------------------------------------------

# WARNING: This function is unsafe! (it doesn't check its arguments)
# Only 2 valid ways to use it:
#   new("BStringViews", subject)
#   new("BStringViews", subject, first, last)
# where 'subject' is a BString (or derived) object,
# and 'first' and 'last' are integer vectors of the same length
# such that 'first <= last'.
setMethod("initialize", "BStringViews",
    function(.Object, subject, first, last)
    {
        .Object@subject <- subject
        if (!missing(first)) {
            .Object@first <- first
            .Object@last <- last
        }
        .Object
    }
)

# The 2 functions above share the following properties:
#   - They are exported (and safe).
#   - First argument is 'subject'. It must be a character vector or a BString
#     (or derived) object.
#   - Passing something else to 'subject' provokes an error.
#   - They return a BStringViews object whose 'subject' slot is the object
#     passed in the 'subject' argument.

# Typical use:
#   dna <- DNAString("AA-CC-GG-TT")
# Just one view:
#   dnav1 <- views(dna, 2, 7)
# 9 views, 3 are out of limits:
#   dnav2 <- views(dna, 6:-2, 6:14)
# 5 out of limits views, all have a width of 6:
#   dnav3 <- views(dna, -5:-1, 0:4)
# Same as doing views(dna, 1, length(dna)):
#   dnav4 <- views(dna)
# A BStringViews object with no view:
#   dnav5 <- views(dna, integer(0), integer(0))
views <- function(subject, first=NA, last=NA)
{
    if (class(subject) == "character")
        subject <- BString(subject)
    ans <- new("BStringViews", subject)
    # Integrity checking
    if (!isLooseNumeric(first) || !isLooseNumeric(last))
        stop("'first' and 'last' must be numerics")
    #if (length(first) != length(last))
    #    stop("'first' and 'last' must have the same length")
    if (!is.integer(first))
        first <- as.integer(first)
    first[is.na(first)] <- as.integer(1)
    if (!is.integer(last))
        last <- as.integer(last)
    last[is.na(last)] <- subject@length
    if (length(first) < length(last))
        first <- recycleVector(first, length(last))
    else if (length(last) < length(first))
        last <- recycleVector(last, length(first))
    # The NA-proof version of 'if (any(last < first))'
    if (!isTRUE(all(first <= last)))
        stop("'first' and 'last' must verify 'first <= last'")
    ans@first <- first
    ans@last <- last
    ans
}

# 'width' is the vector of view widths.
# 'gapwidth' is the vector of inter-view widths (recycled).
adjacentViews <- function(subject, width, gapwidth=0)
{
    if (class(subject) == "character")
        subject <- BString(subject)
    ans <- new("BStringViews", subject)
    if (!is.numeric(width) || !isTRUE(all(width >= 1))) # NA-proof
        stop("'width' must be numerics >= 1")
    if (!is.numeric(gapwidth) || !isTRUE(all(gapwidth >= 0))) # NA-proof
        stop("'gapwidth' must be numerics >= 0")
    lw <- length(width)
    if (lw == 0)
        return(ans)
    if (!is.integer(width))
        width <- as.integer(width)
    lg <- length(gapwidth)
    if (!is.integer(gapwidth))
        gapwidth <- as.integer(gapwidth)
    first <- integer(lw)
    last <- integer(lw)
    one <- as.integer(1)
    first[1] <- one
    last[1] <- width[1]
    if (lw >= 2) {
        j <- 1
        for (i in 2:lw) {
            first[i] <- last[i-1] + one + gapwidth[j]
            last[i] <- first[i] + width[i] - one
            if (j < lg) j <- j + 1 else j <- 1
        }
    }
    ans@first <- first
    ans@last <- last
    ans
}


# The BStringViews() generic function
# -----------------------------------
# 'src' should typically be a character vector but the function will also
# work for other kind of input like numeric or even logical vectors.
# 'subjectClass' must be "BString" or one of its derivations ("DNAString",
# "RNAString" or "AAString").
#
# Benchmarks:
#   n <- 40000
#   src <- sapply(1:n, function(i) {paste(sample(DNA_ALPHABET, 250, replace=TRUE), collapse="")})
#   v <- BStringViews(src, "DNAString")
# Comparing BStringViews() speed vs "old" vectorized DNAString() speed:
#       n  BStringViews  "old" DNAString
#   -----  ------------  ---------------
#    5000        0.26 s           4.15 s
#   10000        0.51 s          16.29 s
#   20000        0.99 s          64.85 s
#   40000        1.69 s         488.43 s
# The quadratic behaviour of "old" DNAString() was first reported
# by Wolfgang.
BStringViews <- function(src, subjectClass, sep="")
{
    if (!is.character(sep))
        sep <- toString(sep)
    collapsed <- paste(src, collapse=sep)
    subject <- new(subjectClass, collapsed)
    adjacentViews(subject, nchar(src), nchar(sep))
}

setGeneric(
    "BStringViews",
    function(src, subjectClass, sep="") standardGeneric("BStringViews")
)

# Only FASTA files are supported for now.
# Typical use:
#   srcpath <- system.file("Exfiles", "someORF.fsa", package="Biostrings")
#   f <- file(srcpath)
#   v <- BStringViews(f, "DNAString")
#   close(f)
setMethod("BStringViews", "file",
    function(src, subjectClass, sep)
    {
        fasta <- readFASTA(src)
        src <- sapply(fasta, function(rec) rec$seq)
        desc <- sapply(fasta, function(rec) rec$desc)
        ans <- BStringViews(src, subjectClass, sep) # call the default method
        ans@desc <- desc
        ans
    }
)

# Called when 'src' is a BString (or derived) object.
# When not missing, 'subjectClass' must be "BString" or one of its
# derivations ("DNAString", "RNAString" or "AAString").
setMethod("BStringViews", "BString",
    function(src, subjectClass, sep)
    {
        if (!missing(sep)) {
            # The semantic is: views are delimited by the occurences of 'sep'
            # in 'src' (a kind of strsplit() for BString objects).
            # Uncomment when normalize() and ! method are ready (see TODO file):
            #return(!normalize(matchPattern(sep, b, fixed=TRUE)))
            stop("'sep' not yet supported when 'src' is a \"BString\" object")
        }
        if (!missing(subjectClass) && subjectClass != class(src))
            src <- new(subjectClass, src)
        new("BStringViews", src, as.integer(1), src@length)
    }
)

# Called when 'src' is a BStringViews object.
# 'subjectClass' must be "BString" or one of its derivations ("DNAString",
# "RNAString" or "AAString").
# The 'sep' arg is ignored.
setMethod("BStringViews", "BStringViews",
    function(src, subjectClass, sep)
    {
        if (!missing(sep))
            stop("'sep' not supported when 'src' is a \"BStringViews\" object")
        src@subject <- new(subjectClass, src@subject)
        src
    }
)
