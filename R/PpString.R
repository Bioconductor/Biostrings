### =========================================================================
### Preprocessed strings
### -------------------------------------------------------------------------

### The most general type of preprocessed string (not much we can say about
### it)
setClass("PpString", representation("VIRTUAL"))

### Now a real one!
### Note that 'base_code1', 'base_code2' and 'base_code3' must be distinct
setClass("3BaseCountPpString",
    contains="PpString",
    representation(
        subject="DNAString",        # TODO: support "RNAString" too
        pattern_nchar="integer",    # A single integer e.g. 36
        base_code1="integer",       # integer belonging to DNA_BASE_CODES
        base_count1="XRaw",         # length(base_count1) must be nchar(subject) - pattern_nchar + 1
        base_code2="integer",
        base_count2="XRaw",
        base_code3="integer",
        base_count3="XRaw"
    )
)

### subject <- DNAString("acag-tcgatgc-NNNDR")
### pp <- new("3BaseCountPpString", subject, 4, c("A", "C", "G"))
setMethod("initialize", "3BaseCountPpString",
    function(.Object, subject, pattern_nchar, base_letters)
    {
        .Object@subject <- subject
        if (!is.numeric(pattern_nchar) || length(pattern_nchar) != 1 || is.na(pattern_nchar))
            stop("'pattern_nchar' must be a single integer")
        pattern_nchar <- as.integer(pattern_nchar)
        if (pattern_nchar <= 0)
            stop("'pattern_nchar' must be positive")
        if (pattern_nchar > nchar(subject))
            stop("'pattern_nchar' must be <= 'nchar(subject)'")
        .Object@pattern_nchar <- pattern_nchar
        base_count_length <- nchar(subject) - pattern_nchar + 1
        if (!is.character(base_letters) || length(base_letters) != 3
         || !all(base_letters %in% names(DNA_BASE_CODES)) || any(duplicated(base_letters)))
            stop("'base_letters' must contain 3 distinct DNA base-letters")
        .Object@base_code1 <- DNA_BASE_CODES[base_letters[1]]
        .Object@base_count1 <- XRaw(base_count_length)
        .Object@base_code2 <- DNA_BASE_CODES[base_letters[2]]
        .Object@base_count2 <- XRaw(base_count_length)
        .Object@base_code3 <- DNA_BASE_CODES[base_letters[3]]
        .Object@base_count3 <- XRaw(base_count_length)
        .Object
    }
)

