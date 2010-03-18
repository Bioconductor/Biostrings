### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "trimLRPatterns" generic.
###

setGeneric("trimLRPatterns", signature = "subject",
    function(Lpattern = "", Rpattern = "", subject,
             max.Lmismatch = 0, max.Rmismatch = 0,
             with.Lindels = FALSE, with.Rindels = FALSE,
             Lfixed = TRUE, Rfixed = TRUE, ranges = FALSE)
        standardGeneric("trimLRPatterns")
)

.XString.XStringSet.trimLRPatterns <-
function(Lpattern = "", Rpattern = "", subject,
         max.Lmismatch = 0, max.Rmismatch = 0,
         with.Lindels = FALSE, with.Rindels = FALSE,
         Lfixed = TRUE, Rfixed = TRUE, ranges = FALSE)
{
    if (nchar(Rpattern) == 0) {
        trim.end <- nchar(subject)
    } else if (length(unique(nchar(subject))) != 1) {
        Rpattern <- normargPattern(Rpattern, subject, argname = "Rpattern")
        reversed.ranges <- 
          .XString.XStringSet.trimLRPatterns(Lpattern = reverse(Rpattern),
                                             subject = reverse(subject),
                                             max.Lmismatch = max.Rmismatch,
                                             with.Lindels = with.Rindels,
                                             Lfixed = Rfixed, ranges = TRUE)
        trim.end <- nchar(subject) - start(reversed.ranges) + 1L
    } else {
        ncharRpattern <- nchar(Rpattern)
        if (length(max.Rmismatch) == 1 && max.Rmismatch >= 0 && max.Rmismatch < 1) {
            max.Rmismatch <- max.Rmismatch * seq_len(ncharRpattern)
        }
        max.Rmismatch <- as.integer(max.Rmismatch)
        if (length(max.Rmismatch) < ncharRpattern) {
            max.Rmismatch <-
              c(rep.int(-1L, ncharRpattern - length(max.Rmismatch)), max.Rmismatch)
        }
        if (any(is.na(max.Rmismatch)) || length(max.Rmismatch) != ncharRpattern) {
            stop("'max.Rmismatch' must be a vector of length 'nchar(Rpattern)'")
        }
        with.Rindels <- normargWithIndels(with.Rindels, argname = "with.Rindels")
        Rfixed <- normargFixed(Rfixed, subject, argname = "Rfixed")

        # test the pattern "from the inside out"
        ncharSubject <- nchar(subject)[1]

        # reverse for the inside-out search
        max.Rmismatch <- rev(max.Rmismatch)

        # get start of first match, or NA
        Rwhich.trim <-
          ncharSubject - ncharRpattern +
            which.isMatchingEndingAt(pattern = Rpattern,
                                     subject = subject,
                                     ending.at = ncharSubject,
                                     max.mismatch = max.Rmismatch,
                                     with.indels = with.Rindels,
                                     fixed = Rfixed,
                                     auto.reduce.pattern = TRUE)

        # last position before best match starts
        trim.end <- ifelse(is.na(Rwhich.trim), ncharSubject, Rwhich.trim - 1L)
    }

    if (nchar(Lpattern) == 0) {
        trim.start <- pmin(1L, nchar(subject))
    } else {
        ncharLpattern <- nchar(Lpattern)
        if (length(max.Lmismatch) == 1 && max.Lmismatch >= 0 && max.Lmismatch < 1) {
            max.Lmismatch <- max.Lmismatch * seq_len(ncharLpattern)
        }
        max.Lmismatch <- as.integer(max.Lmismatch)
        if (length(max.Lmismatch) < ncharLpattern) {
            max.Lmismatch <-
              c(rep.int(-1L, ncharLpattern - length(max.Lmismatch)), max.Lmismatch)
        }
        if (any(is.na(max.Lmismatch)) || length(max.Lmismatch) != ncharLpattern) {
            stop("'max.Lmismatch' must be a vector of length 'nchar(Lpattern)'")
        }
        with.Lindels <- normargWithIndels(with.Lindels, argname = "with.Lindels")
        Lfixed <- normargFixed(Lfixed, subject, argname = "Lfixed")

        # test the pattern "from the inside out"
        # reverse for the inside-out search
        max.Lmismatch <- rev(max.Lmismatch)

        # get end of first match, or NA
        Lwhich.trim <-
          1L + ncharLpattern -
            which.isMatchingStartingAt(pattern = Lpattern,
                                       subject = subject,
                                       starting.at = 1,
                                       max.mismatch = max.Lmismatch,
                                       with.indels = with.Lindels,
                                       fixed = Lfixed,
                                       auto.reduce.pattern = TRUE)

        # next position after best match ends
        trim.start <- ifelse(is.na(Lwhich.trim), pmin(1L, nchar(subject)), Lwhich.trim + 1L)
    }
    trim.start <- ifelse(trim.end >= 1L, trim.start, 2L)
    trim.end <- ifelse(trim.end >= trim.start, trim.end, trim.start - 1L)
    if (ranges) {
        IRanges(start = trim.start, end = trim.end)
    } else if (is(subject, "XString")) {
        subseq(subject, start = trim.start, end = trim.end)
    } else {
        narrow(subject, start = trim.start, end = trim.end)
    }
}

### Dispatch on 'subject' (see signature of generic).
setMethod("trimLRPatterns", "XString",
    function(Lpattern = "", Rpattern = "", subject,
             max.Lmismatch = 0, max.Rmismatch = 0,
             with.Lindels = FALSE, with.Rindels = FALSE,
             Lfixed = TRUE, Rfixed = TRUE, ranges = FALSE)
    .XString.XStringSet.trimLRPatterns(Lpattern = Lpattern,
                                       Rpattern = Rpattern,
                                       subject = subject,
                                       max.Lmismatch = max.Lmismatch,
                                       max.Rmismatch = max.Rmismatch,
                                       with.Lindels = with.Lindels,
                                       with.Rindels = with.Rindels,
                                       Lfixed = Lfixed,
                                       Rfixed = Rfixed,
                                       ranges = ranges))

### Dispatch on 'subject' (see signature of generic).
setMethod("trimLRPatterns", "XStringSet",
    function(Lpattern = "", Rpattern = "", subject,
             max.Lmismatch = 0, max.Rmismatch = 0,
             with.Lindels = FALSE, with.Rindels = FALSE,
             Lfixed = TRUE, Rfixed = TRUE, ranges = FALSE)
    .XString.XStringSet.trimLRPatterns(Lpattern = Lpattern,
                                       Rpattern = Rpattern,
                                       subject = subject,
                                       max.Lmismatch = max.Lmismatch,
                                       max.Rmismatch = max.Rmismatch,
                                       with.Lindels = with.Lindels,
                                       with.Rindels = with.Rindels,
                                       Lfixed = Lfixed,
                                       Rfixed = Rfixed,
                                       ranges = ranges))

### Dispatch on 'subject' (see signature of generic).
setMethod("trimLRPatterns", "character",
    function(Lpattern = "", Rpattern = "", subject,
             max.Lmismatch = 0, max.Rmismatch = 0,
             with.Lindels = FALSE, with.Rindels = FALSE,
             Lfixed = TRUE, Rfixed = TRUE, ranges = FALSE) {
        if (length(subject) == 1)
            subject <- BString(subject)
        else
            subject <- BStringSet(subject)
        ans <-
          .XString.XStringSet.trimLRPatterns(Lpattern = Lpattern,
                                             Rpattern = Rpattern,
                                             subject = subject,
                                             max.Lmismatch = max.Lmismatch,
                                             max.Rmismatch = max.Rmismatch,
                                             with.Lindels = with.Lindels,
                                             with.Rindels = with.Rindels,
                                             Lfixed = Lfixed,
                                             Rfixed = Rfixed,
                                             ranges = ranges)
        if (!ranges)
            ans <- as.character(ans)
        ans
    })
