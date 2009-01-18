### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "trimLRPatterns" new generic.
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
    if (is(subject, "XString")) {
        subjectStringClass <- class(subject)
    } else {
        subjectStringClass <- baseXStringSubtype(subject)
    }
    if (nchar(Rpattern) == 0) {
        trim.end <- nchar(subject)
    } else if (length(unique(nchar(subject))) != 1) {
        if (class(Rpattern) != class(super(subject)))
            Rpattern <- XString(subjectStringClass, Rpattern)
        reversed.ranges <- 
          .XString.XStringSet.trimLRPatterns(Lpattern = reverse(Rpattern),
                                             subject = reverse(subject),
                                             max.Lmismatch = max.Rmismatch,
                                             with.Lindels = with.Rindels,
                                             Lfixed = Rfixed, ranges = TRUE)
        trim.end <- nchar(subject) - start(reversed.ranges) + 1L
    } else {
        ncharRpattern <- nchar(Rpattern)
        if (length(max.Rmismatch) == 1 && max.Rmismatch > 0 && max.Rmismatch < 1) {
            max.Rmismatch <- max.Rmismatch * seq_len(ncharRpattern)
        }
        max.Rmismatch <- as.integer(max.Rmismatch)
        if (length(max.Rmismatch) < ncharRpattern) {
            max.Rmismatch <-
              c(rep.int(-1L, ncharRpattern - length(max.Rmismatch)), max.Rmismatch)
        }
        if (any(is.na(max.Rmismatch)) || length(max.Rmismatch) != ncharRpattern)
            stop("'max.Rmismatch' must be a vector of length 'nchar(Rpattern)'")
        with.Rindels <- normargWithIndels(with.Rindels, name = "with.Rindels")
        Rfixed <- normargFixed(Rfixed, subjectStringClass, name = "Rfixed")

        ncharSubject <- nchar(subject)[1]
        Rstarting.at <- ncharSubject:(ncharSubject - ncharRpattern + 1L)
        Rcandidate <-
          t(t(neditStartingAt(Rpattern, subject, starting.at = Rstarting.at,
                              with.indels = with.Rindels, fixed = Rfixed)) <=
            max.Rmismatch)
        if (is(subject, "XString")) {
            Rwhich.candidate <- which.max(Rcandidate)
            Rwhich.trim <- any(Rcandidate)
        } else {
            Rwhich.candidate <- max.col(Rcandidate, ties.method = "first")
            Rwhich.trim <- rowSums(Rcandidate) > 0
        }
        trim.end <-
          ifelse(Rwhich.trim, ncharSubject - Rwhich.candidate, nchar(subject))
    }

    if (nchar(Lpattern) == 0) {
        trim.start <- pmin(1L, nchar(subject))
    } else {
        ncharLpattern <- nchar(Lpattern)
        if (length(max.Lmismatch) == 1 && max.Lmismatch > 0 && max.Lmismatch < 1) {
            max.Lmismatch <- max.Lmismatch * seq_len(ncharLpattern)
        }
        max.Lmismatch <- as.integer(max.Lmismatch)
        if (length(max.Lmismatch) < ncharLpattern) {
            max.Lmismatch <-
              c(rep.int(-1L, ncharLpattern - length(max.Lmismatch)), max.Lmismatch)
        }
        if (any(is.na(max.Lmismatch)) || length(max.Lmismatch) != ncharLpattern)
            stop("'max.Lmismatch' must be a vector of length 'nchar(Lpattern)'")
        with.Lindels <- normargWithIndels(with.Lindels, name = "with.Lindels")
        Lfixed <- normargFixed(Lfixed, subjectStringClass, name = "Lfixed")
        
        Lending.at <- seq_len(ncharLpattern)
        Lcandidate <-
          t(t(neditEndingAt(Lpattern, subject, ending.at = Lending.at,
                            with.indels = with.Lindels, fixed = Lfixed)) <=
            max.Lmismatch)
        if (is(subject, "XString")) {
            Lwhich.candidate <- which.max(Lcandidate)
            Lwhich.trim <- any(Lcandidate)
        } else {
            Lwhich.candidate <- max.col(Lcandidate, ties.method = "first")
            Lwhich.trim <- rowSums(Lcandidate) > 0
        }
        trim.start <-
          ifelse(Lwhich.trim, Lwhich.candidate + 1L, min(1L, nchar(subject)))
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
