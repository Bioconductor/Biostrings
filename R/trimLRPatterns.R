### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "trimLRPatterns" new generic.
###

setGeneric("trimLRPatterns", signature = "subject",
    function(Lpattern = "", Rpattern = "", subject,
             max.Lmismatch = rep.int(c(-1, 0), c(nchar(Lpattern) - 1, 1)),
             max.Rmismatch = rep.int(c(-1, 0), c(nchar(Rpattern) - 1, 1)),
             with.Lindels = FALSE, with.Rindels = FALSE,
             Lfixed = TRUE, Rfixed = TRUE, ranges = FALSE)
        standardGeneric("trimLRPatterns")
)

.XString.XStringSet.trimLRPatterns <- 
function(Lpattern = "", Rpattern = "", subject,
         max.Lmismatch = rep.int(c(-1, 0), c(nchar(Lpattern) - 1, 1)),
         max.Rmismatch = rep.int(c(-1, 0), c(nchar(Rpattern) - 1, 1)),
         with.Lindels = FALSE, with.Rindels = FALSE,
         Lfixed = TRUE, Rfixed = TRUE, ranges = FALSE)
{
    if (nchar(Rpattern) == 0) {
        trim.end <- nchar(subject)
    } else {
        max.Rmismatch <- as.integer(max.Rmismatch)
        if (length(max.Rmismatch) < nchar(Rpattern)) {
            max.Rmismatch <-
              c(rep.int(-1L, nchar(Rpattern) - length(max.Rmismatch)), max.Rmismatch)
        }
        if (any(is.na(max.Rmismatch)) || is.unsorted(max.Rmismatch) ||
            length(max.Rmismatch) != nchar(Rpattern))
            stop("'max.Rmismatch' must be a non-decreasing vector of length 'nchar(Rpattern)'")
        with.Rindels <- normargWithIndels(with.Rindels, name = "with.Rindels")
        Rfixed <- normargFixed(Rfixed, class(subject), name = "Rfixed")

        if (length(unique(nchar(subject))) != 1)
            stop("'Rpattern' can only be used with equal length 'subject' strings")
        ncharSubject <- nchar(subject)[1]
        ncharRpattern <- nchar(Rpattern)
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
        trim.start <- min(1L, nchar(subject))
    } else {
        max.Lmismatch <- as.integer(max.Lmismatch)
        if (length(max.Lmismatch) < nchar(Lpattern)) {
            max.Lmismatch <-
              c(rep.int(-1L, nchar(Lpattern) - length(max.Lmismatch)), max.Lmismatch)
        }
        if (any(is.na(max.Lmismatch)) || is.unsorted(max.Lmismatch) ||
            length(max.Lmismatch) != nchar(Lpattern))
            stop("'max.Lmismatch' must be a non-decreasing vector of length 'nchar(Lpattern)'")
        with.Lindels <- normargWithIndels(with.Lindels, name = "with.Lindels")
        Lfixed <- normargFixed(Lfixed, class(subject), name = "Lfixed")
        
        ncharLpattern <- nchar(Lpattern)
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
             max.Lmismatch = rep.int(c(-1, 0), c(nchar(Lpattern) - 1, 1)),
             max.Rmismatch = rep.int(c(-1, 0), c(nchar(Rpattern) - 1, 1)),
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
             max.Lmismatch = rep.int(c(-1, 0), c(nchar(Lpattern) - 1, 1)),
             max.Rmismatch = rep.int(c(-1, 0), c(nchar(Rpattern) - 1, 1)),
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
