### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "trimLRPatterns" new generic.
###

setGeneric("trimLRPatterns", signature = "subject",
    function(Lpattern = NULL, Rpattern = NULL, subject,
             max.Lmismatch = 0, max.Rmismatch = 0,
             with.Lindels = FALSE, with.Rindels = FALSE,
             Lfixed = TRUE, Rfixed = TRUE, ranges = FALSE)
        standardGeneric("trimLRPatterns")
)

.XString.XStringSet.trimLRPatterns <- 
function(Lpattern = NULL, Rpattern = NULL, subject,
         max.Lmismatch = 0, max.Rmismatch = 0,
         with.Lindels = FALSE, with.Rindels = FALSE,
         Lfixed = TRUE, Rfixed = TRUE, ranges = FALSE)
{
    whichMin <- function(x) {
        if (is.matrix(x))
            max.col(- x, ties.method = "first")
        else
            which.min(x)
    }

    if (is.null(Rpattern)) {
        trim.end <- nchar(subject)
    } else {
        max.Rmismatch <- normargMaxMismatch(max.Rmismatch, name = "max.Rmismatch")
        with.Rindels <- normargWithIndels(with.Rindels, name = "with.Rindels")
        Rfixed <- normargFixed(Rfixed, class(subject), name = "Rfixed")

        if (length(unique(nchar(subject))) != 1)
            stop("'Rpattern' can only be used with equal length 'subject' strings")
        ncharSubject <- nchar(subject)[1]
        ncharRpattern <- nchar(Rpattern)
        max.Rmismatch <- min(ncharRpattern, max.Rmismatch)
        Rstarting.at <-
          (ncharSubject - ncharRpattern + 1L + max.Rmismatch):
          (ncharSubject - ncharRpattern + 1L)
        Rnedit <-
          neditStartingAt(Rpattern, subject, starting.at = Rstarting.at,
                          with.indels = with.Rindels, fixed = Rfixed)
        Rwhich.min <- whichMin(Rnedit)
        if (is(subject, "XString")) {
            Rwhich.trim <- (Rnedit[Rwhich.min] <= max.Rmismatch)
        } else {
            Rwhich.minVector <-
              seq_len(length(subject)) + (Rwhich.min - 1L) * length(subject)
            Rwhich.trim <- (Rnedit[Rwhich.minVector] <= max.Rmismatch)
        }
        trim.end <-
          ifelse(Rwhich.trim, Rstarting.at[Rwhich.min] - 1L, nchar(subject))
    }

    if (is.null(Lpattern)) {
        trim.start <- min(1L, nchar(subject))
    } else {
        max.Lmismatch <- normargMaxMismatch(max.Lmismatch, name = "max.Lmismatch")
        with.Lindels <- normargWithIndels(with.Lindels, name = "with.Lindels")
        Lfixed <- normargFixed(Lfixed, class(subject), name = "Lfixed")
        
        ncharLpattern <- nchar(Lpattern)
        max.Lmismatch <- min(ncharLpattern, max.Lmismatch)
        Lending.at <- (ncharLpattern - max.Lmismatch):ncharLpattern
        Lnedit <-
          neditEndingAt(Lpattern, subject, ending.at = Lending.at,
                        with.indels = with.Lindels, fixed = Lfixed)
        Lwhich.min <- whichMin(Lnedit)
        if (is(subject, "XString")) {
            Lwhich.trim <- (Lnedit[Lwhich.min] <= max.Lmismatch)
        } else {
            Lwhich.minVector <-
              seq_len(length(subject)) + (Lwhich.min - 1L) * length(subject)
            Lwhich.trim <- (Lnedit[Lwhich.minVector] <= max.Lmismatch)
        }
        trim.start <-
          ifelse(Lwhich.trim, Lending.at[Lwhich.min] + 1L,
                 min(1L, nchar(subject)))
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
    function(Lpattern = NULL, Rpattern = NULL, subject,
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
    function(Lpattern = NULL, Rpattern = NULL, subject,
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
