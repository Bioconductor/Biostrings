### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "matchLRPatterns" new generic.
###

setGeneric(
    "matchLRPatterns", signature="subject",
    function(Lpattern, Rpattern, max.ngaps, subject, Lmismatch=0, Rmismatch=0,
             Lfixed=TRUE, Rfixed=TRUE)
        standardGeneric("matchLRPatterns")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("matchLRPatterns", "BString", 
    function(Lpattern, Rpattern, max.ngaps, subject, Lmismatch=0, Rmismatch=0,
             Lfixed=TRUE, Rfixed=TRUE)
    {
        ans_start <- ans_end <- integer(0)
        Lmatches <- matchPattern(Lpattern, subject, mismatch=Lmismatch, fixed=Lfixed)
        if (length(Lmatches) != 0L) {
            Rmatches <- matchPattern(Rpattern, subject, mismatch=Rmismatch, fixed=Rfixed)
            if (length(Rmatches) != 0L) {
                for (i in seq_len(length(Lmatches))) {
                    ngaps <- start(Rmatches) - end(Lmatches)[i] - 1L
                    jj <- which(0L <= ngaps & ngaps <= max.ngaps)
                    ans_start <- c(ans_start, rep.int(start(Lmatches)[i], length(jj)))
                    ans_end <- c(ans_end, end(Rmatches)[jj])
                }
            }
        }
        views(subject, ans_start, ans_end)
    }
)

