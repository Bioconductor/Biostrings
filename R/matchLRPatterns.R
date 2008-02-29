### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "matchLRPatterns" new generic.
###

setGeneric("matchLRPatterns", signature="subject",
    function(Lpattern, Rpattern, max.ngaps, subject, max.Lmismatch=0, max.Rmismatch=0,
             Lfixed=TRUE, Rfixed=TRUE)
        standardGeneric("matchLRPatterns")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("matchLRPatterns", "BString", 
    function(Lpattern, Rpattern, max.ngaps, subject, max.Lmismatch=0, max.Rmismatch=0,
             Lfixed=TRUE, Rfixed=TRUE)
    {
        ans_start <- ans_end <- integer(0)
        Lmatches <- matchPattern(Lpattern, subject, max.mismatch=max.Lmismatch, fixed=Lfixed)
        if (length(Lmatches) != 0L) {
            Rmatches <- matchPattern(Rpattern, subject, max.mismatch=max.Rmismatch, fixed=Rfixed)
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

### Dispatch on 'subject' (see signature of generic).
### WARNING: Unlike the other "matchLRPatterns" methods, the BStringViews object
### returned by this method is not guaranteed to have its views ordered from
### left to right in general! One important particular case where this is
### guaranteed though is when 'subject' is a normalized BStringViews object
### and 'max.Lmismatch=0' and 'max.Rmismatch=0' so there are no "out of limits"
### matches.
setMethod("matchLRPatterns", "BStringViews",
    function(Lpattern, Rpattern, max.ngaps, subject, max.Lmismatch=0, max.Rmismatch=0,
             Lfixed=TRUE, Rfixed=TRUE)
    {
        ans_start <- ans_nchar <- integer(0)
        for (i in seq_len(length(subject))) {
            pm <- matchLRPatterns(Lpattern, Rpattern, max.ngaps, subject[[i]],
                                   max.Lmismatch=max.Lmismatch, max.Rmismatch=max.Rmismatch,
                                   Lfixed=Lfixed, Rfixed=Rfixed)
            offset <- start(subject)[i] - 1L
            ans_start <- c(ans_start, offset + start(pm))
            ans_nchar <- c(ans_nchar, nchar(pm))
        }
        new("BStringViews", subject=subject(subject),
            start=ans_start, nchar=ans_nchar, check=FALSE)
    }
)

