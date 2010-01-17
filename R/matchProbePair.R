### =========================================================================
### matchProbePair()
### ----------------
###
### -------------------------------------------------------------------------


# NOT USED, NOT EXPORTED
# 'x' and 'y' must be numerics
# Return the _list_ of all (x0, y0) pairs such that
#   d(x0,y0) == min{d(xx,yy), xx in x, yy in y}
closestValues <- function(x, y)
{
    sux <- sort(unique(x))
    suy <- sort(unique(y))
    lsux <- length(sux)
    lsuy <- length(suy)
    x0 <- sux[1]
    y0 <- suy[1]
    d0 <- abs(x0 - y0)
    i <- j <- 1
    repeat {
        if (sux[i] <= suy[j]) i <- i + 1 else j <- j + 1
        if (i > lsux || j > lsuy)
            break
        xx <- sux[i]
        yy <- suy[j]
        d <- abs(xx - yy)
        if (d < d0) {
            x0 <- xx
            y0 <- yy
            d0 <- d
        } else if (d == d0) {
            x0 <- append(x0, xx)
            y0 <- append(y0, yy)
        }
    }
    ans <- list()
    for (k in 1:length(x0))
        ans[[k]] <- c(x0[k], y0[k])
    ans
}

reduceProbePairMatches <- function(start, end)
{
    start <- sort(unique(start))
    end <- sort(unique(end))
    nstart <- length(start)
    nend <- length(end)
    i <- j <- 1
    start0 <- end0 <- integer(0)
    while (i <= nstart && j <= nend) {
        if (end[j] < start[i]) {
            j <- j + 1
            next
        }
        while (i < nstart && start[i + 1] <= end[j])
            i <- i + 1
        start0 <- c(start0, start[i])
        end0 <- c(end0, end[j])
        i <- i + 1
        j <- j + 1
    }
    data.frame(start=start0, end=end0)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "matchProbePair" new generic.
###
### Simulates a PCR experiment by finding the "theoretical amplicons" mapped
### to a given probe pair.
###

setGeneric("matchProbePair", signature="subject",
    function(Fprobe, Rprobe, subject, algorithm="auto",
             logfile=NULL, verbose=FALSE)
        standardGeneric("matchProbePair")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("matchProbePair", "DNAString", 
    function(Fprobe, Rprobe, subject, algorithm="auto",
             logfile=NULL, verbose=FALSE)
    {
        ## This won't copy the data if Fprobe and Rprobe are already DNAString objects
        F <- DNAString(Fprobe)
        R <- DNAString(Rprobe)

        ## F and R hits on the + strand
        Fp_hits <- start(matchPattern(F, subject, algorithm=algorithm))
        Rp_hits <- start(matchPattern(R, subject, algorithm=algorithm))

        ## F and R hits on the - strand
        Fm_hits <- end(matchPattern(reverseComplement(F), subject, algorithm=algorithm))
        Rm_hits <- end(matchPattern(reverseComplement(R), subject, algorithm=algorithm))

        if (verbose) {
            cat("Fp_hits:", Fp_hits, "  Rp_hits:", Rp_hits,
                "  Fm_hits:", Fm_hits, "  Rm_hits:", Rm_hits, "\n")
        }

        matches0 <- reduceProbePairMatches(c(Fp_hits, Rp_hits), c(Fm_hits, Rm_hits))
        ans <- Views(subject, start=matches0$start, end=matches0$end)

        if (!is.null(logfile)) {
            nFp <- length(Fp_hits)
            nRp <- length(Rp_hits)
            nFm <- length(Fm_hits)
            nRm <- length(Rm_hits)
            nmatches0 <- length(ans)
            ## cat("", ..., sep="\t") is a trick to get an extra tab
            cat("", nFp, nRp, nFm, nRm, nmatches0, file=logfile, sep="\t")
        }
        ans
    }
)

### Dispatch on 'subject' (see signature of generic).
### WARNING: Unlike the other "matchProbePair" methods, the XStringViews object
### returned by this method is not guaranteed to have its views ordered from
### left to right in general! One important particular case where this is
### guaranteed though is when 'isNormal(subject)' is TRUE (i.e. 'subject' is
### a normal XStringViews object).
setMethod("matchProbePair", "XStringViews",
    function(Fprobe, Rprobe, subject, algorithm="auto", logfile=NULL, verbose=FALSE)
    {
        ans_start <- ans_width <- integer(0)
        for (i in seq_len(length(subject))) {
            pm <- matchProbePair(Fprobe, Rprobe, subject[[i]],
                                 algorithm=algorithm, logfile=logfile, verbose=verbose)
            offset <- start(subject)[i] - 1L
            ans_start <- c(ans_start, offset + start(pm))
            ans_width <- c(ans_width, width(pm))
        }
        unsafe.newXStringViews(subject(subject), ans_start, ans_width)
    }
)

setMethod("matchProbePair", "MaskedDNAString",
    function(Fprobe, Rprobe, subject, algorithm="auto", logfile=NULL, verbose=FALSE)
        matchProbePair(Fprobe, Rprobe, toXStringViewsOrXString(subject),
                       algorithm=algorithm, logfile=logfile, verbose=verbose)
)

