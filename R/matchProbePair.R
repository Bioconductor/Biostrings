### matchProbePair simulates a PCR experiment by finding the theoretical
### amplicons mapped to a given probe pair.
matchProbePair <- function(Fprobe, Rprobe, subject, logfile=NULL, verbose=FALSE)
{
    if (!is(subject, "DNAString"))
        stop("'subject' must be a DNAString object")
    # This won't copy the data if Fprobe and Rprobe are already DNAString objects
    F <- DNAString(Fprobe)
    R <- DNAString(Rprobe)

    # F and R hits on the + strand
    Fp_hits <- start(matchPattern(F, subject, fixed=TRUE))
    Rp_hits <- start(matchPattern(R, subject, fixed=TRUE))

    # F and R hits on the - strand
    Fm_hits <- end(matchPattern(reverseComplement(F), subject, fixed=TRUE))
    Rm_hits <- end(matchPattern(reverseComplement(R), subject, fixed=TRUE))

    if (verbose) {
        cat("Fp_hits:", Fp_hits, "  Rp_hits:", Rp_hits,
            "  Fm_hits:", Fm_hits, "  Rm_hits:", Rm_hits, "\n")
    }

    starts0 <- sort(unique(c(Fp_hits, Rp_hits)))
    ends0 <- sort(unique(c(Fm_hits, Rm_hits)))
    nstarts0 <- length(starts0)
    nends0 <- length(ends0)
    i <- j <- 1
    start <- end <- integer(0)
    while (i <= nstarts0 && j <= nends0) {
        if (ends0[j] < starts0[i]) {
            j <- j + 1
            next
        }
        while (i < nstarts0 && starts0[i + 1] <= ends0[j])
            i <- i + 1
        start <- c(start, starts0[i])
        end <- c(end, ends0[j])
        i <- i + 1
        j <- j + 1
    }
    ans <- views(subject, start, end)

    if (!is.null(logfile)) {
        nFp <- length(Fp_hits)
        nRp <- length(Rp_hits)
        nFm <- length(Fm_hits)
        nRm <- length(Rm_hits)
        namp <- length(ans)
        # cat("", ..., sep="\t") is a trick to get an extra tab
        cat("", nFp, nRp, nFm, nRm, namp, file=logfile, sep="\t")
    }
    ans
}

