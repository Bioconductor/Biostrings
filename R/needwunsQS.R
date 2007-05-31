### Originally written by VJC for Biostrings 1.
### Adapted for Biostrings 2 by H. Pages.

### Needleman-Wunsch global alignment following Durbin et al
### QS = quadratic space requirement, simple gap penalty 

### 's1' and 's2' must be character vectors of length 1
.needwunsQS <- function(s1, s2, substmat, gappen)
{
    if (length(s1) != 1 || length(s2) != 1)
        stop("currently requires s1, s2 each to be vectors of length 1")
    es1 <- safeExplode(s1)
    es2 <- safeExplode(s2)
    allco <- unique(c(es1, es2))
    allscoco <- colnames(substmat)
    if (!all(allco %in% allscoco)) {
        print("not all symbols used in strings are present in scoring matrix; discrepancies:")
        print(setdiff(allco, allscoco))
        stop("fatal error")
    }
    n1 <- length(es1)
    n2 <- length(es2)
    ## we are going to create matrices to hold score and traceback info
    sco <- tra <- matrix(0, nr=n1+1, nc=n2+1)
    sco[1,] <- -1*(0:n2)*gappen
    sco[,1] <- -1*(0:n1)*gappen
    for (i in seq_len(n1))
      for (j in seq_len(n2)) {
        ## traceback is 1 if we have a match, 2 (3) indicates gap
        ## in aligned string 2 (1)
        tmp <- c(sco[i  ,j  ] + substmat[es1[i], es2[j]],
                 sco[i  ,j+1] - gappen,
                 sco[i+1,j  ] - gappen)
        tra[i+1,j+1] <- which.max(tmp)
        sco[i+1,j+1] <- max(tmp)
      }

    reval1 <- reval2 <- NULL
    fcol <- n2+1
    frow <- n1+1
    while (max(c(fcol, frow)) > 1) {
        crit <- tra[frow, fcol]
        if (crit > 0) {
            tmp1 <- es1[frow-1]
            reval1 <- c(switch(crit, tmp1, tmp1, "-"), reval1)
            tmp2 <- es2[fcol-1]
            reval2 <- c(switch(crit, tmp2, "-", tmp2), reval2)
            frow <- switch(crit, frow-1, frow-1, frow)
            fcol <- switch(crit, fcol-1, fcol, fcol-1)
        } else if (fcol > 1) {
            reval2 <- c(es2[fcol-1], reval2)
            reval1 <- c("-", reval1)
            fcol <- fcol - 1
        } else if (frow > 1) {
            reval1 <- c(es1[frow-1], reval1)
            reval2 <- c("-", reval2)
            frow <- frow - 1
        }
    }
    ans <- c(al1=paste(reval1, collapse=""), al2=paste(reval2, collapse=""))
    attr(ans,"scomat") <- sco
    attr(ans,"score") <- sco[n1+1, n2+1]
    attr(ans,"tramat") <- tra
    ## we will use S3 class concepts here until we come
    ## up with a consensus approach to class architecture
    class(ans) <- "needwunsQS"
    ans
}

setGeneric(
    "needwunsQS", signature=c("s1", "s2"),
    function(s1, s2, substmat, gappen=8) standardGeneric("needwunsQS")
)
setMethod("needwunsQS", signature(s1="character", s2="character"),
    function(s1, s2, substmat, gappen)
        .needwunsQS(s1, s2, substmat, gappen)
)
setMethod("needwunsQS", signature(s1="character", s2="BString"),
    function(s1, s2, substmat, gappen)
        .needwunsQS(s1, as.character(s2), substmat, gappen)
)
setMethod("needwunsQS", signature(s1="BString", s2="character"),
    function(s1, s2, substmat, gappen)
        .needwunsQS(as.character(s1), s2, substmat, gappen)
)
setMethod("needwunsQS", signature(s1="BString", s2="BString"),
    function(s1, s2, substmat, gappen)
        .needwunsQS(as.character(s1), as.character(s2), substmat, gappen)
)

alignScore <- function(x, ...) UseMethod("alignScore")
alignScore.needwunsQS <- function(x, ...) attr(x,"score")

print.needwunsQS <- function(x, ...)
{
    print(matrix(c(x["al1"], x["al2"]),ncol=1))
    cat(paste("Score: ", alignScore(x), collapse=""),"\n")
    invisible(x)
}

