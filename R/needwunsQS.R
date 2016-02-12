### =========================================================================
### The needwunsQS() generic & related functions
### --------------------------------------------
###
### Needleman-Wunsch global alignment following Durbin et al
### QS = quadratic space requirement, simple gap penalty 
###
### -------------------------------------------------------------------------


### NOT USED ANYMORE!
### R-implementation of the Needleman-Wunsch algo originally written by VJC
### for Biostrings 1 and adapted for Biostrings 2 by H. Pag\`es.
### 's1' and 's2' must be character vectors of length 1.
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
    sco <- tra <- matrix(0, nrow=n1+1, ncol=n2+1)
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


print.needwunsQS <- function(x, ...)
{
    print(matrix(c(x["al1"], x["al2"]),ncol=1))
    cat(paste("Score: ", score(x), collapse=""),"\n")
    invisible(x)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### C-implementation of the Needleman-Wunsch algo.
###
### The various "needwunsQS" methods below don't call
### .Call2("align_needwunsQS", ...) directly but call XString.needwunsQS()
### instead which itself calls .Call2("align_needwunsQS", ...).
### Some quick testing shows that, depending on the size of the strings to
### align, this C version is 100 to 1000 times faster than the above
### .needwunsQS().
### 's1' and 's2' must be XString objects of same sequence type.
### Return a PairwiseAlignments object where the "al1" and "al2" slots
### contain the aligned versions of 's1' and 's2'.
XString.needwunsQS <- function(s1, s2, substmat, gappen)
{
    .Deprecated("pairwiseAlignment")
    if (class(s1) != class(s2))
        stop("'s1' and 's2' are not of the same class")
    if (is.character(substmat)) {
        if (length(substmat) != 1)
            stop("'substmat' is a character vector of length != 1")
        matname <- substmat
        substmat <- try(getdata(matname), silent=TRUE)
        if (is(substmat, "try-error"))
            stop("unknown scoring matrix \"", matname, "\"")
    }
    if (!is.matrix(substmat) || !is.integer(substmat))
        stop("'substmat' must be a matrix of integers")
    if (!identical(rownames(substmat), colnames(substmat)))
        stop("row and column names differ for matrix 'substmat'")
    if (is.null(rownames(substmat)))
        stop("matrix 'substmat' must have row and column names")
    if (any(duplicated(rownames(substmat))))
        stop("matrix 'substmat' has duplicated row names")
    if (is.null(xscodec(s1))) {
        codes <- as.integer(charToRaw(paste(rownames(substmat), collapse="")))
        gap_code <- charToRaw("-")
    } else {
        if (!all(rownames(substmat) %in% alphabet(s1)))
            stop("matrix 'substmat' is incompatible with 's1' alphabet")
        letters2codes <- xscodec(s1)@codes
        names(letters2codes) <- xscodec(s1)@letters
        codes <- letters2codes[rownames(substmat)]
        gap_code <- as.raw(letters2codes[["-"]])
    }
    lkup <- buildLookupTable(codes, 0:(nrow(substmat)-1))
    C_ans <- .Call2("align_needwunsQS",
                   s1, s2,
                   substmat, nrow(substmat), lkup,
                   as.integer(gappen), gap_code,
                   PACKAGE="Biostrings")
    PairwiseAlignments(new(class(s1), shared = C_ans$al1, length = length(C_ans$al1)), 
                       new(class(s2), shared = C_ans$al2, length = length(C_ans$al2)),
                       type = "global", substitutionMatrix = substmat,
                       gapOpening = 0, gapExtension = gappen)
}

setGeneric("needwunsQS", signature=c("s1", "s2"),
    function(s1, s2, substmat, gappen=8) standardGeneric("needwunsQS")
)
setMethod("needwunsQS", signature(s1="character", s2="character"),
    function(s1, s2, substmat, gappen)
        XString.needwunsQS(BString(s1), BString(s2), substmat, gappen)
)
setMethod("needwunsQS", signature(s1="character", s2="XString"),
    function(s1, s2, substmat, gappen)
        XString.needwunsQS(XString(seqtype(s2), s1), s2, substmat, gappen)
)
setMethod("needwunsQS", signature(s1="XString", s2="character"),
    function(s1, s2, substmat, gappen)
        XString.needwunsQS(s1, XString(seqtype(s1), s2), substmat, gappen)
)
setMethod("needwunsQS", signature(s1="XString", s2="XString"),
    function(s1, s2, substmat, gappen)
        XString.needwunsQS(s1, s2, substmat, gappen)
)

