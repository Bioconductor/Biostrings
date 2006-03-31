
expl <- function(x) strsplit(as.character(x), "")

setGeneric("needwunsQS", function(s1, s2, substmat, gappen=8)
 standardGeneric("needwunsQS"))

setMethod("needwunsQS", c("character", "character", "matrix", "numeric"),
  function(s1, s2, substmat, gappen=8) {
#
# Needleman-Wunsch global alignment following Durbin et al
# QS = quadratic space requirement, simple gap penalty 
#
 if (length(s1) > 1 | length(s2) > 1) stop("currently requires s1, s2 each to be vectors of length 1")
 alp <- colnames(substmat)
 es1 <- expl(s1)[[1]]
 es2 <- expl(s2)[[1]]
 allco <- unique(c(es1,es2))
 allscoco <- colnames(substmat)
 if (!all(allco %in% allscoco)) {
     print("not all symbols used in strings are present in scoring matrix; discrepancies:")
     print(setdiff(allco, allscoco))
     stop("fatal error")
     }
 n1 <- length(es1)
 n2 <- length(es2)
# we are going to create matrices to hold score and traceback info
 sco <- matrix(0, nr=n1+1, nc=n2+1)
 tra <- matrix(0, nr=n1+1, nc=n2+1)
 sco[1,] <- -1*(0:(n2))*gappen
 sco[,1] <- -1*(0:(n1))*gappen
 for (i in 2:(n1+1))
  for (j in 2:(n2+1))
   {
# traceback is 1 if we have a match, 2 (3) indicates gap in aligned string 2 (1)
   tra[i,j] <- which.max(c(sco[i-1,j-1]+substmat[es1[i-1], es2[j-1]], sco[i-1,j]-gappen,
                sco[i,j-1]-gappen))
   sco[i,j] <- max(c(sco[i-1,j-1]+substmat[es1[i-1], es2[j-1]], sco[i-1,j]-gappen,
                sco[i,j-1]-gappen))
   }
 reval1 <- reval2 <- NULL
 fcol <- n2+1
 frow <- n1+1
 while (max(c(fcol,frow)) > 1) {
    crit <- tra[frow,fcol]
    if (crit>0) {
       reval1 <- c(switch(crit,es1[frow-1],es1[frow-1],"-"),reval1)
       reval2 <- c(switch(crit,es2[fcol-1],"-",es2[fcol-1]),reval2)
       frow <- switch(crit,frow-1,frow-1,frow)
       fcol <- switch(crit,fcol-1,fcol,fcol-1)
       }
    else if (fcol > 1) {
       reval2 <- c(es2[fcol-1], reval2)
       reval1 <- c("-", reval1)
       fcol <- fcol - 1 
       }
    else if (frow > 1) {
       reval1 <- c(es1[frow-1], reval1)
       reval2 <- c("-", reval2)
       frow <- frow - 1 
       }
    }
 ans <- c(al1=paste(reval1,collapse=""),
  al2=paste(reval2,collapse=""))
 attr(ans,"scomat") <- sco
 attr(ans,"score") <- sco[n1,n2]
 attr(ans,"tramat") <- tra
#
# we will use S3 class concepts here until we come
# up with a consensus approach to class architecture
#
 class(ans) <- "needwunsQS"
 ans
})

setMethod("needwunsQS", c("character", "character", "matrix", "missing"),
 function(s1, s2, substmat, gappen=8) needwunsQS(s1, s2, substmat, gappen) )

setMethod("needwunsQS", c("BString", "BString", "matrix", "missing"),
 function(s1, s2, substmat, gappen=8) needwunsQS(toString(s1), toString(s2), substmat, gappen) )

setMethod("needwunsQS", c("BString", "BString", "matrix", "numeric"),
  function(s1, s2, substmat, gappen=8) needwunsQS(toString(s1), toString(s2), substmat, gappen) )

alignScore <- function(x, ...) UseMethod("alignScore")
alignScore.needwunsQS <- function(x, ...) attr(x,"score")

print.needwunsQS <- function(x, ...) {
 print(matrix(c(x["al1"], x["al2"]),ncol=1))
 cat(paste("Score: ", alignScore(x), collapse=""),"\n")
 invisible(x)
}

