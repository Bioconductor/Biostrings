#
# naive implementation of some algorithms in Durbin et al
# includes definition of simpleBioString that is just for R character vbls
#

# note by VJC, 10 nov 05: The aminoacid alphabet doesn't seem to work
# along the same lines as the nucleotide alphabet.  so i've introduced
# a simpleBioString class that is based on R strings.  we can still use
# the alphabet structure if we want, and define, e.g., simpAAString to 
# extend simpleBioString and verify that the alphabet is honored.
# but i haven't done that yet

StandardPeptideAlphabet <- function()
     new("AminoAcidAlphabet",
         letters= c('A', # alanine
                    'B', # aspartatic acid or asparagine
                    'C', # cysteine
                    'D', # aspartic acid
                    'E', # glutamic acid
                    'F', # phenylalanine
                    'G', # glycine
                    'H', # histidine
                    'I', # isoleucine
                    'K', # lysine
                    'L', # leucine
                    'M', # methionine
                    'N', # asparagine
                    'P', # proline
                    'Q', # glutamine
                    'R', # arginine
                    'S', # serine
                    'T', # threonine
                    'U', # selenocysteine
                    'V', # valine
                    'W', # tryptophan
                    'Y', # tyrosine
                    'Z', # glutamaic acid or glutamine
                    '*'  # translation stop
                    ))

setClass("simpleBioString", representation(values="character"),
 contains=c("BioString"), prototype=list(values="", offsets=matrix(NA,1,1), alphabet=StandardPeptideAlphabet(),
  initialized=FALSE))

setMethod("show", "simpleBioString", 
          function (object)
      {
              cat("    Object of class simpleBioString with\n")
              if (is(object@alphabet, "BioPatternAlphabet"))
                  cat("Pattern alphabet: ")
              else cat("Alphabet: ")
              cat(object@alphabet@letters, "\n")
	      print(object@values)
      })


AAString <- function(src, alphabet=StandardPeptideAlphabet(), gap=alphabet@gap) {
 new("simpleBioString", values=src, alphabet=alphabet, gap=gap, initialized=TRUE,
   offsets=matrix(c(0,0),nc=2)) }

setMethod("substr", "simpleBioString", function(x, start, stop) {
 new("simpleBioString", values=substr(x@values, start, stop), alphabet=x@alphabet)
})

setAs("simpleBioString", "character", function(from, to) {
 from@values })

setMethod("[",
          signature(x = "simpleBioString"),
          function (x, i, j, ..., drop)
      {
          if (!missing(j)) stop("only single indexing permitted")
          as(x,"character")[i]
      })

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

setMethod("needwunsQS", c("BioString", "BioString", "matrix", "missing"),
 function(s1, s2, substmat, gappen=8) needwunsQS(s1, s2, substmat, gappen) )

setMethod("needwunsQS", c("BioString", "BioString", "matrix", "numeric"),
  function(s1, s2, substmat, gappen=8) needwunsQS( s1@values, s2@values, substmat, gappen) )

alignScore <- function(x, ...) UseMethod("alignScore")
alignScore.needwunsQS <- function(x, ...) attr(x,"score")

print.needwunsQS <- function(x, ...) {
 print(matrix(c(x["al1"], x["al2"]),ncol=1))
 cat(paste("Score: ", alignScore(x), collapse=""),"\n")
 invisible(x)
}

  


#data(blosum50)
#needwuns( "PAWHEAE", "HEAGAWGHEE", blosum50 )


rpeptide <- function(nchar, base="ARNDCQEGHILKMFPSTWYVBZX" ) {
# this generates a random peptide, the base alphabet is as in BLOSUM50
 lets <- expl(base)[[1]]
 lets <- lets[-c(1,length(lets))]
 paste(sample(lets,size=nchar,replace=TRUE), collapse="")
}
 
