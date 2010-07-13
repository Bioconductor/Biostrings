strs <- c(string1 = "AAGGTCTCCA-GCCTGCCCTTCAGTGTGGAGGCGCTCATG--TCGGACA",
          string2 = "AAGGTCTCCA-GCCTGCCCTTCAGCGTGGAGGCGCTCATG--TCCGACA",
          string3 = "CATTTATATATGGTCCCCCTCCCCCCAAGAAACACACATAGTTTTGACA")

BSet <- BStringSet(strs)
BMSASet <- MultipleAlignment(strs)

## I can't use checkEquals on the BSet or BMSASet objects :(
## checkEquals(BSet, BMSASet@set)  


## check that we are enforcing the lengths being the same
strsBad <- strs
strsBad[1] <- narrow(strs[1], end=9)
checkException(MultipleAlignment(strsBad))


## check my new subset function
strs2 <- c(string1 = "AAGGTCTCCAGCCTGCCCTTCAGTGTGGAGGCGCTCATGTCGGACA",
           string2 = "AAGGTCTCCAGCCTGCCCTTCAGCGTGGAGGCGCTCATGTCCGACA",
           string3 = "CATTTATATAGGTCCCCCTCCCCCCAAGAAACACACATATTTGACA")

checkEquals(as.character(subsetColumns(BMSASet,
                                       start=c(1,12,43),
                                       end=c(10,40,49))),
            as.character(MultipleAlignment(strs2)))


## check consensus methods
conStr <- "AAGGTCTCCA-GCCTGCCCTTCAGCGTGGAGGCGCTCATG--TC?GACA"
checkEquals(as.character(consensusString(BMSASet)),
            as.character(conStr))

conMat <- matrix(c(2,3,0, 1,0,0, 0,0,2, 0,0,1), nrow=4, ncol=3, byrow=TRUE)
rownames(conMat) <- c("A","C","G","T")
checkEquals(consensusMatrix(narrow(BMSASet,end=3)),
            conMat)

