strings_DNAMultipleAlignment <- function()
{
    c(string1 = "AAGGTCTCCA-GCCTGCCCTTCAGTGTGGAGGCGCTCATG--TCGGACA",
      string2 = "AAGGTCTCCA-GCCTGCCCTTCAGCGTGGAGGCGCTCATG--TCCGACA",
      string3 = "CATTTATATATGGTCCCCCTCCCCCCAAGAAACACACATAGTTTTGACA")
}

make_DNAMultipleAlignment <- function()
{
    DNAMultipleAlignment(strings_DNAMultipleAlignment())
}

test_DNAMultipleAlignment_empty <- function()
{
    malign <- DNAMultipleAlignment()
    checkTrue(validObject(malign, test=TRUE, complete=TRUE))
    checkIdentical(as.character(unmasked(malign)), as.character(DNAStringSet()))
    checkIdentical(rownames(malign), NULL)
    checkIdentical(rowmask(malign), new("NormalIRanges"))
    checkIdentical(colmask(malign), new("NormalIRanges"))
    checkIdentical(as.character(maskMotif(malign, "GC")), character())
    checkIdentical(as.character(maskGaps(malign)), character())
    checkIdentical(nrow(malign), 0L)
    checkIdentical(ncol(malign), 0L)
    checkIdentical(dim(malign), c(0L, 0L))
    checkIdentical(maskednrow(malign), 0L)
    checkIdentical(maskedncol(malign), 0L)
    checkIdentical(maskeddim(malign), c(0L, 0L))
    checkIdentical(maskedratio(malign), c(0L/0L, 0L/0L))
    checkIdentical(nchar(malign), 0L)
    checkIdentical(seqtype(malign), "DNA")
    checkIdentical(as.character(malign), character(0))
    checkIdentical(consensusMatrix(malign),
                   matrix(integer(), nrow=length(DNA_ALPHABET), ncol=0,
                          dimnames=list(DNA_ALPHABET, NULL)))
    checkIdentical(consensusString(malign), character())
    checkIdentical(as.character(consensusViews(malign)), character())
    checkIdentical(alphabetFrequency(malign),
                   matrix(integer(), nrow=0, ncol=length(DNA_ALPHABET),
                          dimnames=list(NULL, DNA_ALPHABET)))
    checkIdentical(alphabetFrequency(malign, collapse=TRUE),
                   structure(integer(length(DNA_ALPHABET)), names=DNA_ALPHABET))
}

test_DNAMultipleAlignment_unnamed <- function()
{
    malign <- make_DNAMultipleAlignment()
    rownames(malign) <- NULL
    checkTrue(validObject(malign, test=TRUE))
    checkIdentical(as.character(unmasked(malign)),
                   unname(strings_DNAMultipleAlignment()))
    checkIdentical(rownames(malign), NULL)
    checkIdentical(rowmask(malign), new("NormalIRanges"))
    checkIdentical(colmask(malign), new("NormalIRanges"))
    checkIdentical(as.character(maskMotif(malign, "GC", fixed=FALSE)),
                   c("AAGGTCTCCA-TCCTTCAGTGGAGTCATG--TCGGACA",
                     "AAGGTCTCCA-TCCTTCAGTGGAGTCATG--TCCGACA",
                     "CATTTATATATCCCTCCCCAAGAAACATAGTTTTGACA"))
    checkIdentical(as.character(maskGaps(malign)),
                   unname(strings_DNAMultipleAlignment()))
    checkIdentical(nrow(malign), length(strings_DNAMultipleAlignment()))
    checkIdentical(ncol(malign), nchar(strings_DNAMultipleAlignment())[[1]])
    checkIdentical(dim(malign),
                   c(length(strings_DNAMultipleAlignment()),
                     nchar(strings_DNAMultipleAlignment())[[1]]))
    checkIdentical(maskednrow(malign), 0L)
    checkIdentical(maskedncol(malign), 0L)
    checkIdentical(maskeddim(malign), c(0L, 0L))
    checkIdentical(maskedratio(malign), c(0, 0))
    checkIdentical(nchar(malign), nchar(strings_DNAMultipleAlignment())[[1]])
    checkIdentical(seqtype(malign), "DNA")
    checkIdentical(as.character(malign), unname(strings_DNAMultipleAlignment()))
    checkIdentical(consensusMatrix(malign)[1:4, 1:4],
                   rbind(A=c(2L,3L,0L,0L),
                         C=c(1L,0L,0L,0L),
                         G=c(0L,0L,2L,2L),
                         T=c(0L,0L,1L,1L)))
    checkIdentical(consensusString(malign),
                   "MAKKTMTMYA-GSYYSCCCTYCMSYSWRGARRCRCWCATR--TYBGACA")
    checkIdentical(as.character(consensusViews(malign)),
                   "MAKKTMTMYA-GSYYSCCCTYCMSYSWRGARRCRCWCATR--TYBGACA")
    checkIdentical(alphabetFrequency(malign)[,1:4],
                   cbind(A=c(8L,8L,15L),
                         C=c(14L,16L,16L),
                         G=c(14L,13L,5L),
                         T=c(10L,9L,13L)))
    checkIdentical(alphabetFrequency(malign, collapse=TRUE)[1:4],
                   c(A=31L, C=46L, G=32L, T=32L))
}

test_DNAMultipleAlignment_named <- function()
{
    malign <- make_DNAMultipleAlignment()
    checkTrue(validObject(malign, test=TRUE))
    checkIdentical(as.character(unmasked(malign)), strings_DNAMultipleAlignment())
    checkIdentical(rownames(malign), names(strings_DNAMultipleAlignment()))
    checkIdentical(rowmask(malign), new("NormalIRanges"))
    checkIdentical(colmask(malign), new("NormalIRanges"))
    checkIdentical(as.character(maskMotif(malign, "GC", fixed=FALSE)),
                   c(string1="AAGGTCTCCA-TCCTTCAGTGGAGTCATG--TCGGACA",
                     string2="AAGGTCTCCA-TCCTTCAGTGGAGTCATG--TCCGACA",
                     string3="CATTTATATATCCCTCCCCAAGAAACATAGTTTTGACA"))
    checkIdentical(as.character(maskGaps(malign)),
                   strings_DNAMultipleAlignment())
    checkIdentical(nrow(malign), length(strings_DNAMultipleAlignment()))
    checkIdentical(ncol(malign), nchar(strings_DNAMultipleAlignment())[[1]])
    checkIdentical(dim(malign),
                   c(length(strings_DNAMultipleAlignment()),
                     nchar(strings_DNAMultipleAlignment())[[1]]))
    checkIdentical(maskednrow(malign), 0L)
    checkIdentical(maskedncol(malign), 0L)
    checkIdentical(maskeddim(malign), c(0L, 0L))
    checkIdentical(maskedratio(malign), c(0, 0))
    checkIdentical(nchar(malign), nchar(strings_DNAMultipleAlignment())[[1]])
    checkIdentical(seqtype(malign), "DNA")
    checkIdentical(as.character(malign), strings_DNAMultipleAlignment())
    checkIdentical(consensusMatrix(malign)[1:4, 1:4], 
                   rbind(A=c(2L,3L,0L,0L),
                         C=c(1L,0L,0L,0L),
                         G=c(0L,0L,2L,2L),
                         T=c(0L,0L,1L,1L)))
    checkIdentical(consensusString(malign),
                   "MAKKTMTMYA-GSYYSCCCTYCMSYSWRGARRCRCWCATR--TYBGACA")
    checkIdentical(as.character(consensusViews(malign)),
                   "MAKKTMTMYA-GSYYSCCCTYCMSYSWRGARRCRCWCATR--TYBGACA")
    checkIdentical(alphabetFrequency(malign)[,1:4],
                   cbind(A=c(8L,8L,15L),
                         C=c(14L,16L,16L),
                         G=c(14L,13L,5L),
                         T=c(10L,9L,13L)))
    checkIdentical(alphabetFrequency(malign, collapse=TRUE)[1:4],
                   c(A=31L, C=46L, G=32L, T=32L))
}

### FIXME!!
BROKEN_test_DNAMultipleAlignment_mask_some <- function()
{
    malign <- make_DNAMultipleAlignment()
    rowmask(malign) <- IRanges(2,2)
    colmask(malign) <- IRanges(c(1,21,43), c(10,35,49))
    checkTrue(validObject(malign, test=TRUE))
    checkIdentical(as.character(unmasked(malign)), strings_DNAMultipleAlignment())
    checkIdentical(rownames(malign), names(strings_DNAMultipleAlignment()))
    checkIdentical(rowmask(malign), asNormalIRanges(IRanges(2,2)))
    checkIdentical(colmask(malign), asNormalIRanges(IRanges(c(1,21,43), c(10,35,49))))
    checkIdentical(as.character(maskMotif(malign, "GC", fixed=FALSE)),
                   c(string1="-TCCTTCATG--", string3="TCCCTACATAGT"))
    checkIdentical(as.character(maskGaps(malign, min.block.width=1)),
                   c(string1="GCCTGCCCTTCATG", string3="GGTCCCCCTACATA"))
    checkIdentical(nrow(malign), length(strings_DNAMultipleAlignment()))
    checkIdentical(ncol(malign), nchar(strings_DNAMultipleAlignment())[[1]])
    checkIdentical(dim(malign),
                   c(length(strings_DNAMultipleAlignment()),
                     nchar(strings_DNAMultipleAlignment())[[1]]))
    checkIdentical(maskednrow(malign), 1L)
    checkIdentical(maskedncol(malign), 32L)
    checkIdentical(maskeddim(malign), c(1L, 32L))
    checkIdentical(maskedratio(malign), c(1/3, 32/49))
    checkIdentical(nchar(malign), 17L)
    checkIdentical(seqtype(malign), "DNA")
    checkIdentical(as.character(malign),
                   c(string1="-GCCTGCCCTTCATG--", string3="TGGTCCCCCTACATAGT"))
    checkIdentical(consensusMatrix(malign)[1:4, 1:4], 
                   rbind(A=rep(NA_integer_,4),
                         C=rep(NA_integer_,4),
                         G=rep(NA_integer_,4),
                         T=rep(NA_integer_,4)))
    checkIdentical(consensusString(malign),
                   "##########TGSYYSCCCT###############WCATRGT#######")
    checkIdentical(as.character(consensusViews(malign)),
                   c("TGSYYSCCCT", "WCATRGT"))
    checkIdentical(alphabetFrequency(malign)[,1:4],
                   cbind(A=c(1L,NA,3L),
                         C=c(6L,NA,6L),
                         G=c(3L,NA,3L),
                         T=c(4L,NA,5L)))
    checkIdentical(alphabetFrequency(malign, collapse=TRUE)[1:4],
                   c(A=4L, C=12L, G=6L, T=9L))
}

### FIXME!!
BROKEN_test_DNAMultipleAlignment_mask_all_rows <- function()
{
    malign <- make_DNAMultipleAlignment()
    rowmask(malign) <- IRanges(1,3)
    colmask(malign) <- IRanges(c(1,21,43), c(10,35,49))
    checkTrue(validObject(malign, test=TRUE))
    checkIdentical(as.character(unmasked(malign)), strings_DNAMultipleAlignment())
    checkIdentical(rownames(malign), names(strings_DNAMultipleAlignment()))
    checkIdentical(rowmask(malign), asNormalIRanges(IRanges(1,3)))
    checkIdentical(colmask(malign), asNormalIRanges(IRanges(c(1,21,43), c(10,35,49))))
    checkIdentical(as.character(maskMotif(malign, "GC", fixed=FALSE)), character())
    checkIdentical(as.character(maskGaps(malign)), character())
    checkIdentical(nrow(malign), length(strings_DNAMultipleAlignment()))
    checkIdentical(ncol(malign), nchar(strings_DNAMultipleAlignment())[[1]])
    checkIdentical(dim(malign),
                   c(length(strings_DNAMultipleAlignment()),
                     nchar(strings_DNAMultipleAlignment())[[1]]))
    checkIdentical(maskednrow(malign), 3L)
    checkIdentical(maskedncol(malign), 32L)
    checkIdentical(maskeddim(malign), c(3L, 32L))
    checkIdentical(maskedratio(malign), c(3/3, 32/49))
    checkIdentical(nchar(malign), 17L)
    checkIdentical(seqtype(malign), "DNA")
    checkIdentical(as.character(malign), character(0))
    checkIdentical(consensusMatrix(malign)[1:4, 1:4], 
                   rbind(A=rep(NA_integer_,4),
                         C=rep(NA_integer_,4),
                         G=rep(NA_integer_,4),
                         T=rep(NA_integer_,4)))
    checkIdentical(consensusString(malign),
                   "#################################################")
    checkIdentical(alphabetFrequency(malign)[,1:4],
                   cbind(A=rep(NA_integer_, 3),
                         C=rep(NA_integer_, 3),
                         G=rep(NA_integer_, 3),
                         T=rep(NA_integer_, 3)))
    checkIdentical(alphabetFrequency(malign, collapse=TRUE)[1:4],
                   c(A=0L, C=0L, G=0L, T=0L))
}
