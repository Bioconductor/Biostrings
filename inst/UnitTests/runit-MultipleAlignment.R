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
    checkTrue(validObject(malign, test=TRUE))
    checkIdentical(as.character(unmasked(malign)), as.character(DNAStringSet()))
    checkIdentical(rownames(malign), NULL)
    checkIdentical(rowmask(malign), new("NormalIRanges"))
    checkIdentical(colmask(malign), new("NormalIRanges"))
    checkIdentical(nrow(malign), 0L)
    checkIdentical(ncol(malign), 0L)
    checkIdentical(dim(malign), c(0L, 0L))
    checkIdentical(maskednrow(malign), 0L)
    checkIdentical(maskedncol(malign), 0L)
    checkIdentical(maskeddim(malign), c(0L, 0L))
    checkIdentical(maskedratio(malign), c(0L/0L, 0L/0L))
    checkIdentical(nchar(malign), 0L)
    checkIdentical(xsbasetype(malign), "DNA")
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
    checkIdentical(xsbasetype(malign), "DNA")
}

test_DNAMultipleAlignment_named <- function()
{
    malign <- make_DNAMultipleAlignment()
    checkTrue(validObject(malign, test=TRUE))
    checkIdentical(as.character(unmasked(malign)), strings_DNAMultipleAlignment())
    checkIdentical(rownames(malign), names(strings_DNAMultipleAlignment()))
    checkIdentical(rowmask(malign), new("NormalIRanges"))
    checkIdentical(colmask(malign), new("NormalIRanges"))
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
    checkIdentical(xsbasetype(malign), "DNA")
}

test_DNAMultipleAlignment_mask_some <- function()
{
    malign <- make_DNAMultipleAlignment()
    rowmask(malign) <- IRanges(2,2)
    colmask(malign) <- IRanges(c(1,21,43), c(10,35,49))
    checkTrue(validObject(malign, test=TRUE))
    checkIdentical(as.character(unmasked(malign)), strings_DNAMultipleAlignment())
    checkIdentical(rownames(malign), names(strings_DNAMultipleAlignment()))
    checkIdentical(rowmask(malign), asNormalIRanges(IRanges(2,2)))
    checkIdentical(colmask(malign), asNormalIRanges(IRanges(c(1,21,43), c(10,35,49))))
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
    checkIdentical(xsbasetype(malign), "DNA")
}

test_DNAMultipleAlignment_mask_all_rows <- function()
{
    malign <- make_DNAMultipleAlignment()
    rowmask(malign) <- IRanges(1,3)
    colmask(malign) <- IRanges(c(1,21,43), c(10,35,49))
    checkTrue(validObject(malign, test=TRUE))
    checkIdentical(as.character(unmasked(malign)), strings_DNAMultipleAlignment())
    checkIdentical(rownames(malign), names(strings_DNAMultipleAlignment()))
    checkIdentical(rowmask(malign), asNormalIRanges(IRanges(1,3)))
    checkIdentical(colmask(malign), asNormalIRanges(IRanges(c(1,21,43), c(10,35,49))))
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
    checkIdentical(xsbasetype(malign), "DNA")
}
