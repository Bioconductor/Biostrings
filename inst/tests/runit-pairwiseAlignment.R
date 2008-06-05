test_pairwiseAlignment_emptyString <- function()
{
    string1 <- DNAString("")
    string2 <- DNAString("ACGT")
    alignment <- pairwiseAlignment(string1, string2)
    checkEquals(as.character(aligned(pattern(alignment))), "")
    checkEquals(score(alignment), as.numeric(NA))
}


test_pairwiseAlignment_emptyLocalAlign <- function()
{
    string1 <- DNAString("A")
    string2 <- DNAString("T")
    alignment <- pairwiseAlignment(string1, string2, type = "local")
    checkEquals(as.character(aligned(pattern(alignment))), "")
    checkEquals(as.character(aligned(subject(alignment))), "")
    checkEquals(score(alignment), 0)
}


test_pairwiseAlignment_constantSubstitutionMatrix <- function()
{
    string1 <- DNAString("ACTTCACCAGCTCCCTGGCGGTAAGTTGATCAAAGGAAACGCAAAGTTTTCAAG")
	string2 <- DNAString("GTTTCACTACTTCCTTTCGGGTAAGTAAATATATAAATATATAAAAATATAATTTTCATC")
    mat <- matrix(-3, nrow = 4, ncol = 4)
    diag(mat) <- 1
    rownames(mat) <- colnames(mat) <- DNA_ALPHABET[1:4]
    globalAlign <-
        pairwiseAlignment(string1, string2, substitutionMatrix = mat, gapOpening = -5, gapExtension = -2)
    overlapAlign <-
        pairwiseAlignment(string1, string2, type = "overlap", substitutionMatrix = mat, gapOpening = -5, gapExtension = -2)
    localAlign <-
        pairwiseAlignment(string1, string2, type = "local", substitutionMatrix = mat, gapOpening = -5, gapExtension = -2)
    checkEquals(as.character(pattern(globalAlign)), "ACTTCACCAGCTCCCTGGCGGTAAGTTGATC---AAAGG---AAACGCAAAGTTTTCAAG")
    checkEquals(as.character(subject(globalAlign)), "GTTTCACTACTTCCTTTCGGGTAAGTAAATATATAAATATATAAAAATATAATTTTCATC")
    checkEquals(compareStrings(globalAlign), "??TTCAC?A??TCC?T???GGTAAGT??AT?---AAA??---AAA???A?A?TTTTCA??")
    checkEquals(score(globalAlign), -52)
    checkEquals(as.character(pattern(overlapAlign)), "G")
    checkEquals(as.character(subject(overlapAlign)), "G")
    checkEquals(score(overlapAlign), 1)
    checkEquals(as.character(pattern(localAlign)), "GGTAAGT")
    checkEquals(as.character(subject(localAlign)), "GGTAAGT")
    checkEquals(score(localAlign), 7)
}
