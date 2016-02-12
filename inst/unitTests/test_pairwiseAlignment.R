### FIXME!!
BROKEN_test_pairwiseAlignment_emptyString <- function()
{
    string1 <- DNAStringSet("")
    string2 <- DNAStringSet("ACGT")

    ## Empty pattern.
    alignment <- pairwiseAlignment(string1, string2)
    checkEquals(as.character(aligned(pattern(alignment))), "")
    checkEquals(as.character(aligned(subject(alignment))), "")
    checkEquals(score(alignment), -26)
    checkEquals(pairwiseAlignment(string1, string2, scoreOnly = TRUE), -26)

    ## Empty subject.
    alignment <- pairwiseAlignment(string2, string1)
    checkEquals(as.character(aligned(pattern(alignment))), "")
    checkEquals(as.character(aligned(subject(alignment))), "")
    checkEquals(score(alignment), -26)
    checkEquals(pairwiseAlignment(string2, string1, scoreOnly = TRUE), -26)

    ## Empty pattern and subject.
    alignment <- pairwiseAlignment(string1, string1)
    checkEquals(as.character(aligned(pattern(alignment))), "")
    checkEquals(score(alignment), 0)
    checkEquals(pairwiseAlignment(string1, string1, scoreOnly = TRUE), 0)
}


test_pairwiseAlignment_emptyLocalAlign <- function()
{
    string1 <- DNAString("A")
    string2 <- DNAString("T")
    alignment <- pairwiseAlignment(string1, string2, type = "local")
    checkEquals(as.character(aligned(pattern(alignment))), "")
    checkEquals(as.character(aligned(subject(alignment))), "")
    checkEquals(score(alignment), 0)
    checkEquals(pairwiseAlignment(string1, string2, type = "local", scoreOnly = TRUE), 0)
}


test_pairwiseAlignment_gappedLocalAlign <- function()
{
    string1 <- DNAString("TCAGTTGCCAAACCCGCT")
    string2 <- DNAString("AGGGTTGACATCCGTTTT")
    sigma <- nucleotideSubstitutionMatrix(match = 10, mismatch = -10, baseOnly = TRUE)
    alignment <-
      pairwiseAlignment(string1, string2, substitutionMatrix = sigma,
                        gapOpening = 12,  gapExtension = 3, type="local")
    checkEquals(as.character(aligned(pattern(alignment))), "GTTGCCAAACCCG")
    checkEquals(as.character(aligned(subject(alignment))), "GTTGACAT--CCG")
    checkEquals(score(alignment), 52)
}


test_pairwiseAlignment_backToBackIndel <- function()
{
    mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -10, baseOnly = TRUE)
    string1 <- DNAString("AC")
    string2 <- DNAString("AT")
    alignment <- pairwiseAlignment(string1, string2, gapOpening = 0, substitutionMatrix = mat)
    alignmentScore <- pairwiseAlignment(string1, string2, gapOpening = 0, substitutionMatrix = mat, scoreOnly = TRUE)
    checkEquals(as.character(aligned(pattern(alignment))), "A-")
    checkEquals(as.character(aligned(subject(alignment))), "AT")
    checkEquals(score(alignment), -7)
    checkEquals(alignmentScore, -7)
}


test_pairwiseAlignment_editDistance <- function()
{
    string1 <- DNAString("ACTTCACCAGCTCCCTGGCGGTAAGTTGATCAAAGGAAACGCAAAGTTTTCAAG")
    string2 <- DNAString("GTTTCACTACTTCCTTTCGGGTAAGTAAATATATAAATATATAAAAATATAATTTTCATC")
    mat <- nucleotideSubstitutionMatrix(match = 0, mismatch = -1, baseOnly = TRUE)
    globalAlign <-
        pairwiseAlignment(string1, string2, substitutionMatrix = mat, gapOpening = 0, gapExtension = 1)
    globalAlignScore <-
        pairwiseAlignment(string1, string2, substitutionMatrix = mat, gapOpening = 0, gapExtension = 1, scoreOnly = TRUE)
    checkEquals(as.character(pattern(globalAlign)), "ACTTCACCAGCTCCCTGGCGG-TAAGTTGATC-A-AAGGA-A-ACGCA-A-AGTTTTCAAG")
    checkEquals(as.character(subject(globalAlign)), "GTTTCACTA-CTTCCTTTCGGGTAAGTAAATATATAAATATATAAAAATATAATTTTCATC")
    checkEquals(compareStrings(globalAlign), "??TTCAC?A+CT?CCT??CGG-TAAGT??AT?-A-AA??A-A-A???A-A-A?TTTTCA??")
    checkEquals(score(globalAlign), -25)
    checkEquals(globalAlignScore, -25)
}


test_pairwiseAlignment_zeroOpening <- function()
{
    string1 <- DNAString("ACTTCACCAGCTCCCTGGCGGTAAGTTGATCAAAGGAAACGCAAAGTTTTCAAG")
    string2 <- DNAString("GTTTCACTACTTCCTTTCGGGTAAGTAAATATATAAATATATAAAAATATAATTTTCATC")
    mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
    globalAlign <-
        pairwiseAlignment(string1, string2, substitutionMatrix = mat, gapOpening = 0, gapExtension = 5)
    globalAlignScore <-
        pairwiseAlignment(string1, string2, substitutionMatrix = mat, gapOpening = 0, gapExtension = 5, scoreOnly = TRUE)
    overlapAlign <-
        pairwiseAlignment(string1, string2, type = "overlap", substitutionMatrix = mat, gapOpening = 0, gapExtension = 5)
    overlapAlignScore <-
        pairwiseAlignment(string1, string2, type = "overlap", substitutionMatrix = mat, gapOpening = 0, gapExtension = 5,
                          scoreOnly = TRUE)
    localAlign <-
        pairwiseAlignment(string1, string2, type = "local", substitutionMatrix = mat, gapOpening = 0, gapExtension = 5)
    localAlignScore <-
        pairwiseAlignment(string1, string2, type = "local", substitutionMatrix = mat, gapOpening = 0, gapExtension = 5,
                          scoreOnly = TRUE)
    checkEquals(as.character(pattern(globalAlign)), "ACTTCACCAGCTCCCTGGCGG-TAAGTTGATC-A-AAGGA-A-ACGCA-A-AGTTTTCAAG")
    checkEquals(as.character(subject(globalAlign)), "GTTTCACTA-CTTCCTTTCGGGTAAGTAAATATATAAATATATAAAAATATAATTTTCATC")
    checkEquals(compareStrings(globalAlign), "??TTCAC?A+CT?CCT??CGG-TAAGT??AT?-A-AA??A-A-A???A-A-A?TTTTCA??")
    checkEquals(score(globalAlign), -55)
    checkEquals(globalAlignScore, -55)
    checkEquals(as.character(pattern(overlapAlign)), "G")
    checkEquals(as.character(subject(overlapAlign)), "G")
    checkEquals(score(overlapAlign), 1)
    checkEquals(overlapAlignScore, 1)
    checkEquals(as.character(pattern(localAlign)), "GGTAAGT")
    checkEquals(as.character(subject(localAlign)), "GGTAAGT")
    checkEquals(score(localAlign), 7)
    checkEquals(localAlignScore, 7)
}


test_pairwiseAlignment_fixedSubstitutionMatrix <- function()
{
    string1 <- DNAString("ACTTCACCAGCTCCCTGGCGGTAAGTTGATCAAAGGAAACGCAAAGTTTTCAAG")
    string2 <- DNAString("GTTTCACTACTTCCTTTCGGGTAAGTAAATATATAAATATATAAAAATATAATTTTCATC")
    mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
    globalAlign <-
        pairwiseAlignment(string1, string2, substitutionMatrix = mat, gapOpening = 5, gapExtension = 2)
    globalAlignScore <-
        pairwiseAlignment(string1, string2, substitutionMatrix = mat, gapOpening = 5, gapExtension = 2, scoreOnly = TRUE)
    overlapAlign <-
        pairwiseAlignment(string1, string2, type = "overlap", substitutionMatrix = mat, gapOpening = 5, gapExtension = 2)
    overlapAlignScore <-
        pairwiseAlignment(string1, string2, type = "overlap", substitutionMatrix = mat, gapOpening = 5, gapExtension = 2,
                          scoreOnly = TRUE)
    localAlign <-
        pairwiseAlignment(string1, string2, type = "local", substitutionMatrix = mat, gapOpening = 5, gapExtension = 2)
    localAlignScore <-
        pairwiseAlignment(string1, string2, type = "local", substitutionMatrix = mat, gapOpening = 5, gapExtension = 2,
                          scoreOnly = TRUE)
    checkEquals(as.character(pattern(globalAlign)), "ACTTCACCAGCTCCCTGGCGGTAAGTTGATC---AAAGG---AAACGCAAAGTTTTCAAG")
    checkEquals(as.character(subject(globalAlign)), "GTTTCACTACTTCCTTTCGGGTAAGTAAATATATAAATATATAAAAATATAATTTTCATC")
    checkEquals(compareStrings(globalAlign), "??TTCAC?A??TCC?T???GGTAAGT??AT?---AAA??---AAA???A?A?TTTTCA??")
    checkEquals(score(globalAlign), -52)
    checkEquals(globalAlignScore, -52)
    checkEquals(as.character(pattern(overlapAlign)), "G")
    checkEquals(as.character(subject(overlapAlign)), "G")
    checkEquals(score(overlapAlign), 1)
    checkEquals(overlapAlignScore, 1)
    checkEquals(as.character(pattern(localAlign)), "GGTAAGT")
    checkEquals(as.character(subject(localAlign)), "GGTAAGT")
    checkEquals(score(localAlign), 7)
    checkEquals(localAlignScore, 7)
}


test_pairwiseAlignment_qualityScoring <- function()
{
    string1 <- DNAString("ACTTCACCAGCTCCCTGGCGGTAAGTTGATCAAAGGAAACGCAAAGTTTTCAAG")
    string2 <- DNAString("GTTTCACTACTTCCTTTCGGGTAAGTAAATATATAAATATATAAAAATATAATTTTCATC")
    classes <- c("PhredQuality", "SolexaQuality", "IlluminaQuality")
    for (qualityClass in classes) {
        scoring <- qualitySubstitutionMatrices(qualityClass = qualityClass)["22", "22", c("1", "0")]
        stringQuality <- do.call(qualityClass, list(22L))
        globalAlign <-
          pairwiseAlignment(string1, string2, patternQuality = stringQuality, subjectQuality = stringQuality)
        globalAlignScore <-
          pairwiseAlignment(string1, string2, scoreOnly = TRUE, patternQuality = stringQuality, subjectQuality = stringQuality)
        overlapAlign <-
          pairwiseAlignment(string1, string2, type = "overlap", patternQuality = stringQuality, subjectQuality = stringQuality)
        overlapAlignScore <-
          pairwiseAlignment(string1, string2, type = "overlap", scoreOnly = TRUE, patternQuality = stringQuality, subjectQuality = stringQuality)
        localAlign <-
          pairwiseAlignment(string1, string2, type = "local", patternQuality = stringQuality, subjectQuality = stringQuality)
        localAlignScore <-
          pairwiseAlignment(string1, string2, type = "local", scoreOnly = TRUE, patternQuality = stringQuality, subjectQuality = stringQuality)
        checkEquals(as.character(pattern(globalAlign)), "ACTTCACCAGCTCCCTGGCGGTAAGTTGATC---AAAGG---AAACGCAAAGTTTTCAAG")
        checkEquals(as.character(subject(globalAlign)), "GTTTCACTACTTCCTTTCGGGTAAGTAAATATATAAATATATAAAAATATAATTTTCATC")
        checkEquals(compareStrings(globalAlign), "??TTCAC?A??TCC?T???GGTAAGT??AT?---AAA??---AAA???A?A?TTTTCA??")
        checkEquals(score(globalAlign), sum(c(33, 21) * scoring) - 44, tolerance = 1e-6)
        checkEquals(globalAlignScore, sum(c(33, 21) * scoring) - 44, tolerance = 1e-6)
        checkEquals(as.character(pattern(overlapAlign)), "G")
        checkEquals(as.character(subject(overlapAlign)), "G")
        checkEquals(score(overlapAlign), scoring[[1]], tolerance = 1e-6)
        checkEquals(overlapAlignScore, scoring[[1]], tolerance = 1e-6)
        checkEquals(as.character(pattern(localAlign)), "GGTAAGT")
        checkEquals(as.character(subject(localAlign)), "GGTAAGT")
        checkEquals(score(localAlign), 7 * scoring[[1]], tolerance = 1e-6)
        checkEquals(localAlignScore, 7 * scoring[[1]], tolerance = 1e-6)
    }
    TRUE
}
