## ============================ Utility functions ============================

## Check that the offsets and substrings in Biostring object 'matches'
## are the expected ones.
##   subject: BioString object of length 1
##   expected_roffsets: matrix of the expected relative offsets of the matches
##   matches: BioString object returned by matchDNAPattern("findme", subject)
checkBioStringMatches <- function(subject, expected_roffsets, matches)
{
    expected_aoffsets <- expected_roffsets + subject@offsets[1,1] - 1
    checkEquals(expected_aoffsets, matches@offsets)
    ## We need to check for this particular case because surprinsingly,
    ## substring() doesn't work with 'first' and 'last' of zero length
    if (dim(expected_roffsets)[1] != 0) {
        first <- expected_roffsets[,1]
        last <- expected_roffsets[,2]
        expected_substrings <- substring(as.character(subject), first, last)
        checkEquals(expected_substrings, as.character(matches))
    }
}

checkExactMatches <- function(pattern, subject, expected_pos=numeric(0))
{
    matches <- matchDNAPattern(pattern, subject)
    expected_roffsets <- cbind(expected_pos, expected_pos+nchar(pattern)-1, deparse.level=0)
    checkBioStringMatches(subject, expected_roffsets, matches)
    expected_length <- length(expected_pos)
    checkEquals(logical(expected_length), as.character(matches) != pattern)
}



## =================== Tests that should PASS (fixed bugs) ===================

## Checks that the fix in the matchIndexToBioString() C-function
## for the "wrong offset" bug works for every cases this function
## deals with.
test_matchIndexToBioString <- function()
{
    seq <- c("CCCC", "TTAATT")
    subject <- DNAString(seq)[2]

    ## CASE 'nmatch == 0'
    checkExactMatches("AAA", subject)
    ## CASE 'nmatch == 1'
    checkExactMatches("AA", subject, c(3))
    ## CASE 'nmatchIndex == 2*nmatch'
    checkExactMatches("TT", subject, c(1,5))
    ## DEFAULT CASE
    checkExactMatches("A", subject, c(3,4))
}

## Checks the ShiftOr_matchInternal() C-function.
test_ShiftOr_matchInternal_A1 <- function()
{
    pattern <- DNAString("CAG")
    subject <- DNAString("CAGTTT")
    expected_roffsets <- rbind(c(1, 3))
    matches <- matchDNAPattern(pattern, subject, mis=1)
    checkBioStringMatches(subject, expected_roffsets, matches)
}

test_ShiftOr_matchInternal_A2 <- function()
{
    pattern <- DNAString("AKADAKA")
    subject <- DNAString("ATATGAATAAAGA")
    expected_roffsets <- rbind(c(7, 13))
    matches <- matchDNAPattern(pattern, subject, mis=0)
    checkBioStringMatches(subject, expected_roffsets, matches)
}

test_ShiftOr_matchInternal_A3 <- function()
{
    pattern <- DNAString("AKADAKA")
    subject <- DNAString("ATATGAATAAAGA")
    expected_roffsets <- rbind(c(3, 9),
                               c(7, 13))
    matches <- matchDNAPattern(pattern, subject, mis=1)
    checkBioStringMatches(subject, expected_roffsets, matches)
}

test_ShiftOr_matchInternal_A4 <- function()
{
    pattern <- DNAString("AKADAKA")
    subject <- DNAString("ATATGAATAAAGA")
    expected_roffsets <- rbind(c(1,7),
                               c(3, 9),
                               c(7, 13))
    matches <- matchDNAPattern(pattern, subject, mis=2)
    checkBioStringMatches(subject, expected_roffsets, matches)
}

test_ShiftOr_matchInternal_B1 <- function()
{
    pattern <- DNAString("AAAA")
    subject <- DNAString("AAAAC")
    expected_roffsets <- rbind(c(-1, 2),
                               c( 0, 3),
                               c( 1, 4),
                               c( 2, 5),
                               c( 3, 6))
    matches <- matchDNAPattern(pattern, subject, mis=2)
    checkBioStringMatches(subject, expected_roffsets, matches)
}

test_ShiftOr_matchInternal_B2 <- function()
{
    pattern <- DNAString("ATC")
    subject <- DNAString("CATCACTCA")
    expected_roffsets <- rbind(c(-1, 1),
                               c( 2, 4),
                               c( 4, 6),
                               c( 5, 7),
                               c( 6, 8),
                               c( 9, 11))
    matches <- matchDNAPattern(pattern, subject, mis=2)
    checkBioStringMatches(subject, expected_roffsets, matches)
}

test_ShiftOr_matchInternal_B3 <- function()
{
    pattern <- DNAString("AAAA")
    subject <- DNAString("AAAAC")
    expected_roffsets <- rbind(c(-2, 1),
                               c(-1, 2),
                               c( 0, 3),
                               c( 1, 4),
                               c( 2, 5),
                               c( 3, 6),
                               c( 4, 7))
    ## Provoques an error
    matches <- matchDNAPattern(pattern, subject, mis=3)
    checkBioStringMatches(subject, expected_roffsets, matches)
}

test_ShiftOr_matchInternal_B4 <- function()
{
    pattern <- DNAString("AAAA")
    subject <- DNAString("AAAAC")
    ## Here: mis >= nchar(pattern). So you get a match
    ## whatever the relative position of the pattern is!
    ## However, a reasonable expectation is to have
    ## matchDNAPattern only return the following matches:
    expected_roffsets <- rbind(c(-2, 1),
                               c(-1, 2),
                               c( 0, 3),
                               c( 1, 4),
                               c( 2, 5),
                               c( 3, 6),
                               c( 4, 7),
                               c( 5, 8))
    matches <- matchDNAPattern(pattern, subject, mis=4)
    checkBioStringMatches(subject, expected_roffsets, matches)
}



## ========== Tests that should PASS on 64-bit platform (fixed bugs) =========
## ================= and FAIL on 32-bit platform (not a bug) =================

## This test uses a big subject (10 millions of characters)
## in order to test the speed of the ShiftOr algorithm.
test_ShiftOr_matchInternal_C1 <- function()
{
    f <- file(system.file("extdata", "bigrandomTGCA.txt", package="Biostrings"))
    big <- scan(file=f, what="")
    subject <- DNAString(big)
    
    ## Simple pattern, mis=5
    pattern <- DNAString("TTTTTTTTTTTTTTTTTTTTT")
    expected_roffsets <- rbind(c(2725443, 2725463),
                               c(6535062, 6535082),
                               c(6535064, 6535084),
                               c(7765179, 7765199),
                               c(8491897, 8491917),
                               c(8491898, 8491918),
                               c(8491899, 8491919),
                               c(9233437, 9233457))
    matches <- matchDNAPattern(pattern, subject, mis=5)
    checkBioStringMatches(subject, expected_roffsets, matches)

    ## With a 32 char long pattern, mis=10
    pattern32 <- DNAString(substr(big, 1, 32)) #substr(subject, 1, 32)
    expected_roffsets <- rbind(c(1,32),
                               c(3807182, 3807213),
                               c(9926155, 9926186))
    matches <- matchDNAPattern(pattern32, subject, mis=10)
    checkBioStringMatches(subject, expected_roffsets, matches)

    ## With a 64 char long pattern (works on 64-bit platforms only)
    pattern64 <- DNAString(substr(big, 1, 64)) #substr(subject, 1, 64)
    expected_roffsets <- rbind(c(1,64),
                               c(1049302, 1049365))
    matches <- matchDNAPattern(pattern64, subject, mis=24)
    checkBioStringMatches(subject, expected_roffsets, matches)
}


## =================== Tests that should FAIL (open bugs) ====================
## ===========================================================================

test_openbug1 <- function()
{
    pattern <- DNAString("CAG")
    subject <- DNAString("GTTCA")
    matches <- matchDNAPattern(pattern, subject, mis=2)

    multi <- DNAString(c("TTT", "GTTCA", "TT"))
    subject2 <- multi[2]
    matches2 <- matchDNAPattern(pattern, subject2, mis=2)

    checkEquals(as.character(matches), as.character(matches2))
}
   
# TODO: Move this test to another file. Has nothing to do
# with testing the matchDNAPattern() function...
# Not sure this one is a bug. Are we supposed to always
# be allowed to change a value in a matrix? What "says" R?
test_openbug2 <- function()
{
    s <- DNAString("ACGT")
    s@offsets[1,1] <- 2
    ## Provoques an error
    show(s)
}
