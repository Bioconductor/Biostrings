## ============================ Utility functions ============================

## Check that the offsets and substrings in BioString object 'matches'
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



## =================== Tests that should FAIL (open bugs) ====================

test_openbug1 <- function()
{
    pattern <- DNAString("CAG")
    subject <- DNAString("CAGTTT")
    expected_roffsets <- rbind(c(1, 3))
    matches <- matchDNAPattern(pattern, subject, mis=1)
    checkBioStringMatches(subject, expected_roffsets, matches)
}

test_openbug2 <- function()
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

test_openbug3 <- function()
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

test_openbug4 <- function()
{
    pattern <- DNAString("AAAA")
    subject <- DNAString("AAAAC")
    ## Here: mis >= nchar(pattern). So you get a match
    ## whatever the relative position of the pattern is!
    ## However, a reasonable expectation is to have
    ## matchDNAPattern only return the following matches:
    expected_roffsets <- rbind(c(-3, 0),
                               c(-2, 1),
                               c(-1, 2),
                               c( 0, 3),
                               c( 1, 4),
                               c( 2, 5),
                               c( 3, 6),
                               c( 4, 7),
                               c( 5, 8),
                               c( 6, 9))
    matches <- matchDNAPattern(pattern, subject, mis=4)
    checkBioStringMatches(subject, expected_roffsets, matches)
}

# TODO: Move this test to another file. Has nothing to do
# with testing the matchDNAPattern() function...
# Not sure this one is a bug. Are we supposed to always
# be allowed to change a value in a matrix? What "says" R?
test_openbug5 <- function()
{
    s <- DNAString("ACGT")
    s@offsets[1,1] <- 2
    ## Provoques an error
    show(s)
}
