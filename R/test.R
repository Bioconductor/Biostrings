## Defined for testing purpose only.

setGeneric(
    "eq",
    function(bs1, bs2) standardGeneric("eq")
)

setMethod(
    "eq",
    signature(bs1="BioString", bs2="BioString"),
    function (bs1, bs2)
    {
        ok <- as.character(bs1) == as.character(bs2)
        ok <- ok && bs1@offsets == bs2@offsets
        return(ok)
    }
)

testBiostrings <- function()
{
    require("RUnit", quietly=TRUE) || stop("RUnit package not found")
    testDirs <- system.file("tests", package="Biostrings")
    testFileRegexp <- "^runit.+\.[rR]$"
    testSuite <- defineTestSuite(name="BioStrings Test Suite",
                                 dirs=testDirs,
                                 testFileRegexp=testFileRegexp)
    testData <- runTestSuite(testSuite)
    printTextProtocol(testData, showDetails=FALSE)
}
