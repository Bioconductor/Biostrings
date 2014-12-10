testBiostrings <- function()
{
    require("RUnit", quietly=TRUE) || stop("RUnit package not found")
    dirs <- system.file("unitTests", package="Biostrings")
    testFileRegexp <- "^runit.+\\.[rR]$"
    testsuite <- defineTestSuite("Biotrings Test Suite", dirs,
                                 testFileRegexp=testFileRegexp)
    testresult <- runTestSuite(testsuite)
    printTextProtocol(testresult, showDetails=FALSE)
}

