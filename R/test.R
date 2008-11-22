testBiostrings <- function()
{
    require("RUnit", quietly=TRUE) || stop("RUnit package not found")
    dirs <- file.path(system.file("tests", package="Biostrings"), "UnitTests")
    testFileRegexp <- "^runit.+\\.[rR]$"
    testsuite <- defineTestSuite("Biotrings Test Suite", dirs,
                                 testFileRegexp=testFileRegexp)
    testresult <- runTestSuite(testsuite)
    printTextProtocol(testresult, showDetails=FALSE)
}

