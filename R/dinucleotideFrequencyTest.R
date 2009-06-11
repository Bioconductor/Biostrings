### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "dinucleotideFrequencyTest" generic and methods.
###

.dinucleotideFrequencyTest.XStringSet <-
function(x, i, j, test = c("chisq", "G", "adjG"))
{
    n <- length(x)
    observed <- nucleotideFrequencyAt(x, c(i, j))
    rowTotals <- rowSums(observed)
    colTotals <- colSums(observed)
    expected <- outer(rowTotals, colTotals, FUN = "*")/n
    test <- match.arg(test)
    switch(test,
            chisq = {
                method <-
                  paste("Chi-squared test for independence of positions", i, "and", j)
                ok <- observed != 0 | expected != 0
                statistic <-
                  c("X-squared" = sum((observed[ok] - expected[ok])^2/expected[ok]))
                parameter <- c("df" = (sum(rowTotals > 0) - 1L) * (sum(colTotals > 0) - 1L))
                p.value <- pchisq(statistic, parameter, lower.tail = FALSE)
            },
            G = {
                method <-
                  paste("G-test for independence of positions", i, "and", j)
                ok <- observed != 0
                statistic <-
                  c("G" = 2 * sum(observed[ok] * log(observed[ok]/expected[ok])))
                parameter <- c("df" = (sum(rowTotals > 0) - 1L) * (sum(colTotals > 0) - 1L))
                p.value <- pchisq(statistic, parameter, lower.tail = FALSE)
            },
            adjG = {
                method <-
                  paste("Adjusted G-test for independence of positions", i, "and", j)
                ok <- observed != 0
                a <- (sum(rowTotals > 0) - 1L) * (sum(colTotals > 0) - 1L) - 1L
                denom <- 1 + (a^2 - 1L)/(6 * n * (a - 1L)) 
                statistic <-
                  c("adjG" = 2 * sum(observed[ok] * log(observed[ok]/expected[ok]))/denom)
                parameter <- c("df" = a + 1L)
                p.value <- pchisq(statistic, parameter, lower.tail = FALSE)
            })
    structure(list(statistic = statistic, parameter = parameter,
                   p.value = p.value, method = method,
                   observed = observed, expected = expected),
              class = "htest")
}

setGeneric("dinucleotideFrequencyTest", signature="x",
    function(x, i, j, test = c("chisq", "G", "adjG"))
        standardGeneric("dinucleotideFrequencyTest")
)

setMethod("dinucleotideFrequencyTest", "DNAStringSet",
    .dinucleotideFrequencyTest.XStringSet
)

setMethod("dinucleotideFrequencyTest", "RNAStringSet",
    .dinucleotideFrequencyTest.XStringSet
)
