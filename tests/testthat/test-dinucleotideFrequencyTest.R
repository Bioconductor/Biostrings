## dinucleotideFrequencyTest exports dinucleotideFrequencyTest()
## for the following classes:
##
## - DNAStringSet
## - RNAStringSet

test_that("examples in dinucleotideFrequencyTest.Rd maintain functionality", {
	## TODO: better tests.
	## This is adapted from current functionality, do we know that's correct?
	data(HNF4alpha)
	expect_s3_class(dinucleotideFrequencyTest(HNF4alpha, 1, 2), "htest")

	## why doesn't p.value have a name?
	test1 <- dinucleotideFrequencyTest(HNF4alpha, 1, 2)
	expect_equal(c(round(test1$statistic,4L), test1$parameter, round(test1$p.value, 4L)),
							 c("X-squared"=19.0727, "df"=9, 0.0246))

  test2 <- dinucleotideFrequencyTest(HNF4alpha, 1, 2, test = "G")
  expect_equal(c(round(test2$statistic,4L), test2$parameter, round(test2$p.value,4L)),
  						c("Log likelihood ratio statistic (G)"=17.2609,
  							"X-squared df"=9,
  							"p.value"=0.0448))

  expect_false(grepl("Williams' correction", test2$method))

  test3 <- dinucleotideFrequencyTest(HNF4alpha, 1, 2, test = "adjG")
  expect_equal(c(round(test3$statistic,4L), test3$parameter, round(test3$p.value,4L)),
  						c("Log likelihood ratio statistic (G)"=10.8064,
  							"X-squared df"=9,
  							"p.value"=0.2892))
  expect_true(grepl("Williams' correction", test3$method))
})