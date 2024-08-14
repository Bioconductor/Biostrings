## matchPWM.R exports the following:
## - PWM
## - matchPWM
## - countPWM
## - PWMscoreStartingAt
## - maxWeights
## - minWeights
## - minScore
## - unitScale
## - reverseComplement

## TODO: more tests, this is really just a first pass
test_that("PWMs can be initialized", {
	d <- DNAStringSet(c("AAAATTC",
											"ATTAAAG"))

	p <- PWM(d, type='prob')
	expect_equal(dim(p), c(4L, width(d)[1]))
	## have to use all.equal because of floating point error
	expect_true(all.equal(colSums(p), rep(sum(p[,1L]), ncol(p))))
	expect_true(all(p[c("C", "G"),-7] == 0))
	expect_true(maxScore(p) == 1)
	expect_true(minScore(p) == 0)
})