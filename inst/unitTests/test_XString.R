test_constructors <- function()
{

  dst <- DNAString("ACGTAT")   # just check this does not throw an error
  checkException( dst <- DNAString("ACFTC") , silent=TRUE)

  rst <- RNAString("ACGUAU")   # just check this does not throw an error
  checkException( rst <- RNAString("ACGTC") , silent=TRUE)

  ast <- AAString("EQHIILF")   # just check this does not throw an error
  #checkException( ast <- AAString("ACXTCZ") , silent=TRUE)  # AAString()
                                            # doesn't check its input yet

}

test_alphabet <- function()
{
  s <- BString("ACGTAT")
  checkEquals(NULL, alphabet(s))

  s <- DNAString("ACGTAT")
  checkEquals(DNA_ALPHABET, alphabet(s))
}

test_length <- function()
{
  rs <- "ACGTAT"
  bs <- BString("ACGTAT")
  checkEquals(nchar(rs), length(bs))

  checkEquals(nchar(rs), nchar(bs))
}

test_subsetting <- function()
{
    s <- BString("ACGT")

    ## out of bounds
    checkException(s[10], silent=TRUE)

    checkEquals("AT", as.character(s[c(1,4)]))
    checkEquals("CGT", as.character(s[-1]))

    ## replacement
    s[1] <- "T"   
    checkEquals("TCGT", as.character(s))
}

test_comparison <- function()
{
  ref_str <- "ACGTA"
  bs1_acgta <- BString(ref_str)
  bs2_acgta <- BString(ref_str)
  bs_aatt <- BString("AATT")
  
  checkTrue(bs1_acgta == bs2_acgta)
  checkTrue(ref_str == bs1_acgta)
  checkTrue(bs1_acgta != bs_aatt)
  
  ds1_acgta <- DNAString(ref_str)
  ds2_acgta <- DNAString(ref_str)
  checkTrue(ds1_acgta == ds2_acgta)
  
  ## comparison between instances of different XString subclasses
  ## not possible
  checkException(bs1_acgta == ds1_acgta, silent=TRUE)
  
}

test_subseq <- function()
{
  bs <- BString("ACGTA")
  checkEquals("GTA", as.character(subseq(bs, start = 3))) # test "==" between "character"
                                                          # and "BString" done elsewhere
  checkEquals("TA", as.character(subseq(bs, start = -2)))
  
}
