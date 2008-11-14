

test_alphabet <- function()
{
  s <- DNAString("ACGT")
  checkEquals(DNA_ALPHABET, alphabet(s))
}



test_subsetting <- function()
{
    s <- DNAString("ACGT")
    checkEquals("AT", as.character(s[c(1,4)]))
    
    checkException(s[1] <- "T", silent=TRUE)
}

