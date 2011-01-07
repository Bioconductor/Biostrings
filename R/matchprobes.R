matchprobes <- function(query, records, probepos=FALSE) 
  .Call("MP_matchprobes", toupper(query), toupper(records), probepos, PACKAGE="Biostrings")

complementSeq <- function(seq, start=1, stop=0)
{
  msg <- c("'complementSeq' is defunct.\n",
           "Use 'complement(DNAStringSet(seq))' instead.\n",
           "See '?complementSeq'")
  .Defunct("complement", msg=msg)
}

reverseSeq  <- function(seq)
{
  msg <- c("'reverseSeq' is defunct.\n",
           "Use 'reverse(DNAStringSet(seq))' or 'reverse(BStringSet(seq))' instead.\n",
           "See '?reverseSeq'")
  .Defunct("reverse", msg=msg)
}

revcompDNA <- function(seq)
{
  msg <- c("'revcompDNA' is defunct.\n",
           "Use 'reverseComplement(DNAString(seq))' instead.\n",
           "See '?revcompDNA'")
  .Defunct("reverseComplement", msg=msg)
}

revcompRNA <- function(seq)
{
  msg <- c("'revcompRNA' is defunct.\n",
           "Use 'reverseComplement(RNAString(seq))' instead.\n",
           "See '?revcompRNA'")
  .Defunct("reverseComplement", msg=msg)
}

longestConsecutive <- function(seq, letter) 
  .Call("MP_longestConsecutive", seq, letter, PACKAGE="Biostrings")

basecontent <- function(seq) {
  msg <- c("'basecontent' is defunct.\n",
           "Use 'alphabetFrequency(DNAStringSet(seq), baseOnly=TRUE)' instead.\n",
           "See '?basecontent'")
  .Defunct("alphabetFrequency", msg=msg)
}

countbases <- function(seq, dna=TRUE) {
  msg <- c("'countbases' is defunct.\n",
           "Use 'alphabetFrequency(DNAStringSet(seq), baseOnly=TRUE)' instead.\n",
           "See '?countbases'")
  .Defunct("alphabetFrequency", msg=msg)
}
