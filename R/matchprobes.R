matchprobes <- function(query, records, probepos=FALSE) 
   .Call("MP_matchprobes", toupper(query), toupper(records), probepos, PACKAGE="Biostrings")

complementSeq <- function(seq, start=1, stop=0)
{
  msg <- c("'complementSeq' is deprecated.\n",
           "Use 'complement(DNAStringSet(seq))' instead.\n",
           "See '?complementSeq'")
  .Deprecated("complement", msg=msg)
  .Call("MP_complementSeq", seq, as.integer(start), as.integer(stop), PACKAGE="Biostrings")
}

reverseSeq  <- function(seq)
{
  msg <- c("'reverseSeq' is deprecated.\n",
           "Use 'reverse(BStringSet(seq))' instead.\n",
           "See '?reverseSeq'")
  .Deprecated("reverse", msg=msg)
  .Call("MP_revstring", seq, PACKAGE="Biostrings")
}

revcompDNA <- function(seq)
{
  msg <- c("'revcompDNA' is deprecated.\n",
           "Use 'reverseComplement(DNAString(seq))' instead.\n",
           "See '?revcompDNA'")
  .Deprecated("reverseComplement", msg=msg)
  .Call("MP_dna_revcomp", seq, PACKAGE="Biostrings")
}

revcompRNA <- function(seq)
{
  msg <- c("'revcompRNA' is deprecated.\n",
           "Use 'reverseComplement(RNAString(seq))' instead.\n",
           "See '?revcompRNA'")
  .Deprecated("reverseComplement", msg=msg)
  .Call("MP_rna_revcomp", seq, PACKAGE="Biostrings")
}

longestConsecutive <- function(seq, letter) 
  .Call("MP_longestConsecutive", seq, letter, PACKAGE="Biostrings")

basecontent <- function(seq) {
  msg <- c("'basecontent' is deprecated.\n",
           "Use 'alphabetFrequency(DNAStringSet(seq), baseOnly=TRUE)' instead.\n",
           "See '?basecontent'")
  .Deprecated("alphabetFrequency", msg=msg)
  good   = !is.na(seq)
  havena = !all(good)
  if(havena)
    seq = seq[good]
  
  rv = .Call("MP_basecontent", seq, TRUE, PACKAGE="Biostrings")

  if(havena) {
    z = rv
    rv = matrix(NA, nrow=length(good), ncol=ncol(z))
    colnames(rv) = colnames(z)
    rv[good, ] = z
  }
  
  class(rv) <- c("probetable", class(rv))
  return(rv)
}

countbases <- function(seq, dna=TRUE) {
    msg <- c("'basecontent' is deprecated.\n",
             "Use 'alphabetFrequency(DNAStringSet(seq), baseOnly=TRUE)' instead.\n",
             "See '?countbases'")
    .Deprecated("alphabetFrequency", msg=msg)
    good = !is.na(seq)
    havena = !all(good)
    if(havena)
      seq = seq[good]
    rv = .Call("MP_basecontent", seq, dna, PACKAGE="Biostrings")
    if (havena) {
        z = rv
        rv = matrix(NA, nrow=length(good), ncol=ncol(z))
        colnames(rv) = colnames(z)
        rv[good, ] = z
    }
    rv
}
