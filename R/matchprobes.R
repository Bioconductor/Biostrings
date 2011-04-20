matchprobes <- function(query, records, probepos=FALSE) 
  .Call("MP_matchprobes", toupper(query), toupper(records), probepos, PACKAGE="Biostrings")

longestConsecutive <- function(seq, letter) 
  .Call("MP_longestConsecutive", seq, letter, PACKAGE="Biostrings")

