matchprobes <- function(query, records, probepos=FALSE) 
  .Call2("MP_matchprobes", toupper(query), toupper(records), probepos, PACKAGE="Biostrings")

longestConsecutive <- function(seq, letter) 
  .Call2("MP_longestConsecutive", seq, letter, PACKAGE="Biostrings")

