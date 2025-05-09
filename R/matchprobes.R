matchprobes <- function(query, records, probepos=FALSE) {
  msg <- "matchprobes() is defunct. Please use matchPdict() instead."
  .Defunct(msg=msg)
}

longestConsecutive <- function(seq, letter) {
  .Call2("MP_longestConsecutive", seq, letter, PACKAGE="Biostrings")
}

