matchprobes <- function(query, records, probepos=FALSE) {
  msg <- "matchprobes() is deprecated. Please use matchPdict() instead."
  .Deprecated(msg=msg)
  .Call2("MP_matchprobes", toupper(query), toupper(records), probepos, PACKAGE="Biostrings")
}

longestConsecutive <- function(seq, letter) {
  msg <- "longestConsecutive() is deprecated and won't be replaced"
  .Deprecated(msg=msg)
  .Call2("MP_longestConsecutive", seq, letter, PACKAGE="Biostrings")
}

