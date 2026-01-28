longestConsecutive <- function(seq, letter) {
  .Call2("MP_longestConsecutive", seq, letter, PACKAGE="Biostrings")
}

