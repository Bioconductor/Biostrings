### =========================================================================
### readFASTA() / writeFASTA()
### --------------------------
### This is legacy stuff. People should use read[B|DNA|RNA|AA]StringSet() /
### writeXStringSet() which are much faster.
### -------------------------------------------------------------------------


### Robert's contribution
readFASTA <- function(file, checkComments=TRUE, strip.descs=TRUE)
{
    .Defunct("readDNAStringSet")
}

writeFASTA <- function(x, file="", desc=NULL, append=FALSE, width=80)
{
    .Defunct("writeXStringSet")
}

