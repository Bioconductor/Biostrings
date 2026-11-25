
print.moved_to_pwalign_pkg <- function(x, ...)
{
    msg <- c("All the BLOSUM* and PAM* scoring matrices are now located ",
             "in the pwalign package. Starting with BioC 3.25, the scoring ",
             "matrices located in the Biostrings package are formally ",
             "defunct and will soon be removed from the package. ",
             "Please use the scoring matrices located in the pwalign ",
             "package instead.")
    .Defunct(msg=wmsg(msg))
}

