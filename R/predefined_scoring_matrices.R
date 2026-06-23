
print.moved_to_pwalign_pkg <- function(x, ...)
{
    msg <- c("All the BLOSUM* and PAM* scoring matrices are now located ",
             "in the pwalign package. Starting with BioC 3.24, the scoring ",
             "matrices located in the Biostrings package are formally ",
             "deprecated and will soon be removed from the package. ",
             "Please use the scoring matrices located in the pwalign ",
             "package instead.")
    .Deprecated(msg=wmsg(msg))
}

