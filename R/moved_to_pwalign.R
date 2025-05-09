### Everything in this file has moved to the pwalign package


.call_fun_in_pwalign <- function(fun, ...)
{
    msg <- c(fun, "() has moved from Biostrings to the pwalign package, ",
             "and is formally defunct in Biostrings >= 2.77.1. ",
             "Please call pwalign::", fun, "() to get rid of this error.")
    .Defunct(msg=wmsg(msg))
}

writePairwiseAlignments <-
    function(...) .call_fun_in_pwalign("writePairwiseAlignments", ...)

nucleotideSubstitutionMatrix <-
    function(...) .call_fun_in_pwalign("nucleotideSubstitutionMatrix", ...)

errorSubstitutionMatrices <-
    function(...) .call_fun_in_pwalign("errorSubstitutionMatrices", ...)

qualitySubstitutionMatrices <-
    function(...) .call_fun_in_pwalign("qualitySubstitutionMatrices", ...)

insertion <-
    function(...) .call_fun_in_pwalign("insertion", ...)

deletion <-
    function(...) .call_fun_in_pwalign("deletion", ...)

unaligned <-
    function(...) .call_fun_in_pwalign("unaligned", ...)

aligned <-
    function(...) .call_fun_in_pwalign("aligned", ...)

indel <-
    function(...) .call_fun_in_pwalign("indel", ...)

nindel <-
    function(...) .call_fun_in_pwalign("nindel", ...)

PairwiseAlignments <-
    function(...) .call_fun_in_pwalign("PairwiseAlignments", ...)

pattern <-
    function(...) .call_fun_in_pwalign("pattern", ...)

alignedPattern <-
    function(...) .call_fun_in_pwalign("alignedPattern", ...)

alignedSubject <-
    function(...) .call_fun_in_pwalign("alignedSubject", ...)

PairwiseAlignmentsSingleSubject <-
    function(...) .call_fun_in_pwalign("PairwiseAlignmentsSingleSubject", ...)

nedit <-
    function(...) .call_fun_in_pwalign("nedit", ...)

mismatchTable <-
    function(...) .call_fun_in_pwalign("mismatchTable", ...)

mismatchSummary <-
    function(...) .call_fun_in_pwalign("mismatchSummary", ...)

compareStrings <-
    function(...) .call_fun_in_pwalign("compareStrings", ...)

pid <-
    function(...) .call_fun_in_pwalign("pid", ...)

pairwiseAlignment <-
    function(...) .call_fun_in_pwalign("pairwiseAlignment", ...)

stringDist <-
    function(...) .call_fun_in_pwalign("stringDist", ...)

print.moved_to_pwalign_pkg <- function(x, ...)
    warning(wmsg("all the BLOSUM* and PAM* scoring matrices are now located ",
                 "in the pwalign package and will soon be removed from the ",
                 "Biostrings package"))

