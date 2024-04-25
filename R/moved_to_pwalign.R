### Everything in this file has moved to the pwalign package


### TODO: Move this to S4Vectors (or BiocBaseUtils).
.load_package_gracefully <- function(package, ...)
{
    if (!requireNamespace(package, quietly=TRUE))
        stop("Could not load package ", package, ". Is it installed?\n\n  ",
             wmsg("Note that ", ..., " requires the ", package, " package. ",
                  "Please install it with:"),
             "\n\n    BiocManager::install(\"", package, "\")")
}

.call_fun_in_pwalign <- function(fun, ...)
{
    .load_package_gracefully("pwalign", "starting with BioC 3.19, ",
                             "calling ", fun, "()")
    msg <- c(fun, "() has moved to the pwalign package. Please ",
             "call pwalign::", fun, "() to get rid of this warning.")
    warning(wmsg(msg))
    FUN <- base::get(fun, envir=asNamespace("pwalign"), inherits=FALSE)
    do.call(FUN, list(...))
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

