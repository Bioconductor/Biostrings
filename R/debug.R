###

debug_copy_seq <- function()
{
    invisible(.Call("debug_copy_seq", PACKAGE="Biostrings"))
}

debug_utils <- function()
{
    invisible(.Call("Biostrings_debug_utils", PACKAGE="Biostrings"))
}

debug_bufutils <- function()
{
    invisible(.Call("Biostrings_debug_bufutils", PACKAGE="Biostrings"))
}

debug_IRanges_class <- function()
{
    invisible(.Call("debug_IRanges_class", PACKAGE="Biostrings"))
}

debug_IRanges_utils <- function()
{
    invisible(.Call("debug_IRanges_utils", PACKAGE="Biostrings"))
}

debug_XRaw_class <- function()
{
    invisible(.Call("debug_XRaw_class", PACKAGE="Biostrings"))
}

debug_XRaw_fillread <- function()
{
    invisible(.Call("Biostrings_debug_XRaw_fillread", PACKAGE="Biostrings"))
}

debug_XInteger <- function()
{
    invisible(.Call("Biostrings_debug_XInteger", PACKAGE="Biostrings"))
}

debug_XString_class <- function()
{
    invisible(.Call("debug_XString_class", PACKAGE="Biostrings"))
}

debug_XStringSet_class <- function()
{
    invisible(.Call("debug_XStringSet_class", PACKAGE="Biostrings"))
}

debug_seqs_to_seqs <- function()
{
    invisible(.Call("Biostrings_debug_seqs_to_seqs", PACKAGE="Biostrings"))
}

debug_views_buffer <- function()
{
    invisible(.Call("Biostrings_debug_views_buffer", PACKAGE="Biostrings"))
}

debug_naive <- function()
{
    invisible(.Call("match_naive_debug", PACKAGE="Biostrings"))
}

debug_boyermoore <- function()
{
    invisible(.Call("match_boyermoore_debug", PACKAGE="Biostrings"))
}

debug_shiftor <- function()
{
    invisible(.Call("match_shiftor_debug", PACKAGE="Biostrings"))
}

debug_find_palindromes <- function()
{
    invisible(.Call("find_palindromes_debug", PACKAGE="Biostrings"))
}

debug_BOC <- function()
{
    invisible(.Call("match_BOC_debug", PACKAGE="Biostrings"))
}

debug_BOC2 <- function()
{
    invisible(.Call("match_BOC2_debug", PACKAGE="Biostrings"))
}

debug_TBdna <- function()
{
    invisible(.Call("match_TBdna_debug", PACKAGE="Biostrings"))
}

