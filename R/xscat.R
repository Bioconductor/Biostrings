### =========================================================================
### The xscat() function
### -------------------------------------------------------------------------
###

### Check the types of the arguments and determine the cardinality of the
### answer:
###   if 0 => empty XStringSet
###   if 1 => XString
###   if >= 2 => XStringSet
.get_xscat_ans_cardinality <- function(...)
{
    get_arg_cardinality <- function(arg)
    {
        if (is.character(arg))
            return(length(arg))
        if (is(arg, "XString"))
            return(1L)
        if (is(arg, "XStringSet") || is(arg, "XStringViews")) {
            if (length(arg) == 0L)
                return(0L)
            return(2L)  # yes return 2, even if length(arg) == 1
        }
        stop("xscat() arguments must be character vectors (with no NAs) ",
             "or XString/XStringSet/XStringViews objects")
    }
    arg_cards <- sapply(list(...), get_arg_cardinality)
    if (all(arg_cards == 0L))
        return(0L)
    if (any(arg_cards == 0L))
        stop("xscat() cannot mix arguments made of 0 sequence with arguments ",
             "made of 1 or more sequences")
    if (all(arg_cards == 1L))
        return(1L)
    return(2L)
}

.get_xscat_ans_seqtype <- function(...)
{
    get_arg_seqtype <- function(arg)
    {
        if (is.character(arg)) "B" else seqtype(arg)
    }
    arg_seqtypes <- unique(sapply(list(...), get_arg_seqtype))
    ans_seqtype <- setdiff(arg_seqtypes, "B")
    if (length(ans_seqtype) >= 2L)
        stop("xscat() cannot mix ", ans_seqtype[1L],
             " and ", ans_seqtype[2L], " input")
    if (length(ans_seqtype) == 0L)
        ans_seqtype <- "B"
    ans_seqtype
}

xscat <- function(...)
{
    if (length(list(...)) == 0)
        stop("no input")
    ans_card <- .get_xscat_ans_cardinality(...)
    ans_seqtype <- .get_xscat_ans_seqtype(...)
    if (ans_card == 1L) {
        ans_class <- paste(ans_seqtype, "String", sep="")
    } else {
        ans_class <- paste(ans_seqtype, "StringSet", sep="")
        if (ans_card == 0L)
            return(as(character(0), ans_class))
    }
    args <- lapply(list(...),
                   function(arg)
                   {
                       if (is(arg, ans_class)) arg else as(arg, ans_class)
                   })
    if (ans_card == 1L) {
        .Call2("XString_xscat", args, PACKAGE="Biostrings")
    } else {
        .Call2("XStringSet_xscat", args, PACKAGE="Biostrings")
    }
}

