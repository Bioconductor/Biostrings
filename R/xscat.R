### =========================================================================
### The xscat() function
### -------------------------------------------------------------------------

### All args must be character data i.e. character vectors (with no NAs)
### or XString/XStringSet/XStringViews/MaskedXString objects.
.get_xscat_ans_subtype <- function(...)
{
    getSubtype <- function(arg)
    {
        if (is.character(arg))
            return("BString")
        subtype <- try(xsbaseclass(arg), silent=TRUE)
        if (is(subtype, "try-error"))
            stop("all xscat() arguments must be character data i.e. character ",
                 "vectors or XString/XStringSet/XStringViews/MaskedXString ",
                 "objects")
        subtype
    }
    subtypes <- unique(sapply(list(...), getSubtype))
    ans_subtype <- setdiff(subtypes, "BString")
    if (length(ans_subtype) >= 2)
        stop("xscat() cannot mix ", ans_subtype[1],
             " and ", ans_subtype[2], " input")
    if (length(ans_subtype) == 0)
        ans_subtype <- "BString"
    ans_subtype
}

### Determine the cardinality of the answer:
###   if 0 => empty XStringSet
###   if 1 => XString
###   if >= 2 => XStringSet
.get_xscat_ans_cardinality <- function(...)
{
    getcard <- function(arg)
    {
        if (is.character(arg))
            return(length(arg))
        if (is(arg, "XString") || is(arg, "MaskedXString"))
            return(1)
        ## from here, arg is either XStringSet or XStringViews
        if (length(arg) == 0) 0 else 2  # yes return 2, even if length(arg) == 1
    }
    cards <- sapply(list(...), getcard)
    if (all(cards == 0))
        return(0)
    if (any(cards == 0))
        stop("xscat() cannot mix arguments made of 0 sequence with arguments ",
             "made of 1 or more sequences")
    if (all(cards == 1))
        return(1)
    return(2)
}

### TODO: Calls to as() below will fail because some coercion methods are
### missing! In particular, methods for coercing any subtype of XString
### to any subtype of XStringSet are missing. FIXME: add these methods!
xscat <- function(...)
{
    if (length(list(...)) == 0)
        stop("no input")
    ans_type <- .get_xscat_ans_subtype(...)
    ans_card <- .get_xscat_ans_cardinality(...)
    if (ans_card != 1) {
        ans_type <- paste(ans_type, "Set", sep="")
        if (ans_card == 0)
            return(as(character(0), ans_type))
    }
    args <- lapply(list(...),
                   function(arg)
                       if (is(arg, ans_type)) arg else as(arg, ans_type))
    if (ans_card == 1) {
        .Call("XString_xscat", args, PACKAGE="Biostrings")
    } else {
        .Call("XStringSet_xscat", args, PACKAGE="Biostrings")
    }
}

