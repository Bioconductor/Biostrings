### =========================================================================
### Environment for storing run-time objects
###

RTobjs <- new.env(hash=TRUE, parent=emptyenv())


### =========================================================================
### Serialized objects
###

SERIALIZED_OBJNAMES <- c(
    "BLOSUM45",
    "BLOSUM50",
    "BLOSUM62",
    "BLOSUM80",
    "BLOSUM100",
    "PAM30",
    "PAM40",
    "PAM70",
    "PAM120",
    "PAM250"
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Objects created "on-the-fly" (not serialized)
###
### WARNING: Improper calls to 'getdata' by the 'createObject' function can
### lead to infinite recursive loops!
###

createObject <- function(objname)
{
    # add more here...
    stop("don't know how to create object '", objname, "'")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "getdata" function (NOT exported)
###

datacache <- new.env(hash=TRUE, parent=emptyenv())

getdata <- function(objname)
{
    if (!exists(objname, envir=datacache)) {
        if (objname %in% SERIALIZED_OBJNAMES) {
            data(list=objname, package="Biostrings", envir=datacache)
        } else {
            assign(objname, createObject(objname), envir=datacache)
        }
    }
    get(objname, envir=datacache)
}

