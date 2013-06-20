\name{replaceAt}

\alias{extractAt}
\alias{replaceAt}

\alias{extractAt,XString-method}
\alias{extractAt,XStringSet-method}
\alias{replaceAt,XString-method}
\alias{replaceAt,XStringSet-method}


\title{Extract/replace arbitrary substrings from/in a string or set
       of strings.}

\description{
  \code{extractAt} extracts multiple subsequences from \link{XString}
  object \code{x}, or from the individual sequences of \link{XStringSet}
  object \code{x}, at the ranges of positions specified thru \code{at}.

  \code{replaceAt} performs multiple subsequence replacements (a.k.a.
  substitutions) in \link{XString} object \code{x}, or in the individual
  sequences of \link{XStringSet} object \code{x}, at the ranges of positions
  specified thru \code{at}.
}

\usage{
extractAt(x, at)
replaceAt(x, at, value="")
}

\arguments{
  \item{x}{
    An \link{XString} or \link{XStringSet} object.
  }
  \item{at}{
    Typically a \link[IRanges]{Ranges} object if \code{x} is an
    \link{XString} object, and a \link[IRanges]{RangesList} object
    if \code{x} is an \link{XStringSet} object.

    Alternatively, the ranges can be specified with only 1 number
    per range (its start position), in which case they are considered
    to be empty ranges (a.k.a. zero-width ranges). So if \code{at}
    is a numeric vector, an \link[IRanges]{IntegerList} object,
    or a list of numeric vectors, each number in it is interpreted
    as the start position of a zero-width range.
    This is useful when using \code{replaceAt} to perform insertions.

    The following applies only if \code{x} is an \link{XStringSet} object:

    \code{at} is recycled to the length of \code{x} if necessary.
    If \code{at} is a \link[IRanges]{Ranges} object (or a numeric vector),
    it is first turned into a \link[IRanges]{RangesList} object of length 1
    and then this \link[IRanges]{RangesList} object is recycled to the
    length of \code{x}.
    This is useful for specifying the same ranges across all sequences
    in \code{x}.
    The \emph{effective shape} of \code{at} is described by its length
    together with the lengths of its list elements \emph{after} recycling.

    As a special case, \code{extractAt} accepts \code{at} and \code{value}
    to be both of length 0, in which case it just returns \code{x} unmodified
    (no-op).
  }
  \item{value}{
    The replacement sequences.

    If \code{x} is an \link{XString} object, \code{value} is typically
    a character vector or an \link{XStringSet} object that is recycled
    to the length of \code{at} (if necessary).

    If \code{x} is an \link{XStringSet} object, \code{value} is typically
    a list of character vectors or a \link[IRanges]{CharacterList} or
    \link{XStringSetList} object. If necessary, it is recycled 
    "vertically" first and then "horizontally" to bring it into the
    \emph{effective shape} of \code{at} (see above). "Vertical recycling"
    is the usual recycling whereas "horizontal recycling" recycles the
    individual list elements .

    As a special case, \code{extractAt} accepts \code{at} and \code{value}
    to be both of length 0, in which case it just returns \code{x} unmodified
    (no-op).
  }
}

\value{
  For \code{extractAt}: An \link{XStringSet} object of the same length as
  \code{at} if \code{x} is an \link{XString} object.
  An \link{XStringSetList} object of the same length as \code{x}
  (and same \emph{effective shape} as \code{at}) if \code{x} is an
  \link{XStringSet} object.

  For \code{replaceAt}: An object of the same class as \code{x}.
  If \code{x} is an \link{XStringSet} object, its length and names and
  metadata columns are preserved.
}

\author{H. Pages}

\seealso{
  \itemize{
    \item The \code{\link[XVector]{subseq}} and \code{\link[XVector]{subseq<-}}
          functions in the \pkg{XVector} package for simpler forms of
          subsequence extractions or replacements.

    \item The \code{\link{replaceLetterAt}} function for a DNA-specific
          single-letter replacement functions useful for SNP injections.

    \item The \code{\link{padAndClip}} function for padding and clipping
          strings.

    \item The \link{XString}, \link{XStringSet}, and \link{XStringSetList}
          classes.

    \item The \link[IRanges]{Ranges}, \link[IRanges]{RangesList},
          \link[IRanges]{IntegerList}, and \link[IRanges]{CharacterList}
          classes in the \pkg{IRanges} package.
  }
}

\examples{
## ---------------------------------------------------------------------
## (A) ON AN XString OBJECT
## ---------------------------------------------------------------------
x <- BString("abcdefghijklm")

at1 <- IRanges(5:1, width=3)
extractAt(x, at1)
names(at1) <- LETTERS[22:26]
extractAt(x, at1)

at2 <- IRanges(c(1, 5, 12), c(3, 4, 12), names=c("X", "Y", "Z"))
extractAt(x, at2)
extractAt(x, rev(at2))

value <- c("+", "-", "*")
replaceAt(x, at2, value=value)
replaceAt(x, rev(at2), value=rev(value))

at3 <- IRanges(c(14, 1, 1, 1, 1, 11), c(13, 0, 10, 0, 0, 10))
value <- 1:6
replaceAt(x, at3, value=value)            # "24536klm1"
replaceAt(x, rev(at3), value=rev(value))  # "54236klm1"

## Deletions:
stopifnot(replaceAt(x, at2) == "defghijkm")
stopifnot(replaceAt(x, rev(at2)) == "defghijkm")
stopifnot(replaceAt(x, at3) == "klm")
stopifnot(replaceAt(x, rev(at3)) == "klm")

## Insertions:
at4 <- IRanges(c(6, 10, 2, 5), width=0)
stopifnot(replaceAt(x, at4, value="-") == "a-bcd-e-fghi-jklm")
stopifnot(replaceAt(x, start(at4), value="-") == "a-bcd-e-fghi-jklm")
at5 <- c(5, 1, 6, 5)  # 2 insertions before position 5 
replaceAt(x, at5, value=c("+", "-", "*", "/"))

## No-ops:
stopifnot(replaceAt(x, NULL, value=NULL) == x)
stopifnot(replaceAt(x, at2, value=extractAt(x, at2)) == x)
stopifnot(replaceAt(x, at3, value=extractAt(x, at3)) == x)
stopifnot(replaceAt(x, at4, value=extractAt(x, at4)) == x)
stopifnot(replaceAt(x, at5, value=extractAt(x, at5)) == x)

## The order of successive transformations matters:
##   T1: insert "+" before position 1 and 4
##   T2: insert "-" before position 3

## T1 followed by T2
x2a <- replaceAt(x, c(1, 4), value="+")
x3a <- replaceAt(x2a, 3, value="-")

## T2 followed by T1
x2b <- replaceAt(x, 3, value="-")
x3b <- replaceAt(x2b, c(1, 4), value="+")

## T1 and T2 simultaneously:
x3c <- replaceAt(x, c(1, 3, 4), value=c("+", "-", "+"))

## ==> 'x3a', 'x3b', and 'x3c' are all different!

## Append "**" to 'x3c':
replaceAt(x3c, length(x3c) + 1L, value="**")

## ---------------------------------------------------------------------
## (B) ON AN XStringSet OBJECT
## ---------------------------------------------------------------------
x <- BStringSet(c(seq1="ABCD", seq2="abcdefghijk", seq3="XYZ"))

at5 <- IRangesList(IRanges(c(2, 1), c(3, 0)),
                   IRanges(c(7, 2, 12, 7), c(6, 5, 11, 8)),
                   IRanges(2, 2))
## Set inner names on 'at5'.
unlisted_at5 <- unlist(at5)
names(unlisted_at5) <-
    paste0("rg", sprintf("\%02d", seq_along(unlisted_at5)))
at5 <- relist(unlisted_at5, at5)

extractAt(x, at5)  # same as 'as(mapply(extractAt, x, at5), "List")'
extractAt(x, at5[3])  # same as 'as(mapply(extractAt, x, at5[3]), "List")'

replaceAt(x, at5, value=extractAt(x, at5))  # no-op
replaceAt(x, at5)  # deletions

at6 <- IRangesList(IRanges(1:5, width=0),
                   IRanges(c(6, 8, 10, 7, 2, 5),
                           width=c(0, 2, 0, 0, 0, 0)),
                   IRanges(c(1, 2, 1), width=c(0, 1, 0)))
replaceAt(x, at6, value="-")
value6 <- relist(paste0("[", seq_along(unlist(at6)), "]"), at6)
replaceAt(x, at6, value=value6)
replaceAt(x, at6, value=as(c("+", "-", "*"), "List"))

## Append "**" to all sequences:
replaceAt(x, as(width(x) + 1L, "List"), value="**")

## ---------------------------------------------------------------------
## (C) ADVANCED EXAMPLES
## ---------------------------------------------------------------------
library(hgu95av2probe)
probes <- DNAStringSet(hgu95av2probe)

## Split the probes in 5-mer chunks:
at <- successiveIRanges(rep(5, 5))
extractAt(probes, at)

## Replace base 13 by its complement:
at <- IRanges(13, width=1)
base13 <- extractAt(probes, at)
base13comp <- relist(complement(unlist(base13)), base13)
replaceAt(probes, at, value=base13comp)
## See ?xscat for a more efficient way to do this.

## Replace all the occurences of a given pattern with another pattern:
midx <- vmatchPattern("VCGTT", probes, fixed=FALSE)
matches <- extractAt(probes, midx)
unlist(matches)
unique(unlist(matches))
probes2 <- replaceAt(probes, midx, value="-++-")

## See strings with 2 or more susbtitutions:
probes2[elementLengths(midx) >= 2]

## 2 sanity checks:
stopifnot(all(replaceAt(probes, midx, value=matches) == probes))
probes2b <- gsub("[ACG]CGTT", "-++-", as.character(probes))
stopifnot(identical(as.character(probes2), probes2b))
}

\keyword{methods}
\keyword{manip}
