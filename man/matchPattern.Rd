\name{matchPattern}

\alias{matchPattern}
\alias{matchPattern,character-method}
\alias{matchPattern,XString-method}
\alias{matchPattern,XStringSet-method}
\alias{matchPattern,XStringViews-method}
\alias{matchPattern,MaskedXString-method}
\alias{countPattern}
\alias{countPattern,character-method}
\alias{countPattern,XString-method}
\alias{countPattern,XStringSet-method}
\alias{countPattern,XStringViews-method}
\alias{countPattern,MaskedXString-method}

\alias{vmatchPattern}
\alias{vmatchPattern,character-method}
\alias{vmatchPattern,XString-method}
\alias{vmatchPattern,XStringSet-method}
\alias{vmatchPattern,XStringViews-method}
\alias{vmatchPattern,MaskedXString-method}
\alias{vcountPattern}
\alias{vcountPattern,character-method}
\alias{vcountPattern,XString-method}
\alias{vcountPattern,XStringSet-method}
\alias{vcountPattern,XStringViews-method}
\alias{vcountPattern,MaskedXString-method}


\title{String searching functions}

\description{
  A set of functions for finding all the occurrences (aka "matches" or "hits")
  of a given pattern (typically short) in a (typically long) reference sequence
  or set of reference sequences (aka the subject)
}

\usage{
matchPattern(pattern, subject,
             max.mismatch=0, min.mismatch=0,
             with.indels=FALSE, fixed=TRUE,
             algorithm="auto")

countPattern(pattern, subject,
             max.mismatch=0, min.mismatch=0,
             with.indels=FALSE, fixed=TRUE,
             algorithm="auto")

vmatchPattern(pattern, subject,
              max.mismatch=0, min.mismatch=0,
              with.indels=FALSE, fixed=TRUE,
              algorithm="auto", ...)

vcountPattern(pattern, subject,
              max.mismatch=0, min.mismatch=0,
              with.indels=FALSE, fixed=TRUE,
              algorithm="auto", ...)
}

\arguments{
  \item{pattern}{
    The pattern string.
  }
  \item{subject}{
    An \link{XString}, \link{XStringViews} or \link{MaskedXString}
    object for \code{matchPattern} and \code{countPattern}.

    An \link{XStringSet} or \link{XStringViews} object for
    \code{vmatchPattern} and \code{vcountPattern}.
  }
  \item{max.mismatch, min.mismatch}{
    The maximum and minimum number of mismatching letters allowed (see
    \code{?`\link{lowlevel-matching}`} for the details).
    If non-zero, an algorithm that supports inexact matching is used.
  }
  \item{with.indels}{
    If \code{TRUE} then indels are allowed. In that case, \code{min.mismatch}
    must be \code{0} and \code{max.mismatch} is interpreted as the maximum
    "edit distance" allowed between the pattern and a match.
    Note that in order to avoid pollution by redundant matches,
    only the "best local matches" are returned.
    Roughly speaking, a "best local match" is a match that is locally
    both the closest (to the pattern P) and the shortest.
    More precisely, a substring S' of the subject S is a "best local match" iff:
    \preformatted{
       (a) nedit(P, S') <= max.mismatch
       (b) for every substring S1 of S':
               nedit(P, S1) > nedit(P, S')
       (c) for every substring S2 of S that contains S':
               nedit(P, S2) >= nedit(P, S')
    }
    One nice property of "best local matches" is that their first and last
    letters are guaranteed to match the letters in P that they align with.
  }
  \item{fixed}{
    If \code{TRUE} (the default), an IUPAC ambiguity code in the pattern
    can only match the same code in the subject, and vice versa.
    If \code{FALSE}, an IUPAC ambiguity code in the pattern can match
    any letter in the subject that is associated with the code, and 
    vice versa. See \code{?`\link{lowlevel-matching}`} for more information.
  }
  \item{algorithm}{
    One of the following: \code{"auto"}, \code{"naive-exact"},
    \code{"naive-inexact"}, \code{"boyer-moore"}, \code{"shift-or"}
    or \code{"indels"}.
  }
  \item{...}{
    Additional arguments for methods.
  }
}

\details{
  Available algorithms are: ``naive exact'', ``naive inexact'',
  ``Boyer-Moore-like'', ``shift-or'' and ``indels''.
  Not all of them can be used in all situations: restrictions
  apply depending on the "search criteria" i.e. on the values of
  the \code{pattern}, \code{subject}, \code{max.mismatch},
  \code{min.mismatch}, \code{with.indels} and \code{fixed}
  arguments.

  It is important to note that the \code{algorithm} argument
  is not part of the search criteria. This is because the supported
  algorithms are interchangeable, that is, if 2 different algorithms
  are compatible with a given search criteria, then choosing one or
  the other will not affect the result (but will most likely affect
  the performance). So there is no "wrong choice" of algorithm (strictly
  speaking).

  Using \code{algorithm="auto"} (the default) is recommended because
  then the best suited algorithm will automatically be selected among
  the set of algorithms that are valid for the given search criteria.
}

\value{
  An \link{XStringViews} object for \code{matchPattern}.

  A single integer for \code{countPattern}.

  An \link{MIndex} object for \code{vmatchPattern}.

  An integer vector for \code{vcountPattern}, with each element in
  the vector corresponding to the number of matches in the corresponding
  element of \code{subject}.
}

\note{
  Use \code{\link{matchPDict}} if you need to match a (big) set of patterns
  against a reference sequence.

  Use \code{\link{pairwiseAlignment}} if you need to solve a (Needleman-Wunsch)
  global alignment, a (Smith-Waterman) local alignment, or an (ends-free)
  overlap alignment problem.
}

\seealso{
  \link{lowlevel-matching},
  \code{\link{matchPDict}},
  \code{\link{pairwiseAlignment}},
  \code{\link{mismatch}},
  \code{\link{matchLRPatterns}},
  \code{\link{matchProbePair}},
  \code{\link{maskMotif}},
  \code{\link{alphabetFrequency}},
  \link{XStringViews-class},
  \link{MIndex-class}
}

\examples{
## ---------------------------------------------------------------------
## A. matchPattern()/countPattern()
## ---------------------------------------------------------------------

## A simple inexact matching example with a short subject:
x <- DNAString("AAGCGCGATATG")
m1 <- matchPattern("GCNNNAT", x)
m1
m2 <- matchPattern("GCNNNAT", x, fixed=FALSE)
m2
as.matrix(m2)

## With DNA sequence of yeast chromosome number 1:
data(yeastSEQCHR1)
yeast1 <- DNAString(yeastSEQCHR1)
PpiI <- "GAACNNNNNCTC" # a restriction enzyme pattern
match1.PpiI <- matchPattern(PpiI, yeast1, fixed=FALSE)
match2.PpiI <- matchPattern(PpiI, yeast1, max.mismatch=1, fixed=FALSE)

## With a genome containing isolated Ns:
library(BSgenome.Celegans.UCSC.ce2)
chrII <- Celegans[["chrII"]]
alphabetFrequency(chrII)
matchPattern("N", chrII)
matchPattern("TGGGTGTCTTT", chrII) # no match
matchPattern("TGGGTGTCTTT", chrII, fixed=FALSE) # 1 match

## Using wildcards ("N") in the pattern on a genome containing N-blocks:
library(BSgenome.Dmelanogaster.UCSC.dm3)
chrX <- maskMotif(Dmelanogaster$chrX, "N")
as(chrX, "Views") # 4 non masked regions
matchPattern("TTTATGNTTGGTA", chrX, fixed=FALSE)
## Can also be achieved with no mask:
masks(chrX) <- NULL
matchPattern("TTTATGNTTGGTA", chrX, fixed="subject")

## ---------------------------------------------------------------------
## B. vmatchPattern()/vcountPattern()
## ---------------------------------------------------------------------

Ebox <- DNAString("CANNTG")
subject <- Celegans$upstream5000
mindex <- vmatchPattern(Ebox, subject, fixed=FALSE)
count_index <- countIndex(mindex)  # Get the number of matches per
                                   # subject element.
sum(count_index)  # Total number of matches.
table(count_index)
i0 <- which(count_index == max(count_index))
subject[i0]  # The subject element with most matches.

## The matches in 'subject[i0]' as an IRanges object:
mindex[[i0]]
## The matches in 'subject[i0]' as an XStringViews object:
Views(subject[[i0]], mindex[[i0]])

## ---------------------------------------------------------------------
## C. WITH INDELS
## ---------------------------------------------------------------------

library(BSgenome.Celegans.UCSC.ce2)
subject <- Celegans$chrI
pattern1 <- DNAString("ACGGACCTAATGTTATC")
pattern2 <- DNAString("ACGGACCTVATGTTRTC")

## Allowing up to 2 mismatching letters doesn't give any match:
m1a <- matchPattern(pattern1, subject, max.mismatch=2)

## But allowing up to 2 edit operations gives 3 matches:
system.time(m1b <- matchPattern(pattern1, subject, max.mismatch=2,
                                with.indels=TRUE))
m1b

## pairwiseAlignment() returns the (first) best match only:
if (interactive()) {
  mat <- nucleotideSubstitutionMatrix(match=1, mismatch=0, baseOnly=TRUE)
  ## Note that this call to pairwiseAlignment() will need to
  ## allocate 733.5 Mb of memory (i.e. length(pattern) * length(subject)
  ## * 3 bytes).
  system.time(pwa <- pairwiseAlignment(pattern1, subject, type="local",
                                       substitutionMatrix=mat,
                                       gapOpening=0, gapExtension=1))
  pwa
}

## With IUPAC ambiguities in the pattern:
m2a <- matchPattern(pattern2, subject, max.mismatch=2,
                    fixed="subject")
m2b <- matchPattern(pattern2, subject, max.mismatch=2,
                    with.indels=TRUE, fixed="subject")

## All the matches in 'm1b' and 'm2a' should also appear in 'm2b':
stopifnot(suppressWarnings(all(ranges(m1b) %in% ranges(m2b))))
stopifnot(suppressWarnings(all(ranges(m2a) %in% ranges(m2b))))

## ---------------------------------------------------------------------
## D. WHEN 'with.indels=TRUE', ONLY "BEST LOCAL MATCHES" ARE REPORTED
## ---------------------------------------------------------------------

## With deletions in the subject:
subject <- BString("ACDEFxxxCDEFxxxABCE")
matchPattern("ABCDEF", subject, max.mismatch=2, with.indels=TRUE)
matchPattern("ABCDEF", subject, max.mismatch=2)

## With insertions in the subject:
subject <- BString("AiBCDiEFxxxABCDiiFxxxAiBCDEFxxxABCiDEF")
matchPattern("ABCDEF", subject, max.mismatch=2, with.indels=TRUE)
matchPattern("ABCDEF", subject, max.mismatch=2)

## With substitutions (note that the "best local matches" can introduce
## indels and therefore be shorter than 6):
subject <- BString("AsCDEFxxxABDCEFxxxBACDEFxxxABCEDF")
matchPattern("ABCDEF", subject, max.mismatch=2, with.indels=TRUE)
matchPattern("ABCDEF", subject, max.mismatch=2)
}

\keyword{methods}
