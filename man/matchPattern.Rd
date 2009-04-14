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

\alias{matchDNAPattern}


\title{String searching functions}

\description{
  A set of functions for finding all the occurences (aka "matches" or "hits")
  of a given pattern (typically short) in a (typically long) reference sequence
  or set of reference sequences (aka the subject)
}

\usage{
  matchPattern(pattern, subject, algorithm="auto",
               max.mismatch=0, with.indels=FALSE, fixed=TRUE)
  countPattern(pattern, subject, algorithm="auto",
               max.mismatch=0, with.indels=FALSE, fixed=TRUE)
  vmatchPattern(pattern, subject, algorithm="auto",
                max.mismatch=0, with.indels=FALSE, fixed=TRUE)
  vcountPattern(pattern, subject, algorithm="auto",
                max.mismatch=0, with.indels=FALSE, fixed=TRUE)
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
  \item{algorithm}{
    One of the following: \code{"auto"}, \code{"naive-exact"},
    \code{"naive-inexact"}, \code{"boyer-moore"}, \code{"shift-or"}
    or \code{"indels"}.
  }
  \item{max.mismatch}{
    The maximum number of mismatching letters allowed (see
    \code{\link{isMatchingAt}} for the details).
    If non-zero, an inexact matching algorithm is used.
  }
  \item{with.indels}{
    If \code{TRUE} then indels are allowed. In that case \code{max.mismatch}
    is interpreted as the maximum "edit distance" allowed between the pattern
    and a match. Note that in order to avoid pollution by redundant matches,
    only the "best local matches" are returned.
    Roughly speaking, a "best local match" is a match that is locally
    both the closest (to the pattern P) and the shortest.
    More precisely, a substring S' of the subject S is a "best local match" iff:
    \preformatted{
       (a) nedit(P, S') <= max.mismatch
       (b) for every substring S1 of S':
               nedit(P, S1) > nedit(P, S')
       (c) for every substring S2 of S that contains S':
               nedit(P, S2) <= nedit(P, S')
    }
    One nice property of "best local matches" is that their first and last
    letters are guaranteed to be aligned with letters in P (i.e. they match
    letters in P).
  }
  \item{fixed}{
    If \code{FALSE} then IUPAC extended letters are interpreted as ambiguities
    (see \code{\link{isMatchingAt}} for the details).
  }
}

\details{
  Available algorithms are: ``naive exact'', ``naive inexact'',
  ``Boyer-Moore-like'', ``shift-or'' and ``indels''.
  Not all of them can be used in all situations: restrictions
  depend on the length of the pattern, the class of the subject,
  and the values of \code{max.mismatch}, \code{with.indels}
  and \code{fixed}. All those parameters form the search criteria.

  Note that the choice of an algorithm is not part of the search criteria.
  This is because algorithms are interchangeable, that is, if 2 different
  algorithms are compatible with a given search criteria, then choosing
  one over the other will not affect the result (but will most likely
  affect the performance). So there is no "wrong choice" of algorithm
  (strictly speaking).

  Using \code{algorithm="auto"} is recommended because then the fastest
  algorithm will automatically be picked up among the set of compatible
  algorithms (if there is more than one).
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
  \code{\link{matchPDict}},
  \code{\link{pairwiseAlignment}},
  \code{\link{isMatchingAt}},
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
  as(chrX, "XStringViews") # 4 non masked regions
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
  ## C. With indels
  ## ---------------------------------------------------------------------
  library(BSgenome.Celegans.UCSC.ce2)
  pattern <- DNAString("ACGGACCTAATGTTATC")
  subject <- Celegans$chrI

  ## Allowing up to 2 mismatching letters doesn't give any match:
  matchPattern(pattern, subject, max.mismatch=2)

  ## But allowing up to 2 edit operations gives 3 matches:
  system.time(m <- matchPattern(pattern, subject, max.mismatch=2, with.indels=TRUE))
  m

  ## pairwiseAlignment() returns the (first) best match only:
  mat <- nucleotideSubstitutionMatrix(match=1, mismatch=0, baseOnly=TRUE)
  system.time(pwa <- pairwiseAlignment(pattern, subject, type="local",
                       substitutionMatrix=mat, gapOpening=0, gapExtension=1))
  pwa

  ## Only "best local matches" are reported:
    ## - with deletions in the subject
  subject <- BString("ACDEFxxxCDEFxxxABCE")
  matchPattern("ABCDEF", subject, max.mismatch=2, with.indels=TRUE)
  matchPattern("ABCDEF", subject, max.mismatch=2)
    ## - with insertions in the subject
  subject <- BString("AiBCDiEFxxxABCDiiFxxxAiBCDEFxxxABCiDEF")
  matchPattern("ABCDEF", subject, max.mismatch=2, with.indels=TRUE)
  matchPattern("ABCDEF", subject, max.mismatch=2)
    ## - with substitutions (note that the "best local matches" can introduce
    ##   indels and therefore be shorter than 6)
  subject <- BString("AsCDEFxxxABDCEFxxxBACDEFxxxABCEDF")
  matchPattern("ABCDEF", subject, max.mismatch=2, with.indels=TRUE)
  matchPattern("ABCDEF", subject, max.mismatch=2)
}

\keyword{methods}
